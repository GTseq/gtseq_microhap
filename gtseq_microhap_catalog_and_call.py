#!/usr/bin/env python3
"""
gtseq_microhap_catalog_and_call.py

Alignment-free microhaplotype genotyping for GT-seq amplicon sequencing data.

This script implements a GT-seq amplicon analysis pipeline that resolves
paired-end reads into primer-bounded amplicons, builds a microhaplotype allele
catalog, calls diploid genotypes, exports phased microhaplotype and optional
single-SNP outputs, and supports accessory presence/absence assays such as
sexID and parasite detection.

Key features in version 0.3.1
-----------------------------------------------
1. Quality-aware R1/R2 read resolution
   - Identifies R1 reads beginning with a locus-specific forward primer
   - Confirms the reverse primer at the beginning of R2
   - Resolves standard overlapping pairs using quality-aware overlap consensus
   - Handles short read-through amplicons where both oriented reads span the
     full primer-bounded amplicon
   - Chooses the higher quality base call when R1 and R2 disagree

2. Allele discovery and catalog construction
   - Ranks unique amplicon sequences within each sample and locus
   - Retains up to two candidate alleles per sample under a diploid model
   - Aggregates supported alleles across samples into a catalog
   - Filters catalog outliers using dominant length class, length difference,
     fractional length difference, and pairwise identity to the top allele

3. Genotype inference and microhaplotype export
   - Performs a catalog-only second pass through resolved reads
   - Calls HOM, HET, NC, and LOW genotypes from A1/A2 abundance patterns
   - Exports allele-coded genotypes, phased SNP/variant strings, and optional
     single-SNP matrices

4. Accessory assay support
   - Optional sexID marker detection using p_OT-normalized read abundance
   - Optional parasite/presence-absence marker detection using normalized
     marker abundance categories
   - Accessory marker dashboards and summary tables

Outputs (written to --outdir)
-----------------------------
sample_stats.tsv
microhap_catalog.csv
microhap_catalog_length_pruned.tsv
microhap_consensus.fasta
microhap_a2_metrics.tsv
microhap_genotypes.tsv
microhap_phased_def.tsv
microhap_phased.tsv
microhap_phasedVARs.tsv
microhap_singleSNP_first.tsv / microhap_singleSNP_highestMAF.tsv (optional)
microhap_singleSNP_metadata.tsv (optional)
sexid_calls.tsv / sexid_marker_pot.tsv / sexid_sample_marker_scores.tsv (optional)
parasite_calls.tsv / parasite_summary.tsv / parasite_sample_marker_scores.tsv (optional)
plots/  (optional dashboards and library summaries)

Repository
----------
https://github.com/GTseq/gtseq_microhap

Citation
--------
If you use this software, please cite the associated GT-seq microhaplotype
manuscript/preprint and the original GT-seq method paper.

Version
-------
0.3.1

Author
------
Nathan Campbell
GTseek LLC

"""

__version__ = "0.3.1"

import os, sys, re, gzip, csv, collections, math
import statistics
from math import sqrt
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---- Optional progress bars (tqdm if available; graceful fallback) ----
try:
    from tqdm import tqdm as _tqdm
    def PROG(iterable, desc="", unit="", total=None):
        return _tqdm(iterable, desc=desc, unit=unit, total=total)
except Exception:
    def PROG(iterable, desc="", unit="", total=None):
        print(f"{desc}...", flush=True)
        count = 0
        for x in iterable:
            count += 1
            if count == 1 or (count % 25 == 0) or (total and count == total):
                if total:
                    print(f"{desc}: {count}/{total} {unit}", flush=True)
                else:
                    print(f"{desc}: {count} {unit}", flush=True)
            yield x

# headless plotting -----------------------------------------------------------
os.environ.setdefault("MPLBACKEND", "Agg")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAVE_MPL = True
except Exception as e:
    HAVE_MPL = False
    print(f"Warning: matplotlib import failed: {e}", file=sys.stderr)

# ---------------------------------------------------------------------
# DNA utils
# ---------------------------------------------------------------------
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]

def open_read(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def next_fastq_record(fh):
    h = fh.readline()
    if not h:
        return None
    s = fh.readline()
    p = fh.readline()
    q = fh.readline()
    if not q:
        return None
    return (h.rstrip(), s.rstrip(), p.rstrip(), q.rstrip())


def count_fastq_records(path: Path) -> int:
    """Count FASTQ records in a .fastq or .fastq.gz file."""
    n = 0
    with open_read(path) as fh:
        while True:
            rec = next_fastq_record(fh)
            if rec is None:
                break
            n += 1
    return n

# ---------------------------------------------------------------------
# primers
# ---------------------------------------------------------------------
def load_primers(path: Path) -> Dict[str, Tuple[str, str]]:
    primers = {}
    with open(path) as fh:
        header = fh.readline()
        if not header:
            sys.exit("ERROR: empty primer file.")
        delim = "," if header.count(",") > header.count("\t") else "\t"
        cols = [c.strip() for c in header.strip().split(delim)]
        name_map = {"locus": None, "fwd_primer": None, "rev_primer": None}
        for i, c in enumerate(cols):
            lc = c.lower()
            if lc == "locus":
                name_map["locus"] = i
            elif lc in ("fwd_primer","forward_primer","fprimer","fwd"):
                name_map["fwd_primer"] = i
            elif lc in ("rev_primer","reverse_primer","rprimer","rev"):
                name_map["rev_primer"] = i
        if any(v is None for v in name_map.values()):
            sys.exit("ERROR: primer header must contain columns for locus,fwd_primer,rev_primer (or equivalent names).")

        for line in fh:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.rstrip().split(delim)]
            loc = parts[name_map["locus"]]
            fwd = parts[name_map["fwd_primer"]].upper()
            rev = parts[name_map["rev_primer"]].upper()
            if loc and fwd and rev:
                primers[loc] = (fwd, rev)

    if not primers:
        sys.exit("ERROR: no primers loaded from file.")
    return primers

# ---------------------------------------------------------------------
# sample pairing
# ---------------------------------------------------------------------
R1_TAG = re.compile(r"(^|[_\-\.])R1([_\-\.]|$)")
R2_TAG = re.compile(r"(^|[_\-\.])R2([_\-\.]|$)")

def infer_sample_base(fn: str) -> str:
    base = re.sub(r"(_S\d+)?_L00\d(_R[12])?_001\.fastq(?:\.gz)?$", "", fn)
    base = re.sub(r"_R[12]\.fastq(?:\.gz)?$", "", base)
    base = re.sub(r"\.fastq(?:\.gz)?$", "", base)
    base = re.sub(r"([_\-\.])R[12].*$", "", base)
    return base

def find_fastq_pairs(indir: Path) -> Dict[str, Tuple[Path, Path]]:
    files = [p for p in indir.iterdir()
             if p.name.endswith(".fastq") or p.name.endswith(".fastq.gz")]
    R1s, R2s = {}, {}
    for p in files:
        if R1_TAG.search(p.name):
            base = infer_sample_base(p.name); R1s.setdefault(base, []).append(p)
        elif R2_TAG.search(p.name):
            base = infer_sample_base(p.name); R2s.setdefault(base, []).append(p)
    pairs = {}
    for base in sorted(set(R1s) & set(R2s)):
        r1 = sorted(R1s[base])[0]; r2 = sorted(R2s[base])[0]
        pairs[base] = (r1, r2)
    return pairs

# ---------------------------------------------------------------------
# trim / merge
# ---------------------------------------------------------------------


def trim_oriented_to_amplicon(seq: str, qual: str, fwd: str, rev: str, require_fwd_at_start: bool = False) -> Optional[Tuple[str, str]]:
    """
    Trim an oriented read to the primer-bounded amplicon (FWD ... REV_RC).

    For R1 we normally require the read to begin with FWD.
    For RC(R2), short amplicons can introduce read-through bases ahead of FWD,
    so we allow FWD to occur internally.
    """
    s = seq.upper()
    q = qual.rstrip()
    if len(q) != len(s):
        return None

    fwd = fwd.upper()
    rev_rc = rc(rev)

    if require_fwd_at_start:
        if not s.startswith(fwd):
            return None
        i = 0
    else:
        i = s.find(fwd)
        if i == -1:
            return None

    j = s.find(rev_rc, i + len(fwd))
    if j == -1:
        return None
    end = j + len(rev_rc)
    return (s[i:end], q[i:end])


def consensus_equal_length(a_seq: str, a_qual: str, b_seq: str, b_qual: str) -> Optional[Tuple[str, str]]:
    """Build a quality-aware consensus for two equal-length amplicons."""
    a_seq = a_seq.upper()
    b_seq = b_seq.upper()
    qa = phred_scores(a_qual)
    qb = phred_scores(b_qual)

    if len(a_seq) != len(b_seq) or len(qa) != len(a_seq) or len(qb) != len(b_seq):
        return None

    merged_seq = []
    merged_q = []
    for x, qx, y, qy in zip(a_seq, qa, b_seq, qb):
        if x == y:
            base = x
            q = max(qx, qy)
        elif x == 'N' and y != 'N':
            base = y
            q = qy
        elif y == 'N' and x != 'N':
            base = x
            q = qx
        else:
            if qx >= qy:
                base = x
                q = qx
            else:
                base = y
                q = qy
        merged_seq.append(base)
        merged_q.append(q)
    return (''.join(merged_seq), scores_to_qual(merged_q))

def r2_starts_with_rev(seq: str, rev: str) -> bool:
    return seq.upper().startswith(rev.upper())

def phred_scores(q: str) -> List[int]:
    return [max(0, ord(c) - 33) for c in q.rstrip()]

def scores_to_qual(scores: List[int]) -> str:
    return "".join(chr(max(0, min(93, s)) + 33) for s in scores)

def rev_qual(q: str) -> str:
    return q.rstrip()[::-1]

def overlap_merge(a: str, b: str, qa: str, qb: str, min_ov: int, max_mismatch_frac: float) -> Optional[Tuple[str, str]]:
    """
    Quality-aware suffix/prefix overlap merge.
    `a` and `b` must already be oriented in the same direction.
    Returns (merged_sequence, merged_quality) or None.
    """
    a = a.upper()
    b = b.upper()
    qa_scores = phred_scores(qa)
    qb_scores = phred_scores(qb)

    if len(qa_scores) != len(a) or len(qb_scores) != len(b):
        return None

    max_len = min(len(a), len(b))

    for ov in range(max_len, min_ov - 1, -1):
        aa = a[-ov:]
        bb = b[:ov]

        mism = 0
        limit = max(1, int(max_mismatch_frac * ov))
        for x, y in zip(aa, bb):
            if x == y:
                continue
            if x == "N" or y == "N":
                continue
            mism += 1
            if mism > limit:
                break
        else:
            merged_seq = list(a[:-ov])
            merged_q = qa_scores[:-ov]

            a_q_ov = qa_scores[-ov:]
            b_q_ov = qb_scores[:ov]

            for x, qx, y, qy in zip(aa, a_q_ov, bb, b_q_ov):
                if x == y:
                    base = x
                    q = max(qx, qy)
                elif x == "N" and y != "N":
                    base = y
                    q = qy
                elif y == "N" and x != "N":
                    base = x
                    q = qx
                else:
                    if qx >= qy:
                        base = x
                        q = qx
                    else:
                        base = y
                        q = qy

                merged_seq.append(base)
                merged_q.append(q)

            merged_seq.extend(list(b[ov:]))
            merged_q.extend(qb_scores[ov:])

            return ("".join(merged_seq), scores_to_qual(merged_q))

    return None

def build_merged_amplicon(
    r1_seq: str,
    r1_qual: str,
    r2_seq: str,
    r2_qual: str,
    fwd: str,
    rev: str,
    min_ov: int,
    max_mismatch_frac: float,
) -> Optional[Tuple[str, str]]:
    """
    Return (amplicon_sequence, amplicon_quality).

    We always use both reads, but there are two valid paired-end geometries:

    1) Standard partial-overlap pairs, where R1 and RC(R2) must be stitched.
    2) Short read-through amplicons, where *both* oriented reads independently span
       the full primer-bounded amplicon. In that case we should not force a raw
       suffix/prefix overlap on the untrimmed reads because RC(R2) may contain
       leading read-through sequence ahead of FWD. Instead, trim both oriented reads
       to FWD..REV_RC and build a quality-aware full-length consensus.
    """
    if not r1_seq.upper().startswith(fwd.upper()):
        return None
    if not r2_starts_with_rev(r2_seq, rev):
        return None

    a = r1_seq.upper()
    qa = r1_qual.rstrip()
    b = rc(r2_seq.upper())
    qb = rev_qual(r2_qual)

    # First handle short read-through amplicons by trimming each oriented read to the
    # primer-bounded amplicon and consensusing those directly.
    a_trim = trim_oriented_to_amplicon(a, qa, fwd, rev, require_fwd_at_start=True)
    b_trim = trim_oriented_to_amplicon(b, qb, fwd, rev, require_fwd_at_start=False)
    if a_trim is not None and b_trim is not None:
        ampa, qampa = a_trim
        ampb, qampb = b_trim
        merged = consensus_equal_length(ampa, qampa, ampb, qampb)
        if merged is not None:
            return merged
        # If both reads span the full amplicon but disagree on length, discard.
        return None

    # Otherwise, require a successful paired-end stitch from the full oriented reads.
    merged = overlap_merge(a, b, qa, qb, min_ov, max_mismatch_frac)
    if merged is None:
        return None

    merged_seq, merged_qual = merged
    trimmed = trim_oriented_to_amplicon(merged_seq, merged_qual, fwd, rev, require_fwd_at_start=True)
    if trimmed is None:
        return None
    return trimmed

def _select_top2_for_catalog(cnts: "collections.Counter[str]", hom_min=0.85, het_min=0.30):
    total = sum(cnts.values())
    if total == 0:
        return []
    top = cnts.most_common(2)
    s1, c1 = top[0]
    f1 = c1 / total
    if len(top) == 1:
        return [s1] if f1 >= hom_min else []

    s2, c2 = top[1]
    f2 = c2 / total

    if f1 >= hom_min and f2 < het_min:
        return [s1]
    if f1 >= het_min and f2 >= het_min:
        return [s1, s2]
    return []

def _n_workers(requested: int | None = None) -> int:
    if requested is not None and requested > 0:
        return requested
    try:
        n = os.cpu_count() or 1
    except Exception:
        n = 1
    return max(1, n // 2)

# --- TOP-LEVEL worker for ProcessPoolExecutor (must not be nested) ---
def _worker_resolve(args):
    (r1, r2, sample, primers, resolved_dir,
     min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac,
     force) = args
    return process_sample_pair(
        r1, r2, sample,
        primers=primers,
        min_amplicon=min_amplicon,
        max_amplicon=max_amplicon,
        min_overlap=min_overlap,
        max_ov_mismatch_frac=max_ov_mismatch_frac,
        out_dir=resolved_dir,
        force=force,
    )

def resolve_fastqs_multiprocess(
    paired_samples,
    resolved_dir: Path,
    primers: dict,
    min_amplicon: int = 60,
    max_amplicon: int = 250,
    min_overlap: int = 12,
    max_ov_mismatch_frac: float = 0.05,
    force: bool = False,
    workers: int | None = None,
    out_stats_path: Path | None = None,
):
    resolved_dir.mkdir(parents=True, exist_ok=True)
    if workers is None:
        workers = _n_workers(None)

    if not paired_samples:
        return 0, []

    jobs = [
        (r1, r2, sample, primers, resolved_dir,
         min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac,
         force)
        for (r1, r2, sample) in paired_samples
    ]

    errors: list[str] = []
    stats_by_sample: dict[str, dict] = {}

    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_worker_resolve, args): args for args in jobs}
        iterator = as_completed(futs)
        try:
            from tqdm import tqdm as _tqdm
            iterator = _tqdm(iterator, total=len(futs), desc="Resolving FASTQs (parallel)", unit="sample")
        except Exception:
            pass

        for fut in iterator:
            args = futs[fut]
            sample = args[2]
            try:
                res = fut.result()
                if isinstance(res, dict):
                    stats_by_sample[sample] = res
            except Exception as e:
                errors.append(f"{sample}: {e}")

    if out_stats_path is not None and stats_by_sample:
        with Path(out_stats_path).open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "raw_reads", "primer_bounded_reads"])
            for sample in sorted(stats_by_sample.keys()):
                st = stats_by_sample[sample]
                raw = int(st.get("read_pairs", 0))
                pb  = int(st.get("resolved",   0))
                w.writerow([sample, raw, pb])

    return len(stats_by_sample), errors

def process_sample_pair(
    r1_path,
    r2_path,
    sample_name,
    primers,
    min_amplicon=60,
    max_amplicon=250,
    min_overlap=12,
    max_ov_mismatch_frac=0.05,
    out_dir=Path("resolved_fastqs"),
    force=False,
    **_ignored,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"resolved_{sample_name}.fastq.gz"
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    if out_path.exists() and not force:
        raw_reads = count_fastq_records(Path(r1_path))
        primer_bounded_reads = count_fastq_records(out_path)
        return {
            "sample": sample_name,
            "read_pairs": raw_reads,
            "resolved": primer_bounded_reads,
            "skipped_existing": True,
        }

    fh1 = open_read(r1_path)
    fh2 = open_read(r2_path)
    out_fh = gzip.open(tmp_path, "wt")

    read_pairs = 0
    resolved = 0

    try:
        while True:
            rec1 = next_fastq_record(fh1)
            rec2 = next_fastq_record(fh2)
            if rec1 is None or rec2 is None:
                break
            read_pairs += 1

            h1, s1, _p1, q1 = rec1
            _h2, s2, _p2, q2 = rec2

            for locus, (fwd, rev) in primers.items():
                if not s1.upper().startswith(fwd.upper()):
                    continue

                built = build_merged_amplicon(
                    s1, q1, s2, q2, fwd, rev,
                    min_ov=min_overlap,
                    max_mismatch_frac=max_ov_mismatch_frac,
                )
                if not built:
                    continue

                amplicon, qual = built
                L = len(amplicon)
                if L < min_amplicon or L > max_amplicon:
                    break

                resolved += 1
                out_h = h1.split()[0] + f"|locus={locus}\n"

                out_fh.write(out_h)
                out_fh.write(amplicon + "\n")
                out_fh.write("+\n")
                out_fh.write(qual[:L] + "\n")
                break
            # next pair

    finally:
        fh1.close()
        fh2.close()
        out_fh.close()

    tmp_path.replace(out_path)
    return {"sample": sample_name, "read_pairs": read_pairs, "resolved": resolved, "skipped_existing": False}

# ---------------------------------------------------------------------
# locus readshare stats
# ---------------------------------------------------------------------
def compute_locus_readshare_stats(locus_depths: dict):
    mean_pct = {}
    sd_pct = {}
    for locus, vals in locus_depths.items():
        if not vals:
            mean_pct[locus] = 0.0
            sd_pct[locus] = 0.0
            continue
        m = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((v - m) ** 2 for v in vals) / (len(vals) - 1)
            sd = sqrt(var)
        else:
            sd = 0.0
        mean_pct[locus] = m
        sd_pct[locus] = sd
    sorted_loci_by_mean = sorted(mean_pct.keys(), key=lambda L: mean_pct[L])
    return {"sorted_loci_by_mean": sorted_loci_by_mean, "mean_pct": mean_pct, "sd_pct": sd_pct}

def compute_locus_readshare_stats_from_resolved(resolved_dir: Path) -> dict:
    resolved_paths = sorted(p for p in resolved_dir.glob("resolved_*.fastq.gz") if not str(p).endswith(".part"))
    locus_depths: dict[str, list[float]] = defaultdict(list)

    for fpath in PROG(resolved_paths, desc="Scanning resolved FASTQs", unit="sample", total=len(resolved_paths)):
        per_locus = Counter()
        total = 0
        with gzip.open(fpath, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                fh.readline()
                _ = fh.readline()
                _ = fh.readline()
                if "|locus=" not in h:
                    continue
                loc = h.strip().split("|locus=")[-1]
                per_locus[loc] += 1
                total += 1

        if total > 0:
            for loc, ct in per_locus.items():
                locus_depths[loc].append(ct / total)

    return compute_locus_readshare_stats(locus_depths)



def _support_weighted_median_length(seq_counter: "collections.Counter[str]", lengths: set[int]) -> int:
    """Return the support-weighted median length among sequences whose lengths are in `lengths`."""
    rows = sorted((len(seq), int(support)) for seq, support in seq_counter.items() if len(seq) in lengths)
    total = sum(support for _length, support in rows)
    if total <= 0 or not rows:
        return 0
    halfway = total / 2.0
    running = 0
    for length, support in rows:
        running += support
        if running >= halfway:
            return int(length)
    return int(rows[-1][0])


def _dominant_length_cluster(seq_counter: "collections.Counter[str]", bin_bp: int = 2) -> tuple[int, set[int], int]:
    """Identify the support-weighted dominant allele-length cluster for one locus.

    Lengths within +/- `bin_bp` are treated as one small-indel/rounding cluster.
    Returns (main_length, accepted_cluster_lengths, cluster_support).
    """
    if not seq_counter:
        return (0, set(), 0)

    bin_bp = max(0, int(bin_bp))
    length_support = Counter()
    for seq, support in seq_counter.items():
        length_support[len(seq)] += int(support)

    best_center = None
    best_lengths: set[int] = set()
    best_support = -1

    for center in sorted(length_support):
        cluster_lengths = {L for L in length_support if abs(L - center) <= bin_bp}
        cluster_support = sum(length_support[L] for L in cluster_lengths)
        if (
            cluster_support > best_support
            or (
                cluster_support == best_support
                and best_center is not None
                and length_support[center] > length_support[best_center]
            )
            or (
                cluster_support == best_support
                and best_center is not None
                and length_support[center] == length_support[best_center]
                and center < best_center
            )
        ):
            best_center = center
            best_lengths = cluster_lengths
            best_support = cluster_support

    main_len = _support_weighted_median_length(seq_counter, best_lengths)
    if main_len <= 0:
        main_len = int(best_center or 0)
    return (main_len, best_lengths, int(best_support))


def _catalog_pairwise_identity(ref_seq: str, alt_seq: str) -> float:
    """Global-alignment identity for short amplicon catalog candidates.

    Identity is calculated over non-double-gap alignment columns using the same
    simple Needleman-Wunsch routine used elsewhere in this script. This catches
    primer-bounded off-locus amplicons that are similar in length but clearly not
    sequence-compatible with the dominant allele cluster.
    """
    ref_seq = (ref_seq or "").upper()
    alt_seq = (alt_seq or "").upper()
    if not ref_seq or not alt_seq:
        return 0.0
    try:
        a_ref, a_alt = _nw_global_align(ref_seq, alt_seq)
    except Exception:
        # Extremely conservative fallback: position-wise identity over the longer length.
        denom = max(len(ref_seq), len(alt_seq), 1)
        matches = sum(1 for x, y in zip(ref_seq, alt_seq) if x == y)
        return matches / float(denom)

    compared = 0
    matches = 0
    for x, y in zip(a_ref, a_alt):
        if x == "-" and y == "-":
            continue
        compared += 1
        if x == y and x != "-":
            matches += 1
    return (matches / float(compared)) if compared else 0.0


def _filter_catalog_length_outliers(
    locus_seq_support: dict,
    max_len_diff: int = 20,
    max_len_ratio: float = 0.20,
    length_bin_bp: int = 2,
    min_pid_to_top: float = 0.85,
    enabled: bool = True,
):
    """Prune implausible allele outliers from aggregated catalog candidates.

    The first safeguard removes large length outliers relative to the dominant
    support-weighted length class. The second safeguard removes sequence-divergent
    off-locus amplicons that happen to be close enough in length to pass the length
    filter. A sequence is retained only if it passes BOTH checks:

      1. length difference <= max_len_diff OR fractional difference <= max_len_ratio
      2. global identity to the top supported allele >= min_pid_to_top

    This still allows small indels and normal SNP variation, but keeps stray
    primer-bounded amplicons out of the catalog before phased SNP construction.
    """
    kept_by_locus = defaultdict(Counter)
    pruned_rows = []

    if not enabled:
        for locus, seq_counter in locus_seq_support.items():
            kept_by_locus[locus].update(seq_counter)
        return kept_by_locus, pruned_rows

    max_len_diff = max(0, int(max_len_diff))
    max_len_ratio = max(0.0, float(max_len_ratio))
    min_pid_to_top = max(0.0, min(1.0, float(min_pid_to_top)))
    length_bin_bp = max(0, int(length_bin_bp))

    for locus, seq_counter in locus_seq_support.items():
        if not seq_counter:
            continue

        main_len, _cluster_lengths, _cluster_support = _dominant_length_cluster(seq_counter, bin_bp=length_bin_bp)
        if main_len <= 0:
            kept_by_locus[locus].update(seq_counter)
            continue

        top_seq, _top_support = max(seq_counter.items(), key=lambda kv: kv[1])

        for seq, support in seq_counter.items():
            L = len(seq)
            diff = abs(L - main_len)
            ratio = diff / float(main_len) if main_len else 0.0
            len_ok = (diff <= max_len_diff) or (ratio <= max_len_ratio)
            pid = 1.0 if seq == top_seq else _catalog_pairwise_identity(top_seq, seq)
            pid_ok = (pid >= min_pid_to_top)
            keep = len_ok and pid_ok
            if keep:
                kept_by_locus[locus][seq] += support
            else:
                if not len_ok and not pid_ok:
                    reason = "catalog_length_and_divergence_outlier"
                elif not len_ok:
                    reason = "catalog_length_outlier"
                else:
                    reason = "catalog_sequence_divergence_outlier"
                pruned_rows.append({
                    "locus": locus,
                    "sequence": seq,
                    "support_count": int(support),
                    "length": int(L),
                    "main_length": int(main_len),
                    "len_diff": int(diff),
                    "len_ratio": round(ratio, 6),
                    "pid_to_top": round(pid, 6),
                    "top_sequence": top_seq,
                    "reason": reason,
                })

    return kept_by_locus, pruned_rows

def build_catalog_from_resolved(
    resolved_dir,
    out_catalog_csv,
    out_consensus_fa,
    hom_min=0.85,
    het_min=0.30,
    min_catalog_depth=10,
    catalog_max_len_diff=20,
    catalog_max_len_ratio=0.20,
    catalog_length_bin_bp=2,
    catalog_min_pid_to_top=0.85,
    disable_catalog_length_filter=False,
):
    resolved_dir = Path(resolved_dir)

    locus_seq_support = defaultdict(Counter)
    locus_total_support = Counter()

    fastqs = sorted(p for p in resolved_dir.glob("resolved_*.fastq.gz") if not str(p).endswith(".part"))

    for fq in PROG(fastqs, desc="Building catalog from resolved FASTQs", unit="sample", total=len(fastqs)):
        per_locus_counts = defaultdict(Counter)

        with gzip.open(fq, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().strip()
                _ = fh.readline()
                _ = fh.readline()

                if "|locus=" not in h:
                    continue
                locus = h.strip().split("|locus=")[-1]
                per_locus_counts[locus][s] += 1

        for locus, seq_counts in per_locus_counts.items():
            total_depth = sum(seq_counts.values())
            if total_depth < min_catalog_depth:
                continue

            keep_seqs = _select_top2_for_catalog(seq_counts, hom_min=hom_min, het_min=het_min)
            if not keep_seqs:
                continue

            for seq in keep_seqs:
                c = seq_counts[seq]
                locus_seq_support[locus][seq] += c
                locus_total_support[locus] += c

    out_catalog_csv = Path(out_catalog_csv)
    out_consensus_fa = Path(out_consensus_fa)
    pruned_tsv = out_catalog_csv.with_name(out_catalog_csv.stem + "_length_pruned.tsv")

    filtered_locus_seq_support, pruned_rows = _filter_catalog_length_outliers(
        locus_seq_support,
        max_len_diff=catalog_max_len_diff,
        max_len_ratio=catalog_max_len_ratio,
        length_bin_bp=catalog_length_bin_bp,
        min_pid_to_top=catalog_min_pid_to_top,
        enabled=(not disable_catalog_length_filter),
    )

    if pruned_rows:
        with pruned_tsv.open("w", newline="") as out_pruned:
            w = csv.writer(out_pruned, delimiter="\t")
            w.writerow(["locus", "sequence", "support_count", "length", "main_length", "len_diff", "len_ratio", "pid_to_top", "top_sequence", "reason"])
            for r in sorted(pruned_rows, key=lambda x: (x["locus"], -x["support_count"], x["sequence"])):
                w.writerow([r["locus"], r["sequence"], r["support_count"], r["length"], r["main_length"], r["len_diff"], r["len_ratio"], r.get("pid_to_top", ""), r.get("top_sequence", ""), r["reason"]])
        print(f"Catalog length filter pruned {len(pruned_rows)} allele sequence(s); wrote {pruned_tsv}", file=sys.stderr)
    elif not disable_catalog_length_filter:
        with pruned_tsv.open("w", newline="") as out_pruned:
            w = csv.writer(out_pruned, delimiter="\t")
            w.writerow(["locus", "sequence", "support_count", "length", "main_length", "len_diff", "len_ratio", "pid_to_top", "top_sequence", "reason"])

    with out_catalog_csv.open("w") as out_csv, out_consensus_fa.open("w") as out_fa:
        out_csv.write("locus,allele_code,sequence,support_count,total_support\n")
        for locus in sorted(filtered_locus_seq_support.keys()):
            seq_counter = filtered_locus_seq_support[locus]
            if not seq_counter:
                continue
            sorted_seqs = sorted(seq_counter.items(), key=lambda kv: kv[1], reverse=True)
            total_support = sum(seq_counter.values())
            allele_code = 101
            consensus_seq = sorted_seqs[0][0]
            out_fa.write(f">{locus}\n{consensus_seq}\n")
            for seq, support in sorted_seqs:
                out_csv.write(f"{locus},{allele_code},{seq},{support},{total_support}\n")
                allele_code += 1

def _allele_freq_from_genos(per_sample_gt: dict[str,str]) -> dict[int, float]:
    counts = Counter()
    n_alleles = 0
    for gt in (per_sample_gt or {}).values():
        if not gt or gt == "000000" or len(gt) != 6:
            continue
        a = int(gt[:3]); b = int(gt[3:])
        counts[a] += 1; counts[b] += 1
        n_alleles += 2
    if n_alleles == 0:
        return {}
    freqs = {k: v / float(n_alleles) for k, v in counts.items()}
    s = sum(freqs.values())
    if s and abs(s - 1.0) > 1e-12:
        kmax = max(freqs, key=freqs.get)
        freqs[kmax] += (1.0 - s)
    return freqs

def call_from_resolved_fastqs(
    resolved_dir: Path,
    catalog_csv: Path,
    out_prefix: str,
    min_depth: int,
    locus_stats: dict,
    plots_dir: Path | None = None,
    a2_lo: float = 10.0,
    a2_hi: float = 25.0,
    locus_universe=None,
    assay_exclude_loci=None,
):
    """Second pass: call genotypes from resolved FASTQs against the catalog.

    Also emits optional phased hap strings derived from catalog allele sequences:
      - *_phased_def.tsv : variable-site definitions per locus
      - *_phased.tsv     : per-sample hap strings for each called locus
    """
    # Genotype decision thresholds are provided as A2 percentages on the CLI.
    # Defaults reproduce the original behavior: HOM if A2 <= 10%, HET if A2 >= 25%,
    # and NC in the intermediate zone.
    a2_lo_frac = float(a2_lo) / 100.0
    a2_hi_frac = float(a2_hi) / 100.0
    if not (0.0 <= a2_lo_frac < a2_hi_frac <= 1.0):
        raise ValueError("a2_lo and a2_hi must satisfy 0 <= a2_lo < a2_hi <= 100")

    # ---- load catalog: locus -> {seq -> code} ----
    catalog: dict[str, dict[str, int]] = defaultdict(dict)
    with Path(catalog_csv).open() as fh:
        rdr = csv.DictReader(fh)
        cols = {k.lower(): k for k in (rdr.fieldnames or [])}
        need = ("locus", "sequence", "allele_code")
        if any(k not in cols for k in need):
            raise ValueError(f"Catalog CSV missing required columns: {need}")
        for row in rdr:
            loc = row[cols["locus"]]
            seq = row[cols["sequence"]]
            code = int(row[cols["allele_code"]])
            catalog[loc][seq] = code

    # ---- phased definitions from catalog (purely catalog-derived; no new alleles) ----
    phased_maps, def_rows = build_phased_definitions_from_catalog(catalog)
    phased_def_path = Path(f"{out_prefix}_phased_def.tsv")
    with phased_def_path.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus", "var_index", "type", "ref_pos", "ref_base", "alts"])
        for r in def_rows:
            w.writerow([r["locus"], r["var_index"], r["type"], r["ref_pos"], r["ref_base"], r["alts"]])

    phased_rows: list[dict] = []  # per-sample phased hap strings

    def _fmt_phased(h1: str, h2: str) -> str:
        """Always emit two hap fields separated by '|'.

        Missing/undefined (including loci with 0 SNPs) is represented as '.|.'.
        """
        h1 = (h1 or "").strip()
        h2 = (h2 or "").strip()
        if not h1 or not h2:
            return ".|."
        return f"{h1}|{h2}"

    # ---- iterate resolved FASTQs ----
    resolved_paths = sorted(p for p in Path(resolved_dir).glob("resolved_*.fastq.gz") if not str(p).endswith(".part"))

    per_locus_sample_gt: dict[str, dict[str, str]] = defaultdict(dict)
    a2_rows: list[dict] = []

    # for plots
    locus_points_a1a2 = defaultdict(list)        # (a1_ct, a2_ct, call, code_for_hom_or_None)
    locus_points_depth_pct = defaultdict(list)   # (depth, a2_pct, call, code_for_hom_or_None)
    locus_called_samples = defaultdict(int)
    locus_total_samples  = defaultdict(int)
    locus_background_sum = defaultdict(float)    # sum over samples of background fraction (later / n)

    loci_iter_master = list(locus_universe) if locus_universe else sorted(catalog.keys())

    for fpath in PROG(resolved_paths, desc="Second pass calling", unit="sample", total=len(resolved_paths)):
        sample = fpath.name.replace("resolved_", "").replace(".fastq.gz", "")

        # count catalog vs off-catalog for this sample
        sample_locus_cnts: dict[str, Counter] = defaultdict(Counter)
        sample_locus_offcat: dict[str, int] = defaultdict(int)

        with gzip.open(fpath, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip()
                _ = fh.readline(); _ = fh.readline()

                if "|locus=" not in h:
                    continue
                loc = h.strip().split("|locus=")[-1]

                if loc in catalog and seq in catalog[loc]:
                    sample_locus_cnts[loc][seq] += 1
                else:
                    # still within a locus bin, but not present in catalog for that locus
                    sample_locus_offcat[loc] += 1

        # background fraction per locus for this sample
        per_sample_background_frac: dict[str, float] = {}

        for loc in loci_iter_master:
            cnts = sample_locus_cnts.get(loc, Counter())
            total_cat = sum(cnts.values())
            off_cat = sample_locus_offcat.get(loc, 0)

            if (total_cat + off_cat) > 0:
                locus_total_samples[loc] += 1

            # no catalog reads -> no call
            if total_cat == 0:
                per_locus_sample_gt[loc][sample] = "000000"
                a2_rows.append({
                    "sample": sample, "locus": loc, "depth": 0,
                    "A1_count": 0, "A2_count": 0,
                    "A1_code": 0, "A2_code": 0,
                    "A2_pct": 0.0, "call": "LOW"
                })
                locus_points_a1a2[loc].append((0, 0, "LOW", None))
                locus_points_depth_pct[loc].append((0, 0.0, "LOW", None))

                if (total_cat + off_cat) > 0:
                    per_sample_background_frac[loc] = 1.0 if off_cat > 0 else 0.0
                continue

            top = cnts.most_common(2)
            a1_seq, a1_ct = top[0]
            a1_code = catalog[loc][a1_seq]

            if len(top) > 1:
                a2_seq, a2_ct = top[1]
                a2_code = catalog[loc][a2_seq]
            else:
                a2_ct, a2_code = 0, None

            denom = a1_ct + a2_ct
            frac_a2 = (a2_ct / denom) if denom > 0 else 0.0
            a2_pct = frac_a2 * 100.0

            if denom < min_depth:
                call = "LOW"
                gt = "000000"
            else:
                if a2_code is not None and frac_a2 >= a2_hi_frac:
                    call = "HET"
                    gt = f"{min(a1_code, a2_code)}{max(a1_code, a2_code)}"
                elif frac_a2 <= a2_lo_frac:
                    call = "HOM"
                    gt = f"{a1_code}{a1_code}"
                else:
                    call = "NC"
                    gt = "000000"

            if call in ("HOM", "HET"):
                locus_called_samples[loc] += 1

            per_locus_sample_gt[loc][sample] = gt

            # phased hap string outputs (derived from catalog allele sequences)
            pm = phased_maps.get(loc)
            if call in ("HOM", "HET") and pm is not None:
                if call == "HOM":
                    c1 = int(a1_code); c2 = int(a1_code)
                else:
                    c1 = int(min(a1_code, a2_code)); c2 = int(max(a1_code, a2_code))
                h1_snp = pm["code_to_hap_snp"].get(c1, pm["code_to_hap_snp"].get(str(c1), ""))
                h2_snp = pm["code_to_hap_snp"].get(c2, pm["code_to_hap_snp"].get(str(c2), ""))
                h1_all = pm["code_to_hap_all"].get(c1, pm["code_to_hap_all"].get(str(c1), ""))
                h2_all = pm["code_to_hap_all"].get(c2, pm["code_to_hap_all"].get(str(c2), ""))
                phased_rows.append({
                    "sample": sample, "locus": loc, "gt": gt,
                    "hap_snp": _fmt_phased(h1_snp, h2_snp),
                    "hap_all": _fmt_phased(h1_all, h2_all),
                    "n_snp": len(pm["var_cols_snp"]),
                    "n_var": len(pm["var_cols_all"]),
                })
            else:
                phased_rows.append({
                    "sample": sample, "locus": loc, "gt": gt,
                    "hap_snp": ".|.", "hap_all": ".|.",
                    "n_snp": 0, "n_var": 0,
                })

            a2_rows.append({
                "sample": sample, "locus": loc, "depth": denom,
                "A1_count": a1_ct, "A2_count": a2_ct,
                "A1_code": int(a1_code) if a1_code is not None else 0,
                "A2_code": int(a2_code) if a2_code is not None else 0,
                "A2_pct": round(a2_pct, 2), "call": call
            })

            code_for_hom = int(a1_code) if call == "HOM" else None
            locus_points_a1a2[loc].append((a1_ct, a2_ct, call, code_for_hom))
            locus_points_depth_pct[loc].append((denom, a2_pct, call, code_for_hom))

            # background fraction: (off-catalog + other-catalog) / total_reads_for_loc
            other_catalog = total_cat - denom
            total_reads_for_loc = total_cat + off_cat
            if total_reads_for_loc > 0:
                per_sample_background_frac[loc] = (off_cat + other_catalog) / total_reads_for_loc

        # accumulate background
        for loc, frac in per_sample_background_frac.items():
            locus_background_sum[loc] += float(frac)

    # ---- write outputs ----
    with Path(f"{out_prefix}_a2_metrics.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "locus", "depth", "A1_count", "A2_count", "A1_code", "A2_code", "A2_pct", "call"])
        for r in a2_rows:
            w.writerow([r["sample"], r["locus"], r["depth"], r["A1_count"], r["A2_count"], r["A1_code"], r["A2_code"], r["A2_pct"], r["call"]])

    all_samples = sorted({s for d in per_locus_sample_gt.values() for s in d})
    loci_to_write = sorted(loci_iter_master)

    with Path(f"{out_prefix}_genotypes.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus"] + all_samples)
        for loc in loci_to_write:
            d = per_locus_sample_gt.get(loc, {})
            w.writerow([loc] + [d.get(s, "000000") for s in all_samples])

    with Path(f"{out_prefix}_phased.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "locus", "gt", "hap_snp", "hap_all", "n_snp", "n_var"])
        for r in phased_rows:
            w.writerow([r["sample"], r["locus"], r["gt"], r["hap_snp"], r["hap_all"], r["n_snp"], r["n_var"]])

    # phased VAR hap matrix (SNP + indel sites; locus x sample)
    # Uses hap_all, but removes comma delimiters so each haplotype is compact.
    phased_hap_var: dict[str, dict[str, str]] = defaultdict(dict)

    def compact_hap_all(v: str) -> str:
        """Convert hap_all format from comma-delimited to compact hap strings."""
        if not v or v in {".", "|", ".|", "|."}:
            return ".|."
        return v.replace(",", "")

    for r in phased_rows:
        phased_hap_var[r["locus"]][r["sample"]] = compact_hap_all(r.get("hap_all", ".|."))

    with Path(f"{out_prefix}_phasedVARs.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus"] + all_samples)

        for loc in loci_to_write:
            d = phased_hap_var.get(loc, {})
            row = []
            for s in all_samples:
                v = compact_hap_all(d.get(s, ".|."))
                row.append(v)
            w.writerow([loc] + row)

    # ---- plots (optional) ----
    if plots_dir is not None and HAVE_MPL:
        plots_dir = Path(plots_dir)
        plots_dir.mkdir(parents=True, exist_ok=True)

        loci_sorted = locus_stats.get("sorted_loci_by_mean", [])
        locus_pct_mean = locus_stats.get("mean_pct", {})
        locus_pct_sd = locus_stats.get("sd_pct", {})

        assay_exclude = set(assay_exclude_loci or [])
        for loc in assay_exclude:
            stale = Path(plots_dir) / f"{loc}_dashboard.png"
            if stale.exists():
                try:
                    stale.unlink()
                except Exception:
                    pass
        loci_to_plot = [loc for loc in loci_to_write if loc not in assay_exclude]

        for loc in PROG(loci_to_plot, desc="Plotting per-locus dashboards", unit="locus", total=len(loci_to_plot)):
            allele_freq_codes = _allele_freq_from_genos(per_locus_sample_gt.get(loc, {}))
            tot = locus_total_samples.get(loc, 0)
            called = locus_called_samples.get(loc, 0)
            call_rate = (called / tot * 100.0) if tot else 0.0
            bg_pct = (locus_background_sum.get(loc, 0.0) / tot * 100.0) if tot else 0.0

            plot_locus_dashboard(
                locus=loc,
                points_a1a2=locus_points_a1a2.get(loc, []),
                depths_a2pct=locus_points_depth_pct.get(loc, []),
                loci_sorted_by_share=loci_sorted,
                locus_pct_mean=locus_pct_mean,
                locus_pct_sd=locus_pct_sd,
                allele_freq_codes=allele_freq_codes,
                outdir=plots_dir,
                extra_info={"call_rate": call_rate, "background_pct": bg_pct},
                a2_lo=a2_lo,
                a2_hi=a2_hi,
            )

    return {"per_locus_sample_gt": per_locus_sample_gt, "all_samples": all_samples}

def plot_locus_dashboard(
    locus,
    points_a1a2,
    depths_a2pct,
    loci_sorted_by_share,
    locus_pct_mean,
    locus_pct_sd,
    allele_freq_codes,
    outdir,
    extra_info=None,
    a2_lo=10.0,
    a2_hi=25.0,
):
    if not HAVE_MPL:
        return

    outdir = Path(outdir) if outdir else Path("plots")
    outdir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12, 8.6))
    gs = GridSpec(2, 2, height_ratios=[1.0, 1.4], width_ratios=[1.0, 1.0],
                  hspace=0.35, wspace=0.28)
    fig.subplots_adjust(top=0.88)

    het_style = dict(marker="o", linestyle="none", alpha=0.9, color="C1", rasterized=True)
    nc_style  = dict(marker="x", linestyle="none", alpha=0.9, color="k",  rasterized=True)
    low_style = dict(marker="^", linestyle="none", alpha=0.9, color="red", rasterized=True)

    markers_cycle = ['o','s','D','^','v','P','X','*','<','>','h','H']
    colors_cycle  = ["#2ca02c","#d62728","#9467bd","#1f77b4","#ff7f0e","#17becf",
                     "#8c564b","#e377c2","#7f7f7f","#bcbd22","#aec7e8","#ffbb78"]
    explicit = {101: ('o', "#2ca02c"), 102: ('s', "#d62728"), 103: ('D', "#9467bd"),
                201: ('o', "#1f77b4"), 202: ('s', "#ff7f0e"), 203: ('D', "#17becf")}

    present_codes = sorted(allele_freq_codes.keys()) if allele_freq_codes else []
    palette = [(m, c) for m in markers_cycle for c in colors_cycle]
    code_style = {}
    for i, cd in enumerate(present_codes):
        code_style[cd] = palette[i % len(palette)]
    for cd, sty in explicit.items():
        if cd in code_style:
            code_style[cd] = sty

    def plot_hom(ax, xvals, yvals, codes):
        groups = collections.defaultdict(lambda: ([], []))
        for x, y, cd in zip(xvals, yvals, codes):
            if cd is None:
                continue
            groups[int(cd)][0].append(x)
            groups[int(cd)][1].append(y)
        for cd in sorted(groups):
            m, c = code_style.get(cd, ('o', "#2ca02c"))
            xs, ys = groups[cd]
            ax.scatter(xs, ys, marker=m, c=c, edgecolors='none', linewidths=0,
                       s=36, alpha=0.95, rasterized=True)

    def _iter_points4(points):
        # Accept (x,y,call) or (x,y,call,code)
        for t in (points or []):
            if len(t) == 4:
                yield t
            elif len(t) == 3:
                x, y, call = t
                yield (x, y, call, None)

    # ---------------------------
    # Top-left: A1 vs A2 counts
    # ---------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    pts = list(_iter_points4(points_a1a2))
    xs = [x for (x, y, call, code) in pts if call == "HOM"]
    ys = [y for (x, y, call, code) in pts if call == "HOM"]
    cs = [code for (x, y, call, code) in pts if call == "HOM"]
    if xs:
        plot_hom(ax1, xs, ys, cs)

    xs = [x for (x, y, call, _) in pts if call == "HET"]
    ys = [y for (x, y, call, _) in pts if call == "HET"]
    if xs:
        ax1.plot(xs, ys, **het_style)

    xs = [x for (x, y, call, _) in pts if call == "NC"]
    ys = [y for (x, y, call, _) in pts if call == "NC"]
    if xs:
        ax1.plot(xs, ys, **nc_style)

    xs = [x for (x, y, call, _) in pts if call == "LOW"]
    ys = [y for (x, y, call, _) in pts if call == "LOW"]
    if xs:
        ax1.plot(xs, ys, **low_style)

    ax1.set_xlabel("A1 count")
    ax1.set_ylabel("A2 count")
    ax1.set_title(f"{locus}: A1 vs A2")

    # Use identical x/y limits so A2 noise is not visually exaggerated
    # when loci are mostly or entirely homozygous. The gray x=y line
    # provides a visual reference for balanced heterozygotes.
    max_val = 0
    for x, y, _call, _code in pts:
        if x > max_val:
            max_val = x
        if y > max_val:
            max_val = y
    max_val = int(max_val * 1.05) + 1

    ax1.set_xlim(0, max_val)
    ax1.set_ylim(0, max_val)
    ax1.plot([0, max_val], [0, max_val], color="gray", linestyle="--", linewidth=1, alpha=0.7, zorder=0)

    # ---------------------------
    # Top-right: depth vs %A2
    # ---------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    pts2 = list(_iter_points4(depths_a2pct))
    xs = [d for (d, p, call, code) in pts2 if call == "HOM"]
    ys = [p for (d, p, call, code) in pts2 if call == "HOM"]
    cs = [code for (d, p, call, code) in pts2 if call == "HOM"]
    if xs:
        plot_hom(ax2, xs, ys, cs)

    xs = [d for (d, p, call, _) in pts2 if call == "HET"]
    ys = [p for (d, p, call, _) in pts2 if call == "HET"]
    if xs:
        ax2.plot(xs, ys, **het_style)

    xs = [d for (d, p, call, _) in pts2 if call == "NC"]
    ys = [p for (d, p, call, _) in pts2 if call == "NC"]
    if xs:
        ax2.plot(xs, ys, **nc_style)

    xs = [d for (d, p, call, _) in pts2 if call == "LOW"]
    ys = [p for (d, p, call, _) in pts2 if call == "LOW"]
    if xs:
        ax2.plot(xs, ys, **low_style)

    ax2.axhline(a2_lo, ls="--", lw=1, color="gray")
    ax2.axhline(a2_hi, ls="--", lw=1, color="gray")
    ax2.set_ylim(0, 50)
    ax2.set_xlabel("Total locus depth")
    ax2.set_ylabel("% A2")
    ax2.set_title(f"{locus}: depth vs %A2")
    ax2.set_xlim(left=0)

    # ---------------------------
    # Center annotation
    # ---------------------------
    if extra_info is None:
        extra_info = {}
    cr = extra_info.get("call_rate", None)
    bg = extra_info.get("background_pct", None)
    label = " | ".join(
        [f"Call rate: {cr:0.1f}%" if cr is not None else "",
         f"Background: {bg:0.1f}%" if bg is not None else ""]
    ).strip(" |")
    if label:
        fig.text(
            0.50, 0.545, label,
            ha="center", va="center",
            fontsize=12, color="dimgray", fontweight="semibold",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="round,pad=0.25"),
            transform=fig.transFigure
        )

    # ---------------------------
    # Bottom-left: read distribution among loci
    # ---------------------------
    ax3 = fig.add_subplot(gs[1, 0])

    means = [locus_pct_mean[L] for L in loci_sorted_by_share]
    sds   = [locus_pct_sd[L]   for L in loci_sorted_by_share]
    xs3 = list(range(len(loci_sorted_by_share)))

    n_loci = max(1, len(loci_sorted_by_share))
    total = sum(means) if means else 0.0
    if total <= 5.0:
        means = [m * 100.0 for m in means]
        sds   = [sd * 100.0 for sd in sds]

    avg_share = 100.0 / n_loci

    highlight_idx = None
    if locus in loci_sorted_by_share:
        highlight_idx = loci_sorted_by_share.index(locus)

    # Background bars
    bars = ax3.bar(
        xs3,
        means,
        color="#b3b3b3",
        edgecolor="#b3b3b3",
        linewidth=0.0,
        zorder=1
    )

    # Background error bars (gray)
    if sds:
        ax3.errorbar(
            xs3,
            means,
            yerr=sds,
            fmt="none",
            ecolor="#777777",
            elinewidth=0.8,
            capthick=0.8,
            capsize=1,
            zorder=2
        )

    # Highlighted locus bar (solid red) + black error bar
    if highlight_idx is not None:
        bars[highlight_idx].set_facecolor("red")
        bars[highlight_idx].set_edgecolor("red")
        bars[highlight_idx].set_linewidth(1.0)
        bars[highlight_idx].set_zorder(4)

        if sds:
            ax3.errorbar(
                highlight_idx,
                means[highlight_idx],
                yerr=sds[highlight_idx],
                fmt="none",
                ecolor="black",
                elinewidth=1.2,
                capthick=1.2,
                capsize=2,
                zorder=5
            )

    ax3.axhline(avg_share, color="tab:red", lw=1, zorder=6)

    ymax = max((m + sd) for m, sd in zip(means or [0.0], sds or [0.0]))
    ax3.set_ylim(0, max(ymax * 1.15, avg_share * 1.8, 2.0))
    ax3.set_xlim(-0.5, len(xs3) - 0.5 + 1.0)
    ax3.set_ylabel("Avg % reads per locus (±1 SD)")
    ax3.set_xlabel("Loci (sorted by mean %, least → most)")
    ax3.set_title("Read distribution among loci")
    if len(xs3) > 20:
        ax3.set_xticks([])

    # ---------------------------
    # Bottom-right: allele frequency
    # ---------------------------
    ax4 = fig.add_subplot(gs[1, 1])
    freqs_dict = allele_freq_codes or {}
    codes = sorted(freqs_dict.keys())
    freqs = [freqs_dict[c] for c in codes]
    if codes:
        s = sum(freqs)
        if s and abs(s - 1.0) > 1e-12:
            kmax = max(range(len(freqs)), key=lambda i: freqs[i])
            freqs[kmax] += (1.0 - s)
        ax4.bar(range(len(codes)), freqs, width=0.8)
        ax4.set_xticks(range(len(codes)))
        ax4.set_xticklabels([str(c) for c in codes])

    ax4.set_ylim(0, 1.0)
    ax4.set_ylabel("allele frequency")
    ax4.set_xlabel("Allele code")
    ax4.set_title(f"{locus}: allele frequency")

    outpath = outdir / f"{locus}_dashboard.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


# ----------------------------------------------------------------------
# Phased microhap string definitions (from catalog alleles)
# ----------------------------------------------------------------------










def _nw_global_align(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    """
    Simple Needleman-Wunsch global alignment (no affine gaps).
    Returns (a_aln, b_aln) with '-' gaps. Works well for short amplicons (~50-400bp).
    """
    a = str(a)
    b = str(b)
    n = len(a)
    m = len(b)
    # DP score + traceback
    # tb: 0=diag, 1=up (gap in b), 2=left (gap in a)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        tb[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        tb[0][j] = 2

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            s_diag = score[i - 1][j - 1] + (match if ai == bj else mismatch)
            s_up = score[i - 1][j] + gap
            s_left = score[i][j - 1] + gap
            best = s_diag
            t = 0
            if s_up > best:
                best = s_up
                t = 1
            if s_left > best:
                best = s_left
                t = 2
            score[i][j] = best
            tb[i][j] = t

    # traceback
    i, j = n, m
    a_aln = []
    b_aln = []
    while i > 0 or j > 0:
        t = tb[i][j] if (i > 0 and j > 0) else (1 if i > 0 else 2)
        if i > 0 and j > 0 and t == 0:
            a_aln.append(a[i - 1])
            b_aln.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or t == 1):
            a_aln.append(a[i - 1])
            b_aln.append('-')
            i -= 1
        else:
            a_aln.append('-')
            b_aln.append(b[j - 1])
            j -= 1

    return ''.join(reversed(a_aln)), ''.join(reversed(b_aln))


def _gaps_per_boundary(gapped_ref: str, ref_len: int):
    """
    Given an alignment string for the reference with '-' gaps, compute how many
    gaps occur at each boundary between reference bases.
    Returns list length ref_len+1: gaps before base0, between bases, after last base.
    """
    gaps = [0] * (ref_len + 1)
    ref_i = 0
    # leading gaps are boundary 0
    k = 0
    while k < len(gapped_ref) and gapped_ref[k] == '-':
        gaps[0] += 1
        k += 1
    # now parse columns
    while k < len(gapped_ref):
        c = gapped_ref[k]
        if c == '-':
            # gaps before next ref base (or trailing)
            if ref_i <= ref_len:
                gaps[ref_i] += 1
        else:
            ref_i += 1
        k += 1
    # trailing gaps are already counted into gaps[ref_len] by the loop above
    return gaps


def _project_alt_to_master(gapped_ref: str, gapped_alt: str, ref: str, master_gaps):
    """
    Project a single allele alignment onto the master gap scheme.
    Returns master_alt (same length as master_ref).
    """
    ref_len = len(ref)
    # Build per-boundary insertion strings and per-base aligned char
    insertions = [''] * (ref_len + 1)
    per_base = [''] * ref_len

    ref_i = 0
    for rch, ach in zip(gapped_ref, gapped_alt):
        if rch == '-':
            # insertion relative to reference at current boundary (before ref base ref_i)
            insertions[ref_i] += ach
        else:
            # aligned to a reference base
            if ref_i < ref_len:
                per_base[ref_i] = ach
            ref_i += 1

    # Fill any unset per_base with '-' (shouldn't happen, but safe)
    for i in range(ref_len):
        if per_base[i] == '':
            per_base[i] = '-'

    # Now build master-alt with padding insertions to master_gaps
    out = []
    # boundary 0 insertions
    ins0 = insertions[0]
    out.append(ins0)
    if len(ins0) < master_gaps[0]:
        out.append('-' * (master_gaps[0] - len(ins0)))
    # for each base: base col then boundary after it
    for i in range(ref_len):
        out.append(per_base[i])
        ins = insertions[i + 1]
        out.append(ins)
        if len(ins) < master_gaps[i + 1]:
            out.append('-' * (master_gaps[i + 1] - len(ins)))
    return ''.join(out)

def _right_normalize_pairwise_alignment(gref: str, galt: str) -> tuple[str, str]:
    """
    Right-normalize equivalent indel placements in a pairwise alignment.

    This helps prevent ambiguous left/right placement of indels in repetitive
    sequence from changing phased SNP strings. It shifts gap runs to the right
    whenever doing so preserves the alignment equivalently.
    """
    r = list(gref)
    a = list(galt)
    n = len(r)

    changed = True
    while changed:
        changed = False
        i = 0
        while i < n - 1:

            # Case 1: deletion in alt (gap run in alt, bases in ref)
            if a[i] == "-" and r[i] != "-":
                j = i
                while j < n and a[j] == "-" and r[j] != "-":
                    j += 1

                # shift right while the next alt base matches the first deleted ref base
                while j < n and r[j] != "-" and a[j] == r[i]:
                    block = a[i:j+1]
                    a[i:j+1] = [a[j]] + block[:-1]
                    i += 1
                    j += 1
                    changed = True

                i = j
                continue

            # Case 2: insertion in alt (gap run in ref, bases in alt)
            if r[i] == "-" and a[i] != "-":
                j = i
                while j < n and r[j] == "-" and a[j] != "-":
                    j += 1

                # shift right while the next ref base matches the first inserted alt base
                while j < n and a[j] != "-" and r[j] == a[i]:
                    block = r[i:j+1]
                    r[i:j+1] = [r[j]] + block[:-1]
                    i += 1
                    j += 1
                    changed = True

                i = j
                continue

            i += 1

    return "".join(r), "".join(a)

def build_phased_definitions_from_catalog(catalog):
    """
    Build per-locus phased variant definitions (SNPs + small INDELs) by aligning
    catalog allele sequences for each locus. This is designed to avoid the
    "indel causes shift -> everything after looks like SNP" failure mode.

    Works with two catalog shapes:
      1) Dict-locus -> Dict-seq -> allele_code (what call_from_resolved_fastqs builds)
      2) Iterable of dict/tuple rows with keys: locus, allele_code, sequence
    """
    by_locus = {}

    # Normalize catalog into {locus: [(code, seq), ...]}
    if isinstance(catalog, dict):
        # common in this script: catalog[locus] = {sequence: code, ...}
        for locus, seq_to_code in catalog.items():
            try:
                items = [(int(code), str(seq)) for (seq, code) in seq_to_code.items()]
            except Exception:
                # maybe it's already code->seq
                items = []
                for k, v in seq_to_code.items():
                    try:
                        items.append((int(k), str(v)))
                    except Exception:
                        continue
            if items:
                by_locus[str(locus)] = items
    else:
        for r in catalog:
            if isinstance(r, dict):
                locus = str(r.get("locus", ""))
                seq = str(r.get("sequence", ""))
                try:
                    code = int(r.get("allele_code", 0))
                except Exception:
                    continue
            else:
                # tuple/list
                try:
                    locus = str(r[0]); code = int(r[1]); seq = str(r[2])
                except Exception:
                    continue
            if locus and seq and code:
                by_locus.setdefault(locus, []).append((code, seq))

    phased_maps = {}
    def_rows = []

    for locus, alleles in by_locus.items():
        # Need at least 2 unique sequences to define variants
        uniq = {}
        for code, seq in alleles:
            uniq.setdefault(seq, code)  # keep a code for each seq
        if len(uniq) < 2:
            # still create an empty phased_map so downstream has keys
            phased_maps[locus] = {
                "code_to_hap_snp": {},
                "code_to_hap_all": {},
                "var_cols_snp": [],
                "var_cols_all": [],
            }
            continue

        # Choose reference allele: prefer code 101 if present, else smallest code
        alleles_sorted = sorted(alleles, key=lambda x: x[0])
        ref_seq = None
        for code, seq in alleles_sorted:
            if code == 101:
                ref_seq = seq
                break
        if ref_seq is None:
            ref_seq = alleles_sorted[0][1]

        ref = ref_seq
        ref_len = len(ref)

        # Align each allele to ref and collect gap schemes
        per_allele_alignment = {}  # code -> (gapped_ref, gapped_alt)
        gap_schemes = []  # list of gaps_per_boundary
        for code, seq in alleles_sorted:
            gref, galt = _nw_global_align(ref, seq)
            gref, galt = _right_normalize_pairwise_alignment(gref, galt)
            per_allele_alignment[code] = (gref, galt)
            gap_schemes.append(_gaps_per_boundary(gref, ref_len))

        # Master gaps at each boundary = max across alleles
        master_gaps = [0] * (ref_len + 1)
        for gs in gap_schemes:
            for i, v in enumerate(gs):
                if v > master_gaps[i]:
                    master_gaps[i] = v

        # Build master reference with union of gaps
        master_ref_parts = []
        master_ref_parts.append('-' * master_gaps[0])
        for i in range(ref_len):
            master_ref_parts.append(ref[i])
            master_ref_parts.append('-' * master_gaps[i + 1])
        master_ref = ''.join(master_ref_parts)
        aln_len = len(master_ref)

        # Project all alleles into master alignment columns
        code_to_master_alt = {}
        for code, (gref, galt) in per_allele_alignment.items():
            code_to_master_alt[code] = _project_alt_to_master(gref, galt, ref, master_gaps)

        # Identify variant columns (SNP and INDEL) and build definitions
        var_cols_all = []
        var_cols_snp = []
        # track reference coordinate (1-based) as we walk master columns
        ref_pos = 0
        var_index = 0

        # First determine which columns are variable
        for col in range(aln_len):
            rch = master_ref[col]
            # update ref coordinate
            if rch != '-':
                ref_pos += 1

            # gather alleles at this column
            col_chars = set()
            for alt in code_to_master_alt.values():
                col_chars.add(alt[col])

            # exclude gaps in variability test carefully:
            # - if rch is a base: SNP if non-gap bases differ; INDEL if any gap present and any base present
            # - if rch is '-': insertion column; variant if any non-gap base exists
            is_var = False
            vtype = None
            if rch == '-':
                # insertion relative to ref
                if any(c != '-' for c in col_chars):
                    is_var = True
                    vtype = "INDEL"
            else:
                non_gap = {c for c in col_chars if c != '-'}
                if len(non_gap) > 1:
                    is_var = True
                    vtype = "SNP"
                elif ('-' in col_chars) and len(non_gap) >= 1:
                    is_var = True
                    vtype = "INDEL"

            if not is_var:
                continue

            var_cols_all.append(col)
            if vtype == "SNP":
                var_cols_snp.append(col)


        # Collapse multi-base INDEL runs down to a single "event" column.
        # For an insertion relative to the chosen reference (master_ref has '-'),
        # we keep only the FIRST inserted base as the allele symbol for that event.
        # For a deletion (some alleles have '-' while master_ref has bases),
        # we keep only the FIRST deleted base position as the event anchor.
        #
        # This prevents a small INDEL from "shifting" the rest of the alignment into a
        # cascade of false SNPs, and matches the desired representation where e.g.
        # ATC/- becomes A/- in the phased definitions / outputs.
        def _collapse_variant_columns(master_ref, code_to_master_alt, var_cols_all):
            cols = sorted(var_cols_all)
            if not cols:
                return [], {}

            def gap_code_set(col):
                return {code for code, alt in code_to_master_alt.items() if alt[col] == "-"}

            event_starts = []
            event_runs = {}  # start_col -> list of cols in the run
            i = 0
            while i < len(cols):
                c0 = cols[i]
                ref0 = master_ref[c0]
                gaps0 = gap_code_set(c0)
                is_indel0 = (ref0 == "-") or bool(gaps0)

                if not is_indel0:
                    event_starts.append(c0)
                    event_runs[c0] = [c0]
                    i += 1
                    continue

                # Insertion run: master_ref is '-' across consecutive columns
                if ref0 == "-":
                    j = i + 1
                    while j < len(cols):
                        c = cols[j]
                        if c != cols[j - 1] + 1:
                            break
                        if master_ref[c] != "-":
                            break
                        j += 1
                    run = cols[i:j]
                    event_starts.append(run[0])
                    event_runs[run[0]] = run
                    i = j
                    continue

                # Deletion run: master_ref has bases, but the same set of alleles are gapped
                # across consecutive columns. Collapse to the first deleted base.
                j = i + 1
                while j < len(cols):
                    c = cols[j]
                    if c != cols[j - 1] + 1:
                        break
                    if master_ref[c] == "-":
                        break
                    gaps = gap_code_set(c)
                    if gaps != gaps0:
                        break
                    j += 1
                run = cols[i:j]
                event_starts.append(run[0])
                event_runs[run[0]] = run
                i = j

            return event_starts, event_runs

        # Collapse any INDEL runs in the "all-variants" column set
        event_cols_all, event_runs = _collapse_variant_columns(master_ref, code_to_master_alt, var_cols_all)

        # Build ref-position mapping for the chosen reference alignment
        master_ref_to_refpos = {}
        ref_i = 0
        for col in range(aln_len):
            if master_ref[col] != "-":
                ref_i += 1
                master_ref_to_refpos[col] = ref_i
            else:
                master_ref_to_refpos[col] = None

        def _prev_refpos(col):
            # for insertions (ref '-') anchor to previous real ref base (or 0 if none)
            c = col - 1
            while c >= 0:
                rp = master_ref_to_refpos.get(c)
                if rp is not None:
                    return rp
                c -= 1
            return 0

        # Build def_rows for event columns (SNPs and collapsed INDEL events)
        def_rows_locus = []
        var_index = 0
        for col in event_cols_all:
            ref_base = master_ref[col]
            alts = set()
            any_gap = False
            for code, alt in code_to_master_alt.items():
                b = alt[col]
                if b == "-":
                    any_gap = True
                if b != ref_base:
                    alts.add(b)

            if not alts:
                continue

            var_index += 1
            if ref_base == "-" or any_gap:
                vtype = "INDEL"
                ref_pos = _prev_refpos(col) if ref_base == "-" else (master_ref_to_refpos[col] or 0)
            else:
                vtype = "SNP"
                ref_pos = master_ref_to_refpos[col] or 0

            def_rows.append({
                "locus": locus,
                "var_index": var_index,
                "type": vtype,
                "ref_pos": ref_pos,
                "ref_base": ref_base,
                "alts": ",".join(sorted(alts)),
            })
            # after finishing the locus
            def_rows.extend(def_rows_locus)

        # Build per-allele hap strings for SNP-only columns and ALL-variant event columns
        def _hap_for_cols(code: int, cols, delim: str):
            alt = code_to_master_alt[code]
            toks = []
            for col in cols:
                ach = alt[col]
                if ach == '-':
                    toks.append("-")  # deletion / no insertion
                else:
                    toks.append(ach)
            if not toks:
                return "."
            return delim.join(toks)

        # decide delimiter: if SNP-only and all SNP tokens are single-base, use "" (historic),
        # else use "," to avoid ambiguous variable-length strings.
        snp_delim = ""  # classic
        all_delim = ""
        # if any INDEL present, force comma for ALL and also for SNPs if desired? keep SNP classic.
        if len(var_cols_all) != len(var_cols_snp):
            all_delim = ","
        else:
            all_delim = ""  # SNP-only locus

        code_to_hap_snp = {}
        code_to_hap_all = {}
        for code, _seq in alleles_sorted:
            if code not in code_to_master_alt:
                continue
            hs = _hap_for_cols(code, var_cols_snp, snp_delim) if var_cols_snp else "."
            ha = _hap_for_cols(code, event_cols_all, all_delim) if event_cols_all else "."
            code_to_hap_snp[str(code)] = hs
            if str(code).isdigit():
                code_to_hap_snp[int(code)] = hs
            code_to_hap_all[str(code)] = ha
            if str(code).isdigit():
                code_to_hap_all[int(code)] = ha

        phased_maps[locus] = {
            "code_to_hap_snp": code_to_hap_snp,
            "code_to_hap_all": code_to_hap_all,
            "var_cols_snp": var_cols_snp,
            "var_cols_all": event_cols_all,
            "var_cols_all_raw": var_cols_all,
            "indel_event_runs": event_runs,
        }

    return phased_maps, def_rows



def write_placeholder_png(path: Path, title: str, message: str):
    """Write a simple placeholder PNG with an error message."""
    if not HAVE_MPL:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8, 4.5))
    ax = fig.add_subplot(111)
    ax.axis('off')
    ax.text(0.5, 0.65, title, ha='center', va='center', fontsize=14, fontweight='bold')
    ax.text(0.5, 0.40, message, ha='center', va='center', fontsize=10, wrap=True)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)

# ---------------------------------------------------------------------
# library summary plots (unchanged)
# ---------------------------------------------------------------------
def plot_library_summary(
    sample_stats_tsv,
    a2_metrics_tsv,
    primers_csv,
    out_prefix="microhap",
):
    import matplotlib.pyplot as plt
    from collections import defaultdict

    sample_stats_tsv = Path(sample_stats_tsv)
    a2_metrics_tsv = Path(a2_metrics_tsv)
    primers_csv = Path(primers_csv)

    def _norm_sample(s: str) -> str:
        s = s.strip()
        s = re.sub(r"\.fastq(\.gz)?$", "", s, flags=re.IGNORECASE)
        return s

    with primers_csv.open() as f:
        num_loci = sum(1 for _ in f) - 1

    raw_reads: dict[str, int] = {}
    pb_reads: dict[str, int] = {}
    with sample_stats_tsv.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            s = _norm_sample(row["sample"])
            raw_reads[s] = int(row["raw_reads"])
            pb_reads[s] = int(row["primer_bounded_reads"])

    called = defaultdict(int)
    with a2_metrics_tsv.open() as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r, None)
        if header is None:
            return
        sample_idx = 0
        call_idx = 8
        for row in r:
            if not row:
                continue
            s = _norm_sample(row[sample_idx])
            call = row[call_idx]
            if call in ("HET", "HOM"):
                called[s] += 1

    all_samples = sorted(set(raw_reads.keys()) | set(pb_reads.keys()) | set(called.keys()))
    if not all_samples or num_loci <= 0:
        return

    x_raw = []
    x_pb = []
    y_gt = []
    ratio_pb_raw = []
    raw_list = []
    pb_list = []

    for s in all_samples:
        rr = raw_reads.get(s, 0)
        pr = pb_reads.get(s, 0)
        n_called = called.get(s, 0)
        gt_pct = 100.0 * n_called / num_loci if num_loci > 0 else 0.0

        x_raw.append(rr / 1000.0)
        x_pb.append(pr / 1000.0)
        y_gt.append(gt_pct)
        raw_list.append(rr / 1000.0)
        pb_list.append(pr / 1000.0)
        ratio_pb_raw.append(pr / rr if rr > 0 else 0.0)

    fig1, ax1 = plt.subplots(figsize=(6, 4))
    ax1.hist(ratio_pb_raw, bins=30)
    ax1.set_xlabel("Primer-bounded / raw reads")
    ax1.set_ylabel("Number of samples")
    ax1.set_title("On-target fraction per sample")
    fig1.tight_layout()
    fig1.savefig(f"{out_prefix}_on_target_fraction.png", dpi=300)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.hist(raw_list, bins=30)
    ax2.set_xlabel("Raw reads (K)")
    ax2.set_ylabel("Number of samples")
    ax2.set_title("Raw reads per sample")
    fig2.tight_layout()
    fig2.savefig(f"{out_prefix}_raw_reads_hist.png", dpi=300)
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(6, 4))
    ax3.hist(pb_list, bins=30)
    ax3.set_xlabel("Primer-bounded reads (K)")
    ax3.set_ylabel("Number of samples")
    ax3.set_title("Primer-bounded reads per sample")
    fig3.tight_layout()
    fig3.savefig(f"{out_prefix}_pb_reads_hist.png", dpi=300)
    plt.close(fig3)

    fig4, ax4 = plt.subplots(figsize=(6, 4))
    ax4.scatter(x_raw, y_gt, s=10)
    ax4.set_xlabel("Raw reads (K)")
    ax4.set_ylabel("Genotyping %")
    ax4.set_ylim(0, 105)
    ax4.set_title("GT% vs raw reads")
    fig4.tight_layout()
    fig4.savefig(f"{out_prefix}_gt_vs_raw.png", dpi=300)
    plt.close(fig4)

    fig5, ax5 = plt.subplots(figsize=(6, 4))
    ax5.scatter(x_pb, y_gt, s=10)
    ax5.set_xlabel("Primer-bounded reads (K)")
    ax5.set_ylabel("Genotyping %")
    ax5.set_ylim(0, 105)
    ax5.set_title("GT% vs primer-bounded reads")
    fig5.tight_layout()
    fig5.savefig(f"{out_prefix}_gt_vs_pb.png", dpi=300)
    plt.close(fig5)

# ---------------------------------------------------------------------
# Accessory assays: sexID / parasite detection
# ---------------------------------------------------------------------


def load_locus_list(path: Optional[Path]) -> list[str]:
    if path is None:
        return []
    vals = []
    with Path(path).open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            vals.append(s.split()[0])
    return vals


def scan_resolved_fastq_counts(resolved_dir: Path) -> tuple[dict[str, Counter], dict[str, int]]:
    sample_locus_counts: dict[str, Counter] = {}
    sample_total_pb: dict[str, int] = {}
    resolved_paths = sorted(p for p in Path(resolved_dir).glob('resolved_*.fastq.gz') if not str(p).endswith('.part'))
    for fpath in PROG(resolved_paths, desc='Scanning resolved reads for accessory assays', unit='sample', total=len(resolved_paths)):
        sample = fpath.name.replace('resolved_', '').replace('.fastq.gz', '')
        cnt = Counter()
        total = 0
        with gzip.open(fpath, 'rt') as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                _ = fh.readline(); _ = fh.readline(); _ = fh.readline()
                if '|locus=' not in h:
                    continue
                loc = h.strip().split('|locus=')[-1]
                cnt[loc] += 1
                total += 1
        sample_locus_counts[sample] = cnt
        sample_total_pb[sample] = total
    return sample_locus_counts, sample_total_pb


def scan_resolved_fastq_counts_and_sequences(resolved_dir: Path) -> tuple[dict[str, Counter], dict[str, int], dict[str, dict[str, Counter]]]:
    """Scan resolved FASTQs and retain per-sample/per-locus sequence counts.

    sample_locus_counts[sample][locus] gives all primer-bounded reads assigned to
    that locus. sample_locus_seq_counts[sample][locus][sequence] gives the
    sequence-resolved counts. The latter is used for sexID markers so putative
    off-target primer-bounded products can be excluded from the normalized p_OT
    calculation.
    """
    sample_locus_counts: dict[str, Counter] = {}
    sample_total_pb: dict[str, int] = {}
    sample_locus_seq_counts: dict[str, dict[str, Counter]] = {}
    resolved_paths = sorted(p for p in Path(resolved_dir).glob('resolved_*.fastq.gz') if not str(p).endswith('.part'))
    for fpath in PROG(resolved_paths, desc='Scanning resolved reads for accessory assays', unit='sample', total=len(resolved_paths)):
        sample = fpath.name.replace('resolved_', '').replace('.fastq.gz', '')
        cnt = Counter()
        seq_cnts: dict[str, Counter] = defaultdict(Counter)
        total = 0
        with gzip.open(fpath, 'rt') as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip().upper()
                _ = fh.readline(); _ = fh.readline()
                if '|locus=' not in h:
                    continue
                loc = h.strip().split('|locus=')[-1]
                cnt[loc] += 1
                seq_cnts[loc][seq] += 1
                total += 1
        sample_locus_counts[sample] = cnt
        sample_locus_seq_counts[sample] = seq_cnts
        sample_total_pb[sample] = total
    return sample_locus_counts, sample_total_pb, sample_locus_seq_counts


def count_pid_matched_reads(seq_counter: Counter, consensus_seq: str, min_pid: float = 0.97) -> tuple[int, int, int, float]:
    """Return (matched_reads, off_target_reads, total_reads, best_pid).

    Reads/sequences matching the sexID consensus at >= min_pid are counted as
    legitimate sex-specific signal. Other primer-bounded products for the same
    marker are retained as off-target/noise counts and excluded from p_OT.
    """
    if not seq_counter:
        return (0, 0, 0, 0.0)
    total = int(sum(seq_counter.values()))
    consensus_seq = (consensus_seq or '').upper()
    if not consensus_seq:
        return (0, total, total, 0.0)
    matched = 0
    best_pid = 0.0
    min_pid = max(0.0, min(1.0, float(min_pid)))
    for seq, n in seq_counter.items():
        seq = (seq or '').upper()
        pid = 1.0 if seq == consensus_seq else _catalog_pairwise_identity(consensus_seq, seq)
        if pid > best_pid:
            best_pid = pid
        if pid >= min_pid:
            matched += int(n)
    off_target = total - matched
    return (int(matched), int(off_target), int(total), float(best_pid))


def marker_specific_sexid_class(
    *,
    reads: int,
    total_pb_reads: int,
    observed_p_ot: float,
    expected_p_ot: float,
    ratio: float,
    sex_system: str = 'XY',
    min_total_pb_reads: int = 100,
    min_marker_reads: int = 2,
    min_marker_pot: float = 0.0005,
    ratio_homo_max: float = 0.10,
    ratio_hetero_min: float = 0.50,
):
    """Classify one sexID marker in one sample from that marker's own signal.

    This is intentionally marker-specific and is used for sexID dashboard plotting.
    It is separate from the multi-marker consensus sex call produced by
    summarize_sexid_markers(). For an XY assay, marker-specific heterogametic
    corresponds to XY/male and marker-specific homogametic corresponds to
    XX/female. For a ZW assay, the labels are reversed accordingly.
    """
    if sex_system == 'XY':
        hetero_label, homo_label = 'male', 'female'
    else:
        hetero_label, homo_label = 'female', 'male'

    reads = int(reads or 0)
    total_pb_reads = int(total_pb_reads or 0)
    observed_p_ot = float(observed_p_ot or 0.0)
    expected_p_ot = float(expected_p_ot or 0.0)
    ratio = float(ratio or 0.0)

    if total_pb_reads < min_total_pb_reads:
        return ('no_call', 'no_call', 'LOW_DATA')
    if expected_p_ot <= 0:
        return ('no_call', 'no_call', 'NO_CALIBRATION')

    marker_pass = (
        reads >= min_marker_reads
        and observed_p_ot >= min_marker_pot
        and ratio >= ratio_hetero_min
    )
    if marker_pass:
        return ('heterogametic', hetero_label, 'PASS')

    # Presence/absence sexID markers need a dead zone between absence and
    # confident presence. Signals below ratio_homo_max are marker-level absence;
    # signals between ratio_homo_max and ratio_hetero_min are weak/ambiguous.
    if ratio < ratio_homo_max:
        return ('homogametic', homo_label, 'PASS')

    return ('ambiguous', 'ambiguous', 'WEAK_SIGNAL')

def summarize_sexid_markers(
    sample_locus_counts: dict[str, Counter],
    sample_total_pb: dict[str, int],
    sex_markers: list[str],
    out_prefix: str,
    sex_system: str = 'XY',
    min_total_pb_reads: int = 100,
    min_marker_reads: int = 2,
    min_marker_pot: float = 0.0005,
    provisional_min_markers: int = 1,
    strong_marker_pot: float = 0.005,
    homo_score_max: float = 0.05,
    hetero_score_min: float = 0.50,
    sample_locus_seq_counts: Optional[dict[str, dict[str, Counter]]] = None,
    sex_consensus_min_pid: float = 0.97,
):
    sex_markers = [m for m in sex_markers if m]
    if not sex_markers:
        return None

    all_samples = sorted(sample_total_pb.keys())
    per_marker_rows = []
    provisional_hetero = []

    # pass 1: observed counts/p_ot and provisional heterogametic set
    for sample in all_samples:
        total_pb = int(sample_total_pb.get(sample, 0))
        cnts = sample_locus_counts.get(sample, Counter())
        present_markers = 0
        strong_hit = False
        for marker in sex_markers:
            reads = int(cnts.get(marker, 0))
            p_ot = (reads / total_pb) if total_pb > 0 else 0.0
            present = int(reads >= min_marker_reads and p_ot >= min_marker_pot)
            if present:
                present_markers += 1
            if p_ot >= strong_marker_pot:
                strong_hit = True
            per_marker_rows.append({
                'sample': sample,
                'marker': marker,
                'reads': reads,
                'total_pb_reads': total_pb,
                'observed_p_ot': p_ot,
                'provisional_present': present,
            })
        if total_pb >= min_total_pb_reads and (present_markers >= provisional_min_markers or strong_hit):
            provisional_hetero.append(sample)

    # Build a consensus sexID sequence for each marker from the provisional
    # heterogametic samples. Expected and observed p_OT below are calculated
    # only from reads with >= sex_consensus_min_pid identity to this consensus.
    # Other primer-bounded products at the same locus are reported as off-target
    # reads and excluded from the sex-specific signal.
    marker_consensus_seq = {}
    marker_consensus_support = {}
    for marker in sex_markers:
        seq_support = Counter()
        if sample_locus_seq_counts is not None:
            for sample in provisional_hetero:
                seq_support.update(sample_locus_seq_counts.get(sample, {}).get(marker, Counter()))
        if seq_support:
            seq, support = seq_support.most_common(1)[0]
            marker_consensus_seq[marker] = seq
            marker_consensus_support[marker] = int(support)

    # calibrate expected p_ot from provisional heterogametic samples
    marker_expected = {}
    for marker in sex_markers:
        vals = []
        for sample in provisional_hetero:
            total_pb = int(sample_total_pb.get(sample, 0))
            if total_pb <= 0:
                continue
            if sample_locus_seq_counts is not None and marker in marker_consensus_seq:
                seq_counter = sample_locus_seq_counts.get(sample, {}).get(marker, Counter())
                reads, _off, _total_marker, _best_pid = count_pid_matched_reads(
                    seq_counter, marker_consensus_seq.get(marker, ''), min_pid=sex_consensus_min_pid
                )
            else:
                reads = int(sample_locus_counts.get(sample, Counter()).get(marker, 0))
            p_ot = reads / total_pb
            if reads >= min_marker_reads and p_ot >= min_marker_pot:
                vals.append(p_ot)
        if vals:
            marker_expected[marker] = {
                'n_heterogametic': len(vals),
                'median_p_ot': statistics.median(vals),
                'mean_p_ot': statistics.mean(vals),
                'min_p_ot': min(vals),
                'max_p_ot': max(vals),
                'consensus_seq': marker_consensus_seq.get(marker, ''),
                'consensus_support': marker_consensus_support.get(marker, 0),
            }

    # write marker calibration
    marker_pot_path = Path(f'{out_prefix}_marker_pot.tsv')
    with marker_pot_path.open('w', newline='') as out:
        w = csv.writer(out, delimiter='	')
        w.writerow(['marker', 'n_heterogametic', 'median_p_ot', 'mean_p_ot', 'min_p_ot', 'max_p_ot', 'consensus_min_pid', 'consensus_support', 'consensus_sequence'])
        for marker in sex_markers:
            d = marker_expected.get(marker)
            if d is None:
                w.writerow([marker, 0, 0, 0, 0, 0, f"{sex_consensus_min_pid:.4f}", marker_consensus_support.get(marker, 0), marker_consensus_seq.get(marker, '')])
            else:
                w.writerow([marker, d['n_heterogametic'], f"{d['median_p_ot']:.8f}", f"{d['mean_p_ot']:.8f}", f"{d['min_p_ot']:.8f}", f"{d['max_p_ot']:.8f}", f"{sex_consensus_min_pid:.4f}", d.get('consensus_support', 0), d.get('consensus_seq', '')])

    # pass 2: score all samples
    score_rows = []
    call_rows = []
    sex_row = {}
    gametic_row = {}
    score_row = {}

    if sex_system == 'XY':
        hetero_label, homo_label = 'male', 'female'
    else:
        hetero_label, homo_label = 'female', 'male'

    usable_markers = [m for m in sex_markers if m in marker_expected and marker_expected[m]['median_p_ot'] > 0]

    for sample in all_samples:
        total_pb = int(sample_total_pb.get(sample, 0))
        cnts = sample_locus_counts.get(sample, Counter())
        marker_passes = 0
        ratios = []
        max_ratio = 0.0
        for marker in sex_markers:
            raw_reads = int(cnts.get(marker, 0))
            matched_reads = raw_reads
            off_target_reads = 0
            best_pid = 0.0
            consensus_seq = marker_expected.get(marker, {}).get('consensus_seq', marker_consensus_seq.get(marker, ''))
            if sample_locus_seq_counts is not None and consensus_seq:
                seq_counter = sample_locus_seq_counts.get(sample, {}).get(marker, Counter())
                matched_reads, off_target_reads, _total_marker_reads, best_pid = count_pid_matched_reads(
                    seq_counter, consensus_seq, min_pid=sex_consensus_min_pid
                )
            reads = int(matched_reads)
            observed_p_ot = (reads / total_pb) if total_pb > 0 else 0.0
            raw_observed_p_ot = (raw_reads / total_pb) if total_pb > 0 else 0.0
            off_target_p_ot = (off_target_reads / total_pb) if total_pb > 0 else 0.0
            expected_p_ot = marker_expected.get(marker, {}).get('median_p_ot', 0.0)
            ratio = 0.0
            if expected_p_ot > 0:
                ratio = observed_p_ot / expected_p_ot
                max_ratio = max(max_ratio, ratio)
                ratios.append(min(1.0, ratio))
            pass_marker = int(expected_p_ot > 0 and reads >= min_marker_reads and observed_p_ot >= min_marker_pot and ratio >= hetero_score_min)
            marker_gametic_class, marker_sex_call, marker_status = marker_specific_sexid_class(
                reads=reads,
                total_pb_reads=total_pb,
                observed_p_ot=observed_p_ot,
                expected_p_ot=expected_p_ot,
                ratio=ratio,
                sex_system=sex_system,
                min_total_pb_reads=min_total_pb_reads,
                min_marker_reads=min_marker_reads,
                min_marker_pot=min_marker_pot,
                ratio_homo_max=homo_score_max,
                ratio_hetero_min=hetero_score_min,
            )
            if pass_marker:
                marker_passes += 1
            score_rows.append({
                'sample': sample,
                'marker': marker,
                'reads': reads,
                'raw_reads': raw_reads,
                'off_target_reads': off_target_reads,
                'total_pb_reads': total_pb,
                'observed_p_ot': observed_p_ot,
                'raw_observed_p_ot': raw_observed_p_ot,
                'off_target_p_ot': off_target_p_ot,
                'expected_p_ot': expected_p_ot,
                'ratio': ratio,
                'best_pid_to_consensus': best_pid,
                'consensus_min_pid': sex_consensus_min_pid,
                'pass': pass_marker,
                'marker_gametic_class': marker_gametic_class,
                'marker_sex_call': marker_sex_call,
                'marker_status': marker_status,
            })
        sex_score = (sum(ratios) / len(ratios)) if ratios else 0.0
        if total_pb < min_total_pb_reads or not usable_markers:
            gametic_class = 'no_call'
            sex_call = 'no_call'
            status = 'LOW_DATA' if total_pb < min_total_pb_reads else 'NO_CALIBRATION'
        elif marker_passes >= max(1, provisional_min_markers) and (sex_score >= hetero_score_min or max_ratio >= 0.50):
            gametic_class = 'heterogametic'
            sex_call = hetero_label
            status = 'PASS'
        elif sex_score <= homo_score_max and marker_passes == 0:
            gametic_class = 'homogametic'
            sex_call = homo_label
            status = 'PASS'
        else:
            gametic_class = 'ambiguous'
            sex_call = 'ambiguous'
            status = 'WEAK_SIGNAL'
        call_rows.append({
            'sample': sample,
            'sex_system': sex_system,
            'sex_call': sex_call,
            'gametic_class': gametic_class,
            'sex_score': sex_score,
            'n_markers_passing': marker_passes,
            'n_markers_tested': len(usable_markers),
            'total_pb_reads': total_pb,
            'status': status,
        })
        sex_row[sample] = sex_call
        gametic_row[sample] = gametic_class
        score_row[sample] = f'{sex_score:.4f}'

    score_path = Path(f'{out_prefix}_sample_marker_scores.tsv')
    with score_path.open('w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['sample', 'marker', 'reads', 'raw_reads', 'off_target_reads', 'total_pb_reads', 'observed_p_ot', 'raw_observed_p_ot', 'off_target_p_ot', 'expected_p_ot', 'ratio', 'best_pid_to_consensus', 'consensus_min_pid', 'pass', 'marker_gametic_class', 'marker_sex_call', 'marker_status'])
        for r in score_rows:
            w.writerow([r['sample'], r['marker'], r['reads'], r.get('raw_reads', r['reads']), r.get('off_target_reads', 0), r['total_pb_reads'], f"{r['observed_p_ot']:.8f}", f"{r.get('raw_observed_p_ot', r['observed_p_ot']):.8f}", f"{r.get('off_target_p_ot', 0.0):.8f}", f"{r['expected_p_ot']:.8f}", f"{r['ratio']:.4f}", f"{r.get('best_pid_to_consensus', 0.0):.4f}", f"{r.get('consensus_min_pid', sex_consensus_min_pid):.4f}", r['pass'], r.get('marker_gametic_class', ''), r.get('marker_sex_call', ''), r.get('marker_status', '')])

    calls_path = Path(f'{out_prefix}_calls.tsv')
    with calls_path.open('w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['sample', 'sex_system', 'sex_call', 'gametic_class', 'sex_score', 'n_markers_passing', 'n_markers_tested', 'total_pb_reads', 'status'])
        for r in call_rows:
            w.writerow([r['sample'], r['sex_system'], r['sex_call'], r['gametic_class'], f"{r['sex_score']:.4f}", r['n_markers_passing'], r['n_markers_tested'], r['total_pb_reads'], r['status']])

    return {
        'sex_row': sex_row,
        'gametic_row': gametic_row,
        'score_row': score_row,
        'call_rows': call_rows,
        'score_rows': score_rows,
        'provisional_heterogametic': provisional_hetero,
        'usable_markers': usable_markers,
        'paths': {
            'marker_pot': marker_pot_path,
            'scores': score_path,
            'calls': calls_path,
        },
    }


def classify_parasite_call(
    reads: int,
    total_pb_reads: int,
    *,
    min_total_pb_reads: int = 100,
    min_reads_detect: int = 2,
    low_pot: float = 0.0005,
    med_pot: float = 0.0020,
    high_pot: float = 0.0100,
    very_high_pot: float = 0.0500,
):
    """
    Return (p_ot, call) for a single parasite marker in one sample.

    Updated logic:
      - LOW_DATA if total primer-bounded depth is too low
      - NEGATIVE if parasite read count is below detection OR p_ot is below low_pot
      - LOW/MED/HIGH/VERY_HIGH otherwise by p_ot band
    """
    total_pb_reads = int(total_pb_reads or 0)
    reads = int(reads or 0)

    p_ot = (reads / total_pb_reads) if total_pb_reads > 0 else 0.0

    if total_pb_reads < min_total_pb_reads:
        return p_ot, "LOW_DATA"
    if reads < min_reads_detect or p_ot < low_pot:
        return p_ot, "NEGATIVE"
    if p_ot < med_pot:
        return p_ot, "LOW"
    if p_ot < high_pot:
        return p_ot, "MED"
    if p_ot < very_high_pot:
        return p_ot, "HIGH"
    return p_ot, "VERY_HIGH"


def summarize_parasite_markers(
    sample_locus_counts,
    sample_total_pb_reads,
    parasite_markers,
    outdir,
    *,
    min_total_pb_reads: int = 100,
    min_reads_detect: int = 2,
    low_pot: float = 0.0005,
    med_pot: float = 0.0020,
    high_pot: float = 0.0100,
    very_high_pot: float = 0.0500,
):
    """
    Build:
      1) parasite_calls.tsv    (long form; one row per sample x marker)
      2) parasite_summary.tsv  (matrix style; one row per sample)

    Expects:
      sample_locus_counts[sample][marker] = integer count of primer-bounded reads
      sample_total_pb_reads[sample] = total primer-bounded reads in that sample
    """
    parasite_markers = [m for m in parasite_markers if m]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not parasite_markers:
        return {"marker_rows": [], "summary_rows": []}

    call_rows = []
    summary_rows = []

    samples = sorted(sample_total_pb_reads.keys())

    for sample in samples:
        total_pb_reads = int(sample_total_pb_reads.get(sample, 0) or 0)
        summary_row = {
            "sample": sample,
            "total_pb_reads": total_pb_reads,
        }

        locus_counts = sample_locus_counts.get(sample, {})

        for marker in parasite_markers:
            reads = int(locus_counts.get(marker, 0) or 0)
            p_ot, call = classify_parasite_call(
                reads,
                total_pb_reads,
                min_total_pb_reads=min_total_pb_reads,
                min_reads_detect=min_reads_detect,
                low_pot=low_pot,
                med_pot=med_pot,
                high_pot=high_pot,
                very_high_pot=very_high_pot,
            )

            call_rows.append({
                "sample": sample,
                "marker": marker,
                "reads": reads,
                "total_pb_reads": total_pb_reads,
                "p_ot": f"{p_ot:.8f}",
                "call": call,
            })
            summary_row[marker] = call

        summary_rows.append(summary_row)

    calls_tsv = outdir / "parasite_calls.tsv"
    with open(calls_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "marker", "reads", "total_pb_reads", "p_ot", "call"])
        for r in call_rows:
            w.writerow([
                r["sample"],
                r["marker"],
                r["reads"],
                r["total_pb_reads"],
                r["p_ot"],
                r["call"],
            ])

    summary_tsv = outdir / "parasite_summary.tsv"
    with open(summary_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "total_pb_reads"] + parasite_markers)
        for r in summary_rows:
            w.writerow([r["sample"], r["total_pb_reads"]] + [r.get(m, "NA") for m in parasite_markers])

    return {
        "marker_rows": call_rows,
        "summary_rows": summary_rows,
        "parasite_calls_tsv": calls_tsv,
        "parasite_summary_tsv": summary_tsv,
    }


def append_assay_rows_to_wide_table(path: Path, row_label: str, row_map: dict[str, str]):
    path = Path(path)
    if not path.exists():
        return
    with path.open() as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh]
    if not rows:
        return
    header = rows[0]
    samples = header[1:]
    rows.append([row_label] + [str(row_map.get(s, '.')) for s in samples])
    with path.open('w', newline='') as out:
        for row in rows:
            out.write('\t'.join(row) + '\n')



def _sexid_consensus_codes(call_rows: list[dict], sex_system: str = "XY") -> dict[str, dict[str, str]]:
    """Map final consensus sex calls to GT-seq-compatible output encodings."""
    sex_system = (sex_system or 'XY').upper()
    if sex_system == 'ZW':
        hetero_label = 'female'
        homo_label = 'male'
        hom_phase = 'Z|Z'
        het_phase = 'Z|W'
    else:
        hetero_label = 'male'
        homo_label = 'female'
        hom_phase = 'X|X'
        het_phase = 'X|Y'

    out = {}
    for r in (call_rows or []):
        sample = r.get('sample')
        sex_call = r.get('sex_call', 'no_call')
        if sex_call == hetero_label:
            out[sample] = {'gt': '201202', 'phase': het_phase, 'n_var': '1'}
        elif sex_call == homo_label:
            out[sample] = {'gt': '201201', 'phase': hom_phase, 'n_var': '1'}
        else:
            out[sample] = {'gt': '000000', 'phase': '.|.', 'n_var': '0'}
    return out


def replace_rows_in_wide_table(path: Path, row_map_by_label: dict[str, dict[str, str]]):
    path = Path(path)
    if not path.exists():
        return
    with path.open() as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    if not rows:
        return
    header = rows[0]
    samples = header[1:]
    new_rows = [header]
    replaced = set()
    for row in rows[1:]:
        label = row[0]
        if label in row_map_by_label:
            mapping = row_map_by_label[label]
            new_rows.append([label] + [str(mapping.get(s, row[i+1] if i+1 < len(row) else '.')) for i, s in enumerate(samples)])
            replaced.add(label)
        else:
            new_rows.append(row)
    for label, mapping in row_map_by_label.items():
        if label not in replaced:
            new_rows.append([label] + [str(mapping.get(s, '.')) for s in samples])
    with path.open('w', newline='') as out:
        for row in new_rows:
            out.write("\t".join(row) + "\n")


def replace_rows_in_phased_longform(path: Path, sex_markers: list[str], call_rows: list[dict], sex_system: str = 'XY'):
    path = Path(path)
    if not path.exists() or not sex_markers:
        return
    sex_markers = set(sex_markers)
    code_map = _sexid_consensus_codes(call_rows, sex_system=sex_system)
    with path.open() as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    if not rows:
        return
    header = rows[0]
    out_rows = [header]
    for row in rows[1:]:
        if len(row) < 7:
            out_rows.append(row)
            continue
        sample, locus = row[0], row[1]
        if locus in sex_markers:
            d = code_map.get(sample, {'gt': '000000', 'phase': '.|.', 'n_var': '0'})
            out_rows.append([sample, locus, d['gt'], d['phase'], d['phase'], d['n_var'], d['n_var']])
        else:
            out_rows.append(row)
    with path.open('w', newline='') as out:
        for row in out_rows:
            out.write("\t".join(row) + "\n")


def harmonize_sexid_outputs(outdir: Path, sex_markers: list[str], sexid_results: dict, sex_system: str = 'XY'):
    if not sex_markers or not sexid_results:
        return
    outdir = Path(outdir)
    code_map = _sexid_consensus_codes(sexid_results.get('call_rows', []), sex_system=sex_system)
    gt_rows = {m: {s: d['gt'] for s, d in code_map.items()} for m in sex_markers}
    phase_rows = {m: {s: d['phase'] for s, d in code_map.items()} for m in sex_markers}

    replace_rows_in_wide_table(outdir / 'microhap_genotypes.tsv', gt_rows)
    replace_rows_in_wide_table(outdir / 'microhap_phasedVARs.tsv', phase_rows)
    replace_rows_in_phased_longform(outdir / 'microhap_phased.tsv', sex_markers, sexid_results.get('call_rows', []), sex_system=sex_system)


def _augment_locus_stats_with_zeros(locus_stats: dict, loci: list[str]) -> dict:
    loci = list(loci or [])
    mean_pct = dict((locus_stats or {}).get("mean_pct", {}) or {})
    sd_pct = dict((locus_stats or {}).get("sd_pct", {}) or {})
    for loc in loci:
        mean_pct.setdefault(loc, 0.0)
        sd_pct.setdefault(loc, 0.0)
    sorted_loci_by_mean = sorted(mean_pct.keys(), key=lambda L: mean_pct.get(L, 0.0))
    return {"sorted_loci_by_mean": sorted_loci_by_mean, "mean_pct": mean_pct, "sd_pct": sd_pct}


def plot_sexid_dashboard(
    locus: str,
    score_rows: list[dict],
    call_rows: list[dict],
    locus_stats: dict,
    sex_system: str,
    outdir: Path,
):
    if not HAVE_MPL:
        return

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Plot marker-specific calls from score_rows, not the multi-marker consensus
    # calls from call_rows. This avoids coloring low-signal samples as
    # heterogametic just because another sexID marker called them that way.
    rows = [r for r in (score_rows or []) if r.get("marker") == locus]
    if not rows:
        return

    cls_to_color = {
        "heterogametic": "C1",
        "homogametic": "C2",
        "ambiguous": "k",
        "no_call": "red",
    }

    fig = plt.figure(figsize=(12, 8.6))
    gs = GridSpec(2, 2, height_ratios=[1.0, 1.4], width_ratios=[1.0, 1.0], hspace=0.35, wspace=0.28)
    fig.subplots_adjust(top=0.88)

    # upper left: expected reads (= expected_p_ot * total_pb_reads) vs observed reads
    ax1 = fig.add_subplot(gs[0, 0])
    # upper right: expected reads vs observed/expected ratio
    ax2 = fig.add_subplot(gs[0, 1])

    grouped = {}
    for r in rows:
        gametic = r.get("marker_gametic_class")
        if not gametic:
            gametic, _sex_call, _status = marker_specific_sexid_class(
                reads=int(r.get("reads", 0)),
                total_pb_reads=int(r.get("total_pb_reads", 0)),
                observed_p_ot=float(r.get("observed_p_ot", 0.0)),
                expected_p_ot=float(r.get("expected_p_ot", 0.0)),
                ratio=float(r.get("ratio", 0.0)),
                sex_system=sex_system,
            )
        x = float(r.get("expected_p_ot", 0.0)) * float(r.get("total_pb_reads", 0))
        y_reads = float(r.get("reads", 0))
        y_ratio = float(r.get("ratio", 0.0))
        grouped.setdefault(gametic, {"x": [], "y_reads": [], "y_ratio": []})
        grouped[gametic]["x"].append(x)
        grouped[gametic]["y_reads"].append(y_reads)
        grouped[gametic]["y_ratio"].append(y_ratio)

    for gametic, vals in grouped.items():
        color = cls_to_color.get(gametic, "0.5")
        label = gametic
        ax1.scatter(vals["x"], vals["y_reads"], s=28, alpha=0.9, c=color, label=label, rasterized=True)
        ax2.scatter(vals["x"], vals["y_ratio"], s=28, alpha=0.9, c=color, label=label, rasterized=True)

    max_x = max([0.0] + [float(r.get("expected_p_ot", 0.0)) * float(r.get("total_pb_reads", 0)) for r in rows])
    max_y_reads = max([0.0] + [float(r.get("reads", 0)) for r in rows])
    lim = max(max_x, max_y_reads, 1.0)
    ax1.plot([0, lim], [0, lim], ls="--", lw=1, color="gray")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel("Expected sexID reads (expected p_OT × total on-target reads)")
    ax1.set_ylabel("Observed sexID marker reads")
    ax1.set_title(f"{locus}: expected vs observed reads")

    ax2.axhline(0.1, ls="--", lw=1, color="gray")
    ax2.axhline(0.5, ls="--", lw=1, color="gray")
    ax2.axhline(1.0, ls="--", lw=1, color="gray")
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel("Expected sexID reads (expected p_OT × total on-target reads)")
    ax2.set_ylabel("Observed p_OT / expected p_OT")
    ax2.set_title(f"{locus}: normalized sexID signal")

    # center annotation
    expected_values = [float(r.get("expected_p_ot", 0.0)) for r in rows if float(r.get("expected_p_ot", 0.0)) > 0]
    median_expected = 0.0
    if expected_values:
        vals = sorted(expected_values)
        mid = len(vals) // 2
        median_expected = vals[mid] if len(vals) % 2 == 1 else (vals[mid - 1] + vals[mid]) / 2.0
    label = f"expected p_OT: {median_expected:.6f} | n samples: {len(rows)}"
    fig.text(
        0.50, 0.545, label,
        ha="center", va="center",
        fontsize=12, color="dimgray", fontweight="semibold",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", boxstyle="round,pad=0.25"),
        transform=fig.transFigure
    )

    # lower left: relative abundance among loci
    padded_stats = _augment_locus_stats_with_zeros(locus_stats, [locus])
    loci_sorted_by_share = padded_stats["sorted_loci_by_mean"]
    locus_pct_mean = padded_stats["mean_pct"]
    locus_pct_sd = padded_stats["sd_pct"]

    ax3 = fig.add_subplot(gs[1, 0])
    means = [locus_pct_mean[L] for L in loci_sorted_by_share]
    sds = [locus_pct_sd[L] for L in loci_sorted_by_share]
    xs3 = list(range(len(loci_sorted_by_share)))
    total = sum(means) if means else 0.0
    if total <= 5.0:
        means = [m * 100.0 for m in means]
        sds = [sd * 100.0 for sd in sds]
    n_loci = max(1, len(loci_sorted_by_share))
    avg_share = 100.0 / n_loci
    highlight_idx = loci_sorted_by_share.index(locus) if locus in loci_sorted_by_share else None

    bars = ax3.bar(xs3, means, color="#b3b3b3", edgecolor="#b3b3b3", linewidth=0.0, zorder=1)
    if sds:
        ax3.errorbar(xs3, means, yerr=sds, fmt="none", ecolor="#777777", elinewidth=0.8, capthick=0.8, capsize=1, zorder=2)
    if highlight_idx is not None:
        bars[highlight_idx].set_facecolor("red")
        bars[highlight_idx].set_edgecolor("red")
        bars[highlight_idx].set_linewidth(1.0)
        bars[highlight_idx].set_zorder(4)
        ax3.errorbar(highlight_idx, means[highlight_idx], yerr=sds[highlight_idx], fmt="none", ecolor="black", elinewidth=1.2, capthick=1.2, capsize=2, zorder=5)
    ax3.axhline(avg_share, color="tab:red", lw=1, zorder=6)
    ymax = max((m + sd) for m, sd in zip(means or [0.0], sds or [0.0]))
    ax3.set_ylim(0, max(ymax * 1.15, avg_share * 1.8, 2.0))
    ax3.set_xlim(-0.5, len(xs3) - 0.5 + 1.0)
    ax3.set_ylabel("Avg % reads per locus (±1 SD)")
    ax3.set_xlabel("Loci (sorted by mean %, least → most)")
    ax3.set_title("Read distribution among loci")
    if len(xs3) > 20:
        ax3.set_xticks([])

    # lower right: marker-specific call counts
    ax4 = fig.add_subplot(gs[1, 1])
    labels = ["heterogametic", "homogametic", "ambiguous", "no_call"]
    counts = {lab: 0 for lab in labels}
    for r in rows:
        gc = r.get("marker_gametic_class")
        if not gc:
            gc, _sex_call, _status = marker_specific_sexid_class(
                reads=int(r.get("reads", 0)),
                total_pb_reads=int(r.get("total_pb_reads", 0)),
                observed_p_ot=float(r.get("observed_p_ot", 0.0)),
                expected_p_ot=float(r.get("expected_p_ot", 0.0)),
                ratio=float(r.get("ratio", 0.0)),
                sex_system=sex_system,
            )
        if gc in counts:
            counts[gc] += 1
    display_labels = []
    if sex_system == "XY":
        mapping = {"heterogametic": "XY", "homogametic": "XX", "ambiguous": "ambiguous", "no_call": "no_call"}
    else:
        mapping = {"heterogametic": "ZW", "homogametic": "ZZ", "ambiguous": "ambiguous", "no_call": "no_call"}
    display_labels = [mapping[x] for x in labels]
    ax4.bar(range(len(labels)), [counts[x] for x in labels])
    ax4.set_xticks(range(len(labels)))
    ax4.set_xticklabels(display_labels)
    ax4.set_ylabel("Number of samples")
    ax4.set_title(f"{locus}: marker-specific sexID class counts")

    handles, labels_ = ax1.get_legend_handles_labels()
    if handles:
        ax1.legend(handles, labels_, loc="best", frameon=False)
    outpath = outdir / f"{locus}_sexid_dashboard.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def plot_parasite_dashboard(
    marker: str,
    marker_rows: list[dict],
    locus_stats: dict,
    outdir: Path,
    min_total_pb_reads: int = 100,
    min_reads_detect: int = 2,
    low_pot: float = 0.0005,
    med_pot: float = 0.0020,
    high_pot: float = 0.0100,
    very_high_pot: float = 0.0500,
):
    if not HAVE_MPL:
        return

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in (marker_rows or []) if r.get('marker') == marker]
    if not rows:
        return

    cls_to_color = {
        'LOW_DATA': '#9e9e9e',
        'NEGATIVE': '#66bb6a',
        'LOW': '#2e7d32',
        'MED': '#fb8c00',
        'HIGH': '#e53935',
        'VERY_HIGH': '#8e24aa',
    }

    fig = plt.figure(figsize=(12, 8.6))
    gs = GridSpec(
        2, 2,
        height_ratios=[1.0, 1.4],
        width_ratios=[1.0, 1.0],
        hspace=0.35,
        wspace=0.28
    )
    fig.subplots_adjust(top=0.88)

    # -------------------------
    # upper left: 1% of total on-target reads vs observed parasite reads
    # -------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    grouped = {}
    for r in rows:
        call = r.get('call', 'NEGATIVE')
        x = float(r.get('total_pb_reads', 0)) * 0.01
        y = float(r.get('reads', 0))
        grouped.setdefault(call, {'x': [], 'y': []})
        grouped[call]['x'].append(x)
        grouped[call]['y'].append(y)

    for call, vals in grouped.items():
        color = cls_to_color.get(call, '0.5')
        ax1.scatter(
            vals['x'], vals['y'],
            s=28, alpha=0.9, c=color,
            label=call, rasterized=True
        )

    max_x1 = max([1.0] + [float(r.get('total_pb_reads', 0)) * 0.01 for r in rows])
    max_y1 = max([1.0] + [float(r.get('reads', 0)) for r in rows])
    lim = max(max_x1, max_y1)
    ax1.plot([0, lim], [0, lim], ls='--', lw=1, color='gray')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel('Total primer-bounded reads × 0.01')
    ax1.set_ylabel('Primer-bounded parasite reads')
    ax1.set_title(f'{marker}: parasite reads vs 1% on-target depth')

    # -------------------------
    # upper right: parasite reads vs p_ot (log x, log y)
    # -------------------------
    ax2 = fig.add_subplot(gs[0, 1])

    raw_reads = [float(r.get('reads', 0)) for r in rows]
    plot_reads = [max(v, 1.0) for v in raw_reads]  # log-safe x
    xmax = max([1.0] + plot_reads)

    pot_floor = min(low_pot / 10.0, 1e-6)
    raw_pots = [float(r.get('p_ot', 0.0)) for r in rows]
    plot_pots = [max(v, pot_floor) for v in raw_pots]
    ymax = max(
        very_high_pot * 1.2,
        max([very_high_pot] + plot_pots) * 1.05
    )

    # Negative strip: driven by reads < min_reads_detect
    neg_left = 1.0
    neg_right = max(float(min_reads_detect), 1.0)
    if neg_right > neg_left:
        ax2.axvspan(neg_left, neg_right, color='#c8e6c9', alpha=0.35, zorder=0)

    # Horizontal qualitative bands for p_ot
    bands = [
        (low_pot, med_pot, '#a5d6a7', 'LOW'),
        (med_pot, high_pot, '#ffcc80', 'MED'),
        (high_pot, very_high_pot, '#ef9a9a', 'HIGH'),
        (very_high_pot, ymax, '#ce93d8', 'VERY_HIGH'),
    ]
    for y0, y1, color, label in bands:
        ax2.axhspan(y0, y1, color=color, alpha=0.45, zorder=0)
        ax2.axhline(y1, ls='--', lw=1, color='gray', zorder=1)
        ymid = math.sqrt(y0 * y1)
        ax2.text(
            xmax * 0.97,
            ymid,
            label,
            ha='right',
            va='center',
            fontsize=11,
            color='black'
        )

    # Add the NEGATIVE label in the same right-side label column
    neg_ymid = math.sqrt(pot_floor * low_pot)
    ax2.text(
        xmax * 0.97,
        neg_ymid,
        'NEGATIVE',
        ha='right',
        va='center',
        fontsize=11,
        color='black'
    )

    ax2.axhline(low_pot, ls='--', lw=1, color='gray', zorder=1)

    # Scatter points: reconcile colors with the negative strip rule.
    # Any sample with reads < min_reads_detect is plotted as NEGATIVE here.
    grouped2 = {}
    for r in rows:
        reads = float(r.get('reads', 0))
        pot = float(r.get('p_ot', 0.0))
        call = r.get('call', 'NEGATIVE')

        plot_call = 'NEGATIVE' if reads < min_reads_detect else call
        x = max(reads, 1.0)
        y = max(pot, pot_floor)

        grouped2.setdefault(plot_call, {'x': [], 'y': []})
        grouped2[plot_call]['x'].append(x)
        grouped2[plot_call]['y'].append(y)

    for call, vals in grouped2.items():
        color = cls_to_color.get(call, '0.5')
        ax2.scatter(
            vals['x'],
            vals['y'],
            s=28,
            alpha=0.9,
            c=color,
            rasterized=True,
            zorder=3
        )

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1.0, xmax * 1.05)
    ax2.set_ylim(pot_floor, ymax)
    ax2.set_xlabel('Primer-bounded parasite reads')
    ax2.set_ylabel('p_OT (log scale)')
    ax2.set_title(f'{marker}: parasite read count vs p_OT')

    positive = sum(1 for r in rows if r.get('call') not in ('NEGATIVE', 'LOW_DATA'))
    fig.text(
        0.50, 0.545,
        f'min_reads={min_reads_detect} | low={low_pot:.4g} med={med_pot:.4g} high={high_pot:.4g} very_high={very_high_pot:.4g} | positives={positive}/{len(rows)}',
        ha='center',
        va='center',
        fontsize=11,
        color='dimgray',
        fontweight='semibold',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', boxstyle='round,pad=0.25'),
        transform=fig.transFigure
    )

    # -------------------------
    # lower left: read distribution among loci
    # -------------------------
    ax3 = fig.add_subplot(gs[1, 0])

    loci_sorted = locus_stats.get("sorted_loci_by_mean", [])
    locus_pct_mean = locus_stats.get("mean_pct", {})
    locus_pct_sd = locus_stats.get("sd_pct", {})

    means = [locus_pct_mean[L] for L in loci_sorted]
    sds = [locus_pct_sd[L] for L in loci_sorted]
    xs3 = list(range(len(loci_sorted)))

    if means:
        total = sum(means)
        if total <= 5.0:
            means = [m * 100.0 for m in means]
            sds = [sd * 100.0 for sd in sds]

    avg_share = 100.0 / max(1, len(loci_sorted))

    bars = ax3.bar(
        xs3,
        means,
        color="#b3b3b3",
        edgecolor="#b3b3b3",
        linewidth=0.0,
        zorder=1
    )

    if sds:
        ax3.errorbar(
            xs3,
            means,
            yerr=sds,
            fmt="none",
            ecolor="#777777",
            elinewidth=0.8,
            capthick=0.8,
            capsize=1,
            zorder=2
        )

    marker_idx = None
    if marker in loci_sorted:
        marker_idx = loci_sorted.index(marker)
        bars[marker_idx].set_facecolor("red")
        bars[marker_idx].set_edgecolor("red")
        bars[marker_idx].set_linewidth(1.0)
        bars[marker_idx].set_zorder(4)

        if sds:
            ax3.errorbar(
                marker_idx,
                means[marker_idx],
                yerr=sds[marker_idx],
                fmt="none",
                ecolor="black",
                elinewidth=1.2,
                capthick=1.2,
                capsize=2,
                zorder=5
            )

    ax3.axhline(avg_share, color="tab:red", lw=1, zorder=6)

    ymax3 = 1.0
    if means:
        ymax3 = max((m + sd) for m, sd in zip(means, sds)) if sds else max(means)
    ax3.set_ylim(0, max(ymax3 * 1.15, avg_share * 1.8, 2.0))
    ax3.set_xlim(-0.5, len(xs3) - 0.5 + 1.0)
    ax3.set_ylabel("Avg % reads per locus (±1 SD)")
    ax3.set_xlabel("Loci (sorted by mean %, least → most)")
    ax3.set_title("Read distribution among loci")
    if len(xs3) > 20:
        ax3.set_xticks([])

    # -------------------------
    # lower right: parasite class counts
    # -------------------------
    ax4 = fig.add_subplot(gs[1, 1])

    class_order = ['NEGATIVE', 'LOW', 'MED', 'HIGH', 'VERY_HIGH', 'LOW_DATA']
    class_counts = {k: 0 for k in class_order}
    for r in rows:
        call = r.get('call', 'NEGATIVE')
        if call not in class_counts:
            class_counts[call] = 0
        class_counts[call] += 1

    x4 = list(range(len(class_order)))
    y4 = [class_counts.get(k, 0) for k in class_order]
    c4 = [cls_to_color.get(k, '0.5') for k in class_order]

    ax4.bar(x4, y4, color=c4)
    ax4.set_xticks(x4)
    ax4.set_xticklabels(class_order, rotation=20)
    ax4.set_ylabel('Number of samples')
    ax4.set_title(f'{marker}: parasite class counts')

    handles, labels_ = ax1.get_legend_handles_labels()
    if handles:
        ax1.legend(handles, labels_, loc="best", frameon=False)

    outpath = outdir / f"{marker}_parasite_dashboard.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)
    
def export_single_snp_outputs(
    catalog_csv: Path,
    microhap_genotypes_tsv: Path,
    microhap_phased_def_tsv: Path,
    out_first_tsv: Path,
    out_highest_maf_tsv: Path,
    out_metadata_tsv: Path,
    mode: str = "both",
    exclude_loci=None,
):
    """
    Export optional single-SNP genotype matrices derived from catalog-based microhap outputs.

    Modes
    -----
    first
        For each locus, export the leftmost valid biallelic SNP.
    highest_maf
        For each locus, export the valid biallelic SNP with highest minor allele frequency
        in the current dataset.
    both
        Write both outputs.
    none
        Do nothing.

    Inputs
    ------
    catalog_csv
        microhap catalog CSV containing locus, sequence, allele_code
    microhap_genotypes_tsv
        Wide genotype matrix from call_from_resolved_fastqs() with 6-digit allele-code genotypes
    microhap_phased_def_tsv
        Per-locus phased variant definitions written by call_from_resolved_fastqs()

    Outputs
    -------
    out_first_tsv
        Wide SNP matrix using first valid SNP per locus
    out_highest_maf_tsv
        Wide SNP matrix using highest-MAF valid SNP per locus
    out_metadata_tsv
        Per-locus metadata describing valid SNPs and selected positions
    """
    import csv
    from pathlib import Path
    from collections import defaultdict

    mode = (mode or "none").lower()
    if mode == "none":
        return
    if mode not in {"first", "highest_maf", "both"}:
        raise ValueError(f"Unsupported single-SNP mode: {mode}")

    exclude_loci = set(exclude_loci or [])

    catalog_csv = Path(catalog_csv)
    microhap_genotypes_tsv = Path(microhap_genotypes_tsv)
    microhap_phased_def_tsv = Path(microhap_phased_def_tsv)
    out_first_tsv = Path(out_first_tsv)
    out_highest_maf_tsv = Path(out_highest_maf_tsv)
    out_metadata_tsv = Path(out_metadata_tsv)

    def _get_code_hap(pm, code):
        if code in pm["code_to_hap_snp"]:
            return pm["code_to_hap_snp"][code]
        return pm["code_to_hap_snp"].get(str(code), None)

    def _is_called_gt(gt: str) -> bool:
        if not gt:
            return False
        gt = gt.strip()
        return (gt != "000000") and (len(gt) == 6) and gt.isdigit()

    def _split_gt_codes(gt: str):
        return int(gt[:3]), int(gt[3:])

    def _canonical_bases(chars):
        return sorted({c for c in chars if c in {"A", "C", "G", "T"}})

    def _format_snp_gt_from_codes(gt: str, pm, snp_idx: int) -> str:
        if not _is_called_gt(gt):
            return "./."
        c1, c2 = _split_gt_codes(gt)
        h1 = _get_code_hap(pm, c1)
        h2 = _get_code_hap(pm, c2)
        if not h1 or not h2 or h1 == "." or h2 == ".":
            return "./."
        if snp_idx >= len(h1) or snp_idx >= len(h2):
            return "./."
        b1 = h1[snp_idx]
        b2 = h2[snp_idx]
        if b1 not in {"A", "C", "G", "T"} or b2 not in {"A", "C", "G", "T"}:
            return "./."
        a, b = sorted([b1, b2])
        return f"{a}/{b}"

    def _compute_maf_for_snp_idx(sample_gt_map, pm, snp_idx: int):
        counts = defaultdict(int)
        total = 0
        for _sample, gt in sample_gt_map.items():
            if not _is_called_gt(gt):
                continue
            c1, c2 = _split_gt_codes(gt)
            h1 = _get_code_hap(pm, c1)
            h2 = _get_code_hap(pm, c2)
            if not h1 or not h2 or h1 == "." or h2 == ".":
                continue
            if snp_idx >= len(h1) or snp_idx >= len(h2):
                continue
            b1 = h1[snp_idx]
            b2 = h2[snp_idx]
            if b1 not in {"A", "C", "G", "T"} or b2 not in {"A", "C", "G", "T"}:
                continue
            counts[b1] += 1
            counts[b2] += 1
            total += 2

        if total == 0:
            return None, 0, {}
        freqs = {b: counts[b] / float(total) for b in sorted(counts)}
        if len(freqs) != 2:
            return None, total, freqs
        maf = min(freqs.values())
        return maf, total, freqs

    # ------------------------------------------------------------------
    # Load catalog and rebuild phased SNP maps from catalog alleles
    # ------------------------------------------------------------------
    catalog = defaultdict(dict)   # locus -> seq -> allele_code
    with catalog_csv.open() as fh:
        rdr = csv.DictReader(fh)
        cols = {k.lower(): k for k in (rdr.fieldnames or [])}
        need = ("locus", "sequence", "allele_code")
        if any(k not in cols for k in need):
            raise ValueError(f"Catalog CSV missing required columns: {need}")
        for row in rdr:
            loc = row[cols["locus"]]
            seq = row[cols["sequence"]]
            code = int(row[cols["allele_code"]])
            catalog[loc][seq] = code

    phased_maps, _def_rows = build_phased_definitions_from_catalog(catalog)

    # ------------------------------------------------------------------
    # Load phased SNP definitions from file so metadata uses written defs
    # ------------------------------------------------------------------
    snp_defs_by_locus = defaultdict(list)
    if microhap_phased_def_tsv.exists():
        with microhap_phased_def_tsv.open() as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                loc = row["locus"]
                if loc in exclude_loci:
                    continue
                if row["type"] != "SNP":
                    continue
                snp_defs_by_locus[loc].append({
                    "var_index": int(row["var_index"]),
                    "ref_pos": int(row["ref_pos"]),
                    "ref_base": row["ref_base"],
                    "alts": row["alts"],
                })

        for loc in snp_defs_by_locus:
            snp_defs_by_locus[loc].sort(key=lambda r: r["var_index"])
    else:
        # Fallback: build lightweight rows from the in-memory def_rows
        for row in _def_rows:
            loc = row["locus"]
            if loc in exclude_loci:
                continue
            if row["type"] != "SNP":
                continue
            snp_defs_by_locus[loc].append({
                "var_index": int(row["var_index"]),
                "ref_pos": int(row["ref_pos"]),
                "ref_base": row["ref_base"],
                "alts": row["alts"],
            })
        for loc in snp_defs_by_locus:
            snp_defs_by_locus[loc].sort(key=lambda r: r["var_index"])

    # ------------------------------------------------------------------
    # Load wide microhap genotype matrix
    # ------------------------------------------------------------------
    sample_names = []
    per_locus_sample_gt = {}
    with microhap_genotypes_tsv.open() as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        sample_names = header[1:]
        for row in rdr:
            if not row:
                continue
            loc = row[0]
            if loc in exclude_loci:
                continue
            per_locus_sample_gt[loc] = dict(zip(sample_names, row[1:]))

    # ------------------------------------------------------------------
    # Per-locus SNP selection
    # ------------------------------------------------------------------
    metadata_rows = []
    first_matrix = {}
    maf_matrix = {}

    loci = sorted([loc for loc in per_locus_sample_gt.keys() if loc not in exclude_loci])

    for loc in loci:
        sample_gt_map = per_locus_sample_gt.get(loc, {})
        pm = phased_maps.get(loc)
        snp_defs = snp_defs_by_locus.get(loc, [])

        row_meta = {
            "locus": loc,
            "n_samples": len(sample_names),
            "n_called_microhap": sum(1 for gt in sample_gt_map.values() if _is_called_gt(gt)),
            "n_catalog_alleles": len(catalog.get(loc, {})),
            "n_snp_defs": len(snp_defs),
            "n_valid_biallelic_snps": 0,
            "first_snp_idx": "",
            "first_var_index": "",
            "first_ref_pos": "",
            "first_ref_base": "",
            "first_alts": "",
            "highest_maf_snp_idx": "",
            "highest_maf_var_index": "",
            "highest_maf_ref_pos": "",
            "highest_maf_ref_base": "",
            "highest_maf_alts": "",
            "highest_maf_value": "",
            "highest_maf_allele_counts": "",
            "note": "",
        }

        if pm is None:
            row_meta["note"] = "locus_missing_from_phased_maps"
            metadata_rows.append(row_meta)
            continue

        # Collect per-allele SNP hap strings
        code_to_hap = {}
        for _seq, code in catalog.get(loc, {}).items():
            hap = _get_code_hap(pm, code)
            if hap and hap != ".":
                code_to_hap[int(code)] = hap

        if not code_to_hap:
            row_meta["note"] = "no_snp_haps_for_locus"
            metadata_rows.append(row_meta)
            continue

        # Number of SNP positions implied by phased SNP haps
        hap_lens = sorted({len(h) for h in code_to_hap.values()})
        if len(hap_lens) != 1:
            row_meta["note"] = "inconsistent_snp_hap_lengths"
            metadata_rows.append(row_meta)
            continue

        n_snp = hap_lens[0]
        if n_snp == 0:
            row_meta["note"] = "no_snp_positions"
            metadata_rows.append(row_meta)
            continue

        # Identify valid biallelic SNP indices
        valid_idx = []
        for i in range(n_snp):
            bases = _canonical_bases(h[i] for h in code_to_hap.values() if i < len(h))
            if len(bases) == 2:
                valid_idx.append(i)

        row_meta["n_valid_biallelic_snps"] = len(valid_idx)

        if not valid_idx:
            row_meta["note"] = "no_valid_biallelic_snp"
            metadata_rows.append(row_meta)
            continue

        # FIRST SNP
        first_idx = valid_idx[0]
        if first_idx < len(snp_defs):
            d = snp_defs[first_idx]
            row_meta["first_snp_idx"] = first_idx
            row_meta["first_var_index"] = d["var_index"]
            row_meta["first_ref_pos"] = d["ref_pos"]
            row_meta["first_ref_base"] = d["ref_base"]
            row_meta["first_alts"] = d["alts"]
        else:
            row_meta["first_snp_idx"] = first_idx
            row_meta["note"] = (row_meta["note"] + ";first_idx_out_of_range").strip(";")

        # HIGHEST MAF SNP
        best_idx = None
        best_maf = None
        best_counts = None
        for i in valid_idx:
            maf, total_alleles, freqs = _compute_maf_for_snp_idx(sample_gt_map, pm, i)
            if maf is None:
                continue
            if (best_maf is None) or (maf > best_maf) or (maf == best_maf and i < best_idx):
                best_maf = maf
                best_idx = i
                if freqs:
                    best_counts = ",".join(f"{b}:{freqs[b]:.4f}" for b in sorted(freqs))
                else:
                    best_counts = ""

        if best_idx is not None:
            if best_idx < len(snp_defs):
                d = snp_defs[best_idx]
                row_meta["highest_maf_snp_idx"] = best_idx
                row_meta["highest_maf_var_index"] = d["var_index"]
                row_meta["highest_maf_ref_pos"] = d["ref_pos"]
                row_meta["highest_maf_ref_base"] = d["ref_base"]
                row_meta["highest_maf_alts"] = d["alts"]
            else:
                row_meta["highest_maf_snp_idx"] = best_idx
                row_meta["note"] = (row_meta["note"] + ";highest_maf_idx_out_of_range").strip(";")
            row_meta["highest_maf_value"] = f"{best_maf:.6f}"
            row_meta["highest_maf_allele_counts"] = best_counts or ""

        # Build output rows
        if mode in {"first", "both"}:
            first_matrix[loc] = {
                s: _format_snp_gt_from_codes(sample_gt_map.get(s, "000000"), pm, first_idx)
                for s in sample_names
            }

        if mode in {"highest_maf", "both"}:
            if best_idx is None:
                maf_matrix[loc] = {s: "./." for s in sample_names}
            else:
                maf_matrix[loc] = {
                    s: _format_snp_gt_from_codes(sample_gt_map.get(s, "000000"), pm, best_idx)
                    for s in sample_names
                }

        metadata_rows.append(row_meta)

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    if mode in {"first", "both"}:
        with out_first_tsv.open("w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["locus"] + sample_names)
            for loc in loci:
                row = first_matrix.get(loc, {s: "./." for s in sample_names})
                w.writerow([loc] + [row.get(s, "./.") for s in sample_names])

    if mode in {"highest_maf", "both"}:
        with out_highest_maf_tsv.open("w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["locus"] + sample_names)
            for loc in loci:
                row = maf_matrix.get(loc, {s: "./." for s in sample_names})
                w.writerow([loc] + [row.get(s, "./.") for s in sample_names])

    with out_metadata_tsv.open("w", newline="") as out:
        fieldnames = [
            "locus",
            "n_samples",
            "n_called_microhap",
            "n_catalog_alleles",
            "n_snp_defs",
            "n_valid_biallelic_snps",
            "first_snp_idx",
            "first_var_index",
            "first_ref_pos",
            "first_ref_base",
            "first_alts",
            "highest_maf_snp_idx",
            "highest_maf_var_index",
            "highest_maf_ref_pos",
            "highest_maf_ref_base",
            "highest_maf_alts",
            "highest_maf_value",
            "highest_maf_allele_counts",
            "note",
        ]
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in metadata_rows:
            w.writerow(row)

# ---------------------------------------------------------------------
# CLI + main
# ---------------------------------------------------------------------



# ----------------------------------------------------------------------
# Reference-anchored VCF export from catalog microhap calls
# ----------------------------------------------------------------------

def _vcf_check_bwa_and_index(reference_fasta: Path, bwa_path: str = "bwa", build_index: bool = False):
    """Validate BWA availability and BWA index files for optional VCF export."""
    import shutil
    import subprocess

    bwa_exe = shutil.which(str(bwa_path)) if os.sep not in str(bwa_path) else str(bwa_path)
    if not bwa_exe or not Path(bwa_exe).exists():
        raise RuntimeError(
            "VCF export requested, but BWA was not found. Install BWA or pass --bwa-path /path/to/bwa."
        )

    reference_fasta = Path(reference_fasta)
    if not reference_fasta.exists():
        raise RuntimeError(f"VCF export requested, but reference FASTA does not exist: {reference_fasta}")

    required = [reference_fasta.with_suffix(reference_fasta.suffix + ext) for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    # BWA index files are normally named genome.fa.amb etc. Path.with_suffix would make genome.fa.amb.
    required = [Path(str(reference_fasta) + ext) for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    if not all(p.exists() for p in required):
        if build_index:
            print(f"VCF export: BWA index not found; running: {bwa_exe} index {reference_fasta}", file=sys.stderr)
            subprocess.run([bwa_exe, "index", str(reference_fasta)], check=True)
        else:
            missing = ", ".join(str(p.name) for p in required if not p.exists())
            raise RuntimeError(
                "VCF export requested, but BWA index files are missing for the reference FASTA. "
                f"Missing: {missing}. Run 'bwa index <reference.fa>' or use --vcf-build-bwa-index."
            )
    return bwa_exe


def _vcf_write_catalog_fasta(catalog_by_locus: dict, fasta_path: Path):
    """Write catalog alleles as FASTA records named locus|allele_code."""
    with Path(fasta_path).open("w") as out:
        for locus in sorted(catalog_by_locus):
            for code in sorted(catalog_by_locus[locus]):
                seq = catalog_by_locus[locus][code]
                out.write(f">{locus}|{code}\n{seq}\n")


def _vcf_load_catalog_by_locus(catalog_csv: Path) -> dict[str, dict[int, str]]:
    catalog_by_locus: dict[str, dict[int, str]] = defaultdict(dict)
    with Path(catalog_csv).open() as fh:
        rdr = csv.DictReader(fh)
        cols = {c.lower(): c for c in (rdr.fieldnames or [])}
        for req in ("locus", "allele_code", "sequence"):
            if req not in cols:
                raise RuntimeError(f"VCF export: catalog CSV missing required column '{req}'")
        for row in rdr:
            locus = row[cols["locus"]]
            code = int(row[cols["allele_code"]])
            seq = row[cols["sequence"]].upper()
            catalog_by_locus[locus][code] = seq
    return catalog_by_locus


def _vcf_parse_cigar(cigar: str):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar or "")]


def _vcf_cigar_ref_query_lengths(cigar: str) -> tuple[int, int, int]:
    ref_len = query_aln_len = matches = 0
    for n, op in _vcf_parse_cigar(cigar):
        if op in "M=X":
            ref_len += n
            query_aln_len += n
            matches += n if op == "=" else 0
        elif op in "DN":
            ref_len += n
        elif op == "I":
            query_aln_len += n
        elif op in "SHP":
            pass
    return ref_len, query_aln_len, matches


class _VcfIndexedFasta:
    """Tiny faidx reader for fetching reference bases without requiring samtools/pysam."""

    def __init__(self, fasta_path: Path):
        self.path = Path(fasta_path)
        self.fai = Path(str(self.path) + ".fai")
        if not self.fai.exists():
            self._build_fai()
        self.index = {}
        with self.fai.open() as fh:
            for line in fh:
                if not line.strip():
                    continue
                name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(line_bases), int(line_width))
        self._fh = open(self.path, "rb")

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass

    def _build_fai(self):
        with self.path.open("rb") as fh, self.fai.open("w") as out:
            name = None
            length = 0
            seq_offset = None
            line_bases = None
            line_width = None
            while True:
                pos = fh.tell()
                line = fh.readline()
                if not line:
                    if name is not None:
                        out.write(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}\n")
                    break
                if line.startswith(b">"):
                    if name is not None:
                        out.write(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}\n")
                    name = line[1:].decode(errors="ignore").strip().split()[0]
                    length = 0
                    seq_offset = None
                    line_bases = None
                    line_width = None
                else:
                    raw = line.rstrip(b"\r\n")
                    if seq_offset is None:
                        seq_offset = pos
                        line_bases = len(raw)
                        line_width = len(line)
                    length += len(raw)

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Fetch 1-based inclusive interval."""
        if contig not in self.index:
            raise KeyError(f"contig not found in reference FASTA index: {contig}")
        length, offset, line_bases, line_width = self.index[contig]
        start = max(1, int(start))
        end = min(length, int(end))
        if end < start:
            return ""
        out = []
        pos = start
        while pos <= end:
            line_idx = (pos - 1) // line_bases
            in_line = (pos - 1) % line_bases
            n = min(end - pos + 1, line_bases - in_line)
            byte_offset = offset + line_idx * line_width + in_line
            self._fh.seek(byte_offset)
            out.append(self._fh.read(n).decode().upper())
            pos += n
        return "".join(out)

    def base(self, contig: str, pos: int) -> str:
        return self.fetch(contig, pos, pos).upper()


def _vcf_run_bwa_mem(catalog_fasta: Path, reference_fasta: Path, sam_path: Path, bwa_exe: str, threads: int = 1):
    import subprocess
    threads = max(1, int(threads or 1))
    cmd = [bwa_exe, "mem", "-t", str(threads), str(reference_fasta), str(catalog_fasta)]
    print("VCF export: running " + " ".join(cmd), file=sys.stderr)
    with Path(sam_path).open("w") as out:
        subprocess.run(cmd, stdout=out, stderr=sys.stderr, check=True)


def _vcf_sam_to_alignments(sam_path: Path, catalog_by_locus: dict, min_mapq: int, min_aln_frac: float, min_pid: float):
    """Parse primary BWA SAM alignments into per-allele placement rows."""
    best_by_allele = {}
    seq_lookup = {(loc, int(code)): seq for loc, d in catalog_by_locus.items() for code, seq in d.items()}
    with Path(sam_path).open() as fh:
        for line in fh:
            if not line or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            qname, flag_s, rname, pos_s, mapq_s, cigar = parts[:6]
            flag = int(flag_s)
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            if "|" not in qname:
                continue
            locus, code_s = qname.rsplit("|", 1)
            try:
                code = int(code_s)
            except Exception:
                continue
            seq = seq_lookup.get((locus, code), parts[9].upper())
            ref_len, query_aln_len, _matches_eq = _vcf_cigar_ref_query_lengths(cigar)
            nm = 0
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    try:
                        nm = int(tag.split(":")[-1])
                    except Exception:
                        nm = 0
                    break
            aln_bases = max(1, query_aln_len)
            pid = max(0.0, 1.0 - (nm / float(aln_bases)))
            aln_frac = query_aln_len / float(max(1, len(seq)))
            mapq = int(mapq_s)
            strand = "-" if (flag & 0x10) else "+"
            start = int(pos_s)
            end = start + ref_len - 1
            status = "PASS"
            reasons = []
            if mapq < min_mapq:
                reasons.append("low_mapq")
            if aln_frac < min_aln_frac:
                reasons.append("low_aln_frac")
            if pid < min_pid:
                reasons.append("low_pid")
            if reasons:
                status = ";".join(reasons)
            row = {
                "locus": locus,
                "allele_code": code,
                "contig": rname,
                "start": start,
                "end": end,
                "strand": strand,
                "mapq": mapq,
                "cigar": cigar,
                "pid": pid,
                "aln_frac": aln_frac,
                "flag": flag,
                "seq": seq,
                "status": status,
            }
            key = (locus, code)
            prev = best_by_allele.get(key)
            if prev is None or (row["mapq"], row["pid"], row["aln_frac"]) > (prev["mapq"], prev["pid"], prev["aln_frac"]):
                best_by_allele[key] = row
    return list(best_by_allele.values())


def _vcf_write_placements(rows: list[dict], out_tsv: Path):
    cols = ["locus", "allele_code", "contig", "start", "end", "strand", "mapq", "pid", "aln_frac", "cigar", "status"]
    with Path(out_tsv).open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for r in sorted(rows, key=lambda x: (x.get("locus", ""), int(x.get("allele_code", 0)))):
            w.writerow([r.get(c, "") if c not in {"pid", "aln_frac"} else f"{float(r.get(c, 0.0)):.6f}" for c in cols])


def _vcf_project_alignment(row: dict, ref: _VcfIndexedFasta):
    """Return per-base projection and indel events for one allele alignment."""
    contig = row["contig"]
    seq = row["seq"].upper()
    if row["strand"] == "-":
        seq = rc(seq)
    ref_pos = int(row["start"])
    qpos = 0
    bases = {}       # ref_pos -> allele base
    insertions = {}  # anchor ref_pos -> inserted query string after anchor
    deletions = {}   # deletion start ref_pos -> deleted reference sequence

    for n, op in _vcf_parse_cigar(row["cigar"]):
        if op in "M=X":
            for i in range(n):
                if qpos + i < len(seq):
                    bases[ref_pos + i] = seq[qpos + i]
            ref_pos += n
            qpos += n
        elif op == "I":
            ins = seq[qpos:qpos+n]
            anchor = ref_pos - 1
            if anchor >= int(row["start"]):
                insertions[anchor] = insertions.get(anchor, "") + ins
            qpos += n
        elif op in "DN":
            deleted = ref.fetch(contig, ref_pos, ref_pos + n - 1)
            deletions[ref_pos] = deleted
            ref_pos += n
        elif op == "S":
            qpos += n
        elif op in "HP":
            pass
    return {"bases": bases, "insertions": insertions, "deletions": deletions}


def _vcf_load_microhap_genotypes(microhap_genotypes_tsv: Path):
    with Path(microhap_genotypes_tsv).open() as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        samples = header[1:]
        by_locus = {}
        for row in rdr:
            if not row:
                continue
            loc = row[0]
            by_locus[loc] = {s: (row[i+1] if i + 1 < len(row) else "000000") for i, s in enumerate(samples)}
    return samples, by_locus


def _vcf_is_called_microhap_gt(gt: str) -> bool:
    gt = (gt or "").strip()
    return len(gt) == 6 and gt.isdigit() and gt != "000000"


def _vcf_split_microhap_gt(gt: str):
    return int(gt[:3]), int(gt[3:])


def _vcf_build_variant_records_for_locus(locus: str, allele_rows: list[dict], catalog_codes: set[int], ref: _VcfIndexedFasta):
    """Discover SNP/indel records for one locus and map each allele code to VCF allele indexes."""
    passing = [r for r in allele_rows if r.get("status") == "PASS"]
    if not passing:
        return [], {"locus": locus, "status": "no_passing_alignments", "n_records": 0}

    contigs = {r["contig"] for r in passing}
    strands = {r["strand"] for r in passing}
    if len(contigs) != 1:
        return [], {"locus": locus, "status": "alleles_map_to_multiple_contigs", "n_records": 0}
    if len(strands) != 1:
        return [], {"locus": locus, "status": "alleles_map_to_both_strands", "n_records": 0}

    contig = next(iter(contigs))
    projections = {int(r["allele_code"]): _vcf_project_alignment(r, ref) for r in passing}
    mapped_codes = set(projections)
    if not mapped_codes:
        return [], {"locus": locus, "status": "no_projected_alleles", "n_records": 0}

    snp_sites = set()
    insertion_sites = set()
    deletion_sites = set()
    for code, proj in projections.items():
        for pos, base in proj["bases"].items():
            rb = ref.base(contig, pos)
            if base in {"A", "C", "G", "T"} and rb in {"A", "C", "G", "T"} and base != rb:
                snp_sites.add(pos)
        insertion_sites.update(proj["insertions"].keys())
        deletion_sites.update(proj["deletions"].keys())

    records = []

    # SNP records
    for pos in sorted(snp_sites):
        ref_base = ref.base(contig, pos)
        alt_bases = []
        for code in sorted(mapped_codes):
            b = projections[code]["bases"].get(pos)
            if b in {"A", "C", "G", "T"} and b != ref_base and b not in alt_bases:
                alt_bases.append(b)
        if not alt_bases:
            continue
        allele_index_by_code = {}
        for code in sorted(mapped_codes):
            b = projections[code]["bases"].get(pos)
            if b == ref_base:
                allele_index_by_code[code] = 0
            elif b in alt_bases:
                allele_index_by_code[code] = alt_bases.index(b) + 1
            else:
                allele_index_by_code[code] = None
        records.append({
            "locus": locus,
            "contig": contig,
            "pos": pos,
            "id": f"{locus}:SNP:{pos}",
            "ref": ref_base,
            "alts": alt_bases,
            "type": "SNP",
            "allele_index_by_code": allele_index_by_code,
        })

    # Insertion records: VCF anchor is the preceding reference base.
    for anchor in sorted(insertion_sites):
        anchor_base = ref.base(contig, anchor)
        alt_alleles = []
        for code in sorted(mapped_codes):
            ins = projections[code]["insertions"].get(anchor, "")
            if ins:
                alt = anchor_base + ins
                if alt not in alt_alleles:
                    alt_alleles.append(alt)
        if not alt_alleles:
            continue
        allele_index_by_code = {}
        for code in sorted(mapped_codes):
            ins = projections[code]["insertions"].get(anchor, "")
            if not ins:
                allele_index_by_code[code] = 0
            else:
                alt = anchor_base + ins
                allele_index_by_code[code] = alt_alleles.index(alt) + 1 if alt in alt_alleles else None
        records.append({
            "locus": locus,
            "contig": contig,
            "pos": anchor,
            "id": f"{locus}:INS:{anchor}",
            "ref": anchor_base,
            "alts": alt_alleles,
            "type": "INS",
            "allele_index_by_code": allele_index_by_code,
        })

    # Deletion records: VCF anchor is the base immediately before deletion start.
    for del_start in sorted(deletion_sites):
        if del_start <= 1:
            continue
        anchor = del_start - 1
        alt_alleles = []
        ref_allele = None
        for code in sorted(mapped_codes):
            deleted = projections[code]["deletions"].get(del_start, "")
            if deleted:
                full_ref = ref.base(contig, anchor) + deleted
                ref_allele = full_ref if ref_allele is None else ref_allele
                alt = ref.base(contig, anchor)
                if alt not in alt_alleles:
                    alt_alleles.append(alt)
        if not ref_allele or not alt_alleles:
            continue
        allele_index_by_code = {}
        for code in sorted(mapped_codes):
            deleted = projections[code]["deletions"].get(del_start, "")
            if deleted:
                alt = ref.base(contig, anchor)
                allele_index_by_code[code] = alt_alleles.index(alt) + 1 if alt in alt_alleles else None
            else:
                allele_index_by_code[code] = 0
        records.append({
            "locus": locus,
            "contig": contig,
            "pos": anchor,
            "id": f"{locus}:DEL:{del_start}",
            "ref": ref_allele,
            "alts": alt_alleles,
            "type": "DEL",
            "allele_index_by_code": allele_index_by_code,
        })

    records.sort(key=lambda r: (r["contig"], int(r["pos"]), r["type"], r["id"]))
    return records, {"locus": locus, "status": "PASS", "n_records": len(records)}


def export_reference_anchored_vcf(
    catalog_csv: Path,
    microhap_genotypes_tsv: Path,
    reference_fasta: Path,
    out_vcf: Path,
    out_placements_tsv: Path,
    out_variant_map_tsv: Path,
    out_excluded_loci_tsv: Path,
    bwa_path: str = "bwa",
    build_bwa_index: bool = False,
    min_mapq: int = 20,
    min_pid: float = 0.95,
    min_aln_frac: float = 0.80,
    workers: int = 1,
    exclude_loci=None,
):
    """Export a reference-positioned VCF by projecting catalog microhap genotypes through BWA alignments.

    This output treats the existing microhap caller as the genotype engine, then decomposes each
    catalog haplotype into reference-positioned SNP/indel states. It is intended as an optional
    interoperability export, not as a replacement for the microhap genotype tables.
    """
    import tempfile

    exclude_loci = set(exclude_loci or [])
    catalog_by_locus = _vcf_load_catalog_by_locus(catalog_csv)
    if exclude_loci:
        catalog_by_locus = {loc: d for loc, d in catalog_by_locus.items() if loc not in exclude_loci}
    samples, gt_by_locus = _vcf_load_microhap_genotypes(microhap_genotypes_tsv)
    bwa_exe = _vcf_check_bwa_and_index(Path(reference_fasta), bwa_path=bwa_path, build_index=build_bwa_index)

    Path(out_vcf).parent.mkdir(parents=True, exist_ok=True)
    Path(out_placements_tsv).parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="gtseq_microhap_vcf_") as td:
        td = Path(td)
        cat_fa = td / "catalog_alleles.fa"
        sam = td / "catalog_alleles.sam"
        _vcf_write_catalog_fasta(catalog_by_locus, cat_fa)
        _vcf_run_bwa_mem(cat_fa, Path(reference_fasta), sam, bwa_exe=bwa_exe, threads=max(1, int(workers or 1)))
        placement_rows = _vcf_sam_to_alignments(sam, catalog_by_locus, min_mapq=min_mapq, min_aln_frac=min_aln_frac, min_pid=min_pid)

    _vcf_write_placements(placement_rows, out_placements_tsv)

    rows_by_locus = defaultdict(list)
    for r in placement_rows:
        rows_by_locus[r["locus"]].append(r)

    all_records = []
    excluded = []
    ref = _VcfIndexedFasta(Path(reference_fasta))
    try:
        for locus in sorted(catalog_by_locus):
            records, status = _vcf_build_variant_records_for_locus(
                locus=locus,
                allele_rows=rows_by_locus.get(locus, []),
                catalog_codes=set(catalog_by_locus[locus]),
                ref=ref,
            )
            if status.get("status") != "PASS" or not records:
                excluded.append({
                    "locus": locus,
                    "status": status.get("status", "no_variant_records"),
                    "n_records": int(status.get("n_records", 0)),
                })
            all_records.extend(records)
    finally:
        ref.close()

    # Write variant map before VCF for debugging/provenance.
    with Path(out_variant_map_tsv).open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus", "contig", "pos", "id", "type", "ref", "alts", "allele_code_to_vcf_index"])
        for rec in sorted(all_records, key=lambda r: (r["contig"], int(r["pos"]), r["id"])):
            code_map = ";".join(
                f"{code}:{'.' if idx is None else idx}" for code, idx in sorted(rec["allele_index_by_code"].items())
            )
            w.writerow([rec["locus"], rec["contig"], rec["pos"], rec["id"], rec["type"], rec["ref"], ",".join(rec["alts"]), code_map])

    with Path(out_excluded_loci_tsv).open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus", "status", "n_records"])
        for r in sorted(excluded, key=lambda x: x["locus"]):
            w.writerow([r["locus"], r["status"], r["n_records"]])

    # Write VCF.
    all_records.sort(key=lambda r: (r["contig"], int(r["pos"]), r["id"]))
    with Path(out_vcf).open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=gtseq_microhap_catalog_and_call.py_{__version__}\n")
        out.write(f"##reference={Path(reference_fasta)}\n")
        out.write("##INFO=<ID=LOCUS,Number=1,Type=String,Description=\"GT-seq microhap locus name\">\n")
        out.write("##INFO=<ID=VTYPE,Number=1,Type=String,Description=\"Variant type projected from catalog haplotypes\">\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype projected from called catalog microhap alleles\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        if samples:
            out.write("\t" + "\t".join(samples))
        out.write("\n")

        for rec in all_records:
            loc = rec["locus"]
            alt = ",".join(rec["alts"])
            info = f"LOCUS={loc};VTYPE={rec['type']}"
            row = [rec["contig"], str(rec["pos"]), rec["id"], rec["ref"], alt, ".", "PASS", info, "GT"]
            sample_gts = []
            for sample in samples:
                mh_gt = gt_by_locus.get(loc, {}).get(sample, "000000")
                if not _vcf_is_called_microhap_gt(mh_gt):
                    sample_gts.append("./.")
                    continue
                c1, c2 = _vcf_split_microhap_gt(mh_gt)
                i1 = rec["allele_index_by_code"].get(c1)
                i2 = rec["allele_index_by_code"].get(c2)
                if i1 is None or i2 is None:
                    sample_gts.append("./.")
                else:
                    sample_gts.append(f"{i1}/{i2}")
            out.write("\t".join(row + sample_gts) + "\n")

    return {
        "vcf": Path(out_vcf),
        "placements": Path(out_placements_tsv),
        "variant_map": Path(out_variant_map_tsv),
        "excluded_loci": Path(out_excluded_loci_tsv),
        "n_records": len(all_records),
        "n_excluded_loci": len(excluded),
    }

def main():

    import argparse

    class _DefaultsForOptionalOnly(argparse.ArgumentDefaultsHelpFormatter):
        def _get_help_string(self, action):
            help_text = action.help or ""
            if action.default is not argparse.SUPPRESS and not action.required:
                if "%(default)" not in help_text and "default:" not in help_text.lower():
                    help_text += " (default: %(default)s)"
            return help_text

    ap = argparse.ArgumentParser(
        description="GT-seq microhap catalog, genotype calling, sexID, parasite detection, single-SNP export, and optional reference-anchored VCF export",
        formatter_class=_DefaultsForOptionalOnly,
    )
    ap.add_argument("--version", action="version", version=f"gtseq_microhap_catalog_and_call.py {__version__}")

    ap.add_argument("--primers", required=True, help="CSV/TSV with locus,fwd_primer,rev_primer")
    ap.add_argument("--indir", required=True, help="Directory of FASTQ(.gz) pairs")
    ap.add_argument("--outdir", default="microhap_out", help="Output directory")
    ap.add_argument("--resolved-dir", default=None, help="Override resolved_fastqs directory; if omitted, uses <outdir>/resolved_fastqs")
    ap.add_argument("--min-amplicon", type=int, default=60, help="Minimum accepted primer-bounded amplicon length")
    ap.add_argument("--max-amplicon", type=int, default=250, help="Maximum accepted primer-bounded amplicon length")
    ap.add_argument("--min-overlap", type=int, default=12, help="Minimum R1/RC(R2) overlap length for paired-end stitching")
    ap.add_argument("--max-ov-mismatch-frac", type=float, default=0.05, help="Maximum mismatch fraction allowed in the R1/RC(R2) overlap")

    ap.add_argument("--workers", type=int, default=0, help="Parallel workers (0 = auto)")
    ap.add_argument("--force", action="store_true", help="Rebuild resolved FASTQs even if they exist")

    ap.add_argument("--min-catalog-depth", type=int, default=10, help="Minimum per-sample locus depth used when contributing alleles to the catalog")
    ap.add_argument("--catalog-max-len-diff", type=int, default=20,
                    help="Catalog safeguard: retain candidate alleles within this many bp of the dominant length class")
    ap.add_argument("--catalog-max-len-ratio", type=float, default=0.20,
                    help="Catalog safeguard: retain candidate alleles within this fractional length difference of the dominant length class")
    ap.add_argument("--catalog-length-bin-bp", type=int, default=2,
                    help="Catalog safeguard: cluster allele lengths within +/- this many bp when choosing the dominant length class")
    ap.add_argument("--catalog-min-pid-to-top", type=float, default=0.85,
                    help="Catalog safeguard: minimum global identity to the top supported allele for retaining close-length candidates")
    ap.add_argument("--disable-catalog-length-filter", action="store_true",
                    help="Disable final catalog pruning of length/divergence outliers")
    ap.add_argument("--hom-min", type=float, default=0.85, help="Catalog-building threshold: top allele fraction required to retain a single homozygous candidate allele")
    ap.add_argument("--het-min", type=float, default=0.30, help="Catalog-building threshold: minimum fraction required for both top alleles to retain a heterozygous candidate pair")
    ap.add_argument("--min-depth", type=int, default=10, help="Min A1+A2 depth for genotype calling")

    ap.add_argument(
        "--use-catalog",
        default=None,
        help="Use a pre-existing microhap catalog CSV for genotype calling instead of building a new one from resolved FASTQs"
    )

    ap.add_argument("--no-plots", action="store_true", help="Disable all plotting even if matplotlib is installed")
    ap.add_argument(
        "--a2-lo",
        type=float,
        default=10.0,
        help=(
            "Lower A2%% genotype threshold. Samples meeting --min-depth are called HOM "
            "when the second catalog allele is <= this percentage of A1+A2 reads "
            "(default: 10.0). Also drawn as a dashboard reference line."
        ),
    )
    ap.add_argument(
        "--a2-hi",
        type=float,
        default=25.0,
        help=(
            "Upper A2%% genotype threshold. Samples meeting --min-depth are called HET "
            "when the second catalog allele is >= this percentage of A1+A2 reads; "
            "values between --a2-lo and --a2-hi are called NC (default: 25.0). "
            "Also drawn as a dashboard reference line."
        ),
    )

    ap.add_argument("--sexid-markers", default=None, help="Optional text file with one sexID locus name per line")
    ap.add_argument("--sex-system", choices=["XY", "ZW"], default="XY", help="Interpret heterogametic calls as male/female (XY) or female/male (ZW)")
    ap.add_argument("--parasite-markers", default=None, help="Optional text file with one parasite locus name per line")

    ap.add_argument(
        "--single-snp-output",
        choices=["none", "first", "highest_maf", "both"],
        default="none",
        help=(
            "Optional single-SNP export mode derived from microhap catalog alleles: "
            "'first' = first valid biallelic SNP per locus, "
            "'highest_maf' = valid biallelic SNP with highest minor allele frequency in this dataset, "
            "'both' = write both outputs, "
            "'none' = disable single-SNP export"
        ),
    )

    ap.add_argument("--emit-vcf", action="store_true",
                    help="Optional reference-anchored VCF export. Requires --reference-fasta and BWA.")
    ap.add_argument("--reference-fasta", default=None,
                    help="Reference genome FASTA used with BWA for --emit-vcf placement and variant coordinates")
    ap.add_argument("--bwa-path", default="bwa",
                    help="BWA executable for --emit-vcf")
    ap.add_argument("--vcf-build-bwa-index", action="store_true",
                    help="If BWA index files are missing, run 'bwa index' on --reference-fasta")
    ap.add_argument("--vcf-min-mapq", type=int, default=20,
                    help="Minimum BWA MAPQ for catalog allele placements used in VCF export")
    ap.add_argument("--vcf-min-pid", type=float, default=0.95,
                    help="Minimum approximate alignment identity for catalog allele placements used in VCF export")
    ap.add_argument("--vcf-min-aln-frac", type=float, default=0.80,
                    help="Minimum fraction of catalog allele sequence aligned for VCF export")
    ap.add_argument("--vcf-out", default=None,
                    help="Output VCF path for --emit-vcf; if omitted, uses <outdir>/microhap_variants.vcf")
    ap.add_argument("--vcf-include-assay-markers", action="store_true",
                    help="Include sexID/parasite accessory marker loci in VCF export; by default they are excluded")
    ap.add_argument("--vcf-strict", action="store_true",
                    help="Raise an error if VCF export fails instead of warning and continuing")

    args = ap.parse_args()

    # Validate dependent command-line options before any analysis begins.
    # VCF export requires the original reference FASTA so catalog alleles can be
    # placed at genomic coordinates. Without this check, the pipeline can spend
    # hours resolving reads and calling genotypes before failing at the final VCF
    # export step. argparse.parser.error() prints usage, exits with status 2, and
    # prevents the run from starting.
    if args.emit_vcf and not args.reference_fasta:
        ap.error(
            "--emit-vcf requires --reference-fasta. "
            "Please provide the reference FASTA used for genomic placement of catalog alleles."
        )

    if not (0.0 <= args.a2_lo < args.a2_hi <= 100.0):
        ap.error("--a2-lo and --a2-hi must satisfy 0 <= --a2-lo < --a2-hi <= 100")

    # Hard-coded sexID defaults (kept intentionally off the CLI to simplify usage)
    args.sex_min_total_pb = 100
    args.sex_min_marker_reads = 2
    args.sex_min_marker_pot = 0.0005
    args.sex_provisional_min_markers = 1
    args.sex_strong_marker_pot = 0.005
    args.sex_homo_score_max = 0.10
    args.sex_hetero_score_min = 0.50
    args.sex_consensus_min_pid = 0.97

    # Hard-coded parasite defaults (parasite marker list is exposed on the CLI; thresholds remain simplified).
    args.parasite_min_total_pb = 100
    args.parasite_min_reads = 2
    args.parasite_low_pot = 0.0005
    args.parasite_med_pot = 0.0020
    args.parasite_high_pot = 0.0100
    args.parasite_very_high_pot = 0.0500

    print(f"gtseq_microhap_catalog_and_call.py version {__version__}", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    resolved_dir = Path(args.resolved_dir) if args.resolved_dir else (outdir / "resolved_fastqs")
    sample_stats_tsv = outdir / "sample_stats.tsv"
    built_catalog_csv = outdir / "microhap_catalog.csv"
    consensus_fa = outdir / "microhap_consensus.fasta"
    plots_dir = None if (args.no_plots or (not HAVE_MPL)) else (outdir / "plots")

    # New optional SNP-style outputs
    single_snp_first_tsv = outdir / "microhap_singleSNP_first.tsv"
    single_snp_highest_maf_tsv = outdir / "microhap_singleSNP_highestMAF.tsv"
    single_snp_metadata_tsv = outdir / "microhap_singleSNP_metadata.tsv"

    # New optional reference-anchored VCF outputs
    vcf_path = Path(args.vcf_out) if args.vcf_out else (outdir / "microhap_variants.vcf")
    vcf_placements_tsv = outdir / "microhap_vcf_locus_placements.tsv"
    vcf_variant_map_tsv = outdir / "microhap_vcf_variant_map.tsv"
    vcf_excluded_loci_tsv = outdir / "microhap_vcf_excluded_loci.tsv"

    if args.use_catalog:
        catalog_csv = Path(args.use_catalog).resolve()
        if not catalog_csv.exists():
            sys.exit(f"ERROR: --use-catalog file not found: {catalog_csv}")
    else:
        catalog_csv = built_catalog_csv

    primers = load_primers(Path(args.primers))
    pairs = find_fastq_pairs(Path(args.indir))
    if not pairs:
        sys.exit("ERROR: no FASTQ pairs found.")

    paired_samples = [(r1, r2, sample) for sample, (r1, r2) in pairs.items()]
    workers = args.workers if args.workers and args.workers > 0 else _n_workers(None)

    n_done, errors = resolve_fastqs_multiprocess(
        paired_samples=paired_samples,
        resolved_dir=resolved_dir,
        primers=primers,
        min_amplicon=args.min_amplicon,
        max_amplicon=args.max_amplicon,
        min_overlap=args.min_overlap,
        max_ov_mismatch_frac=args.max_ov_mismatch_frac,
        force=args.force,
        workers=workers,
        out_stats_path=sample_stats_tsv,
    )
    if errors:
        print("Some samples failed during resolve:", file=sys.stderr)
        for e in errors[:50]:
            print("  " + e, file=sys.stderr)
        if len(errors) > 50:
            print(f"  ... ({len(errors)-50} more)", file=sys.stderr)

    if args.use_catalog:
        print(f"Using existing catalog: {catalog_csv}")
    else:
        build_catalog_from_resolved(
            resolved_dir=resolved_dir,
            out_catalog_csv=catalog_csv,
            out_consensus_fa=consensus_fa,
            hom_min=args.hom_min,
            het_min=args.het_min,
            min_catalog_depth=args.min_catalog_depth,
            catalog_max_len_diff=args.catalog_max_len_diff,
            catalog_max_len_ratio=args.catalog_max_len_ratio,
            catalog_length_bin_bp=args.catalog_length_bin_bp,
            catalog_min_pid_to_top=args.catalog_min_pid_to_top,
            disable_catalog_length_filter=args.disable_catalog_length_filter,
        )

    sex_markers = load_locus_list(Path(args.sexid_markers)) if args.sexid_markers else []
    parasite_markers = load_locus_list(Path(args.parasite_markers)) if args.parasite_markers else []

    # Use the full attempted primer set as the locus-share histogram universe.
    # This adds zero-height bars for attempted loci that produced no primer-
    # bounded reads, making autosomal and accessory-marker dashboards directly
    # comparable and avoiding an overly optimistic distribution.
    locus_stats = compute_locus_readshare_stats_from_resolved(resolved_dir)
    locus_stats = _augment_locus_stats_with_zeros(locus_stats, list(primers.keys()))

    call_from_resolved_fastqs(
        resolved_dir=resolved_dir,
        catalog_csv=catalog_csv,
        out_prefix=str(outdir / "microhap"),
        min_depth=args.min_depth,
        locus_stats=locus_stats,
        plots_dir=plots_dir,
        a2_lo=args.a2_lo,
        a2_hi=args.a2_hi,
        locus_universe=list(primers.keys()),
        assay_exclude_loci=(list(sex_markers) + list(parasite_markers)),
    )

    # Optional SNP-style export derived from the microhap catalog and genotype calls.
    # Expected helper behavior:
    #   - mode="first" writes microhap_singleSNP_first.tsv
    #   - mode="highest_maf" writes microhap_singleSNP_highestMAF.tsv
    #   - mode="both" writes both files
    #   - always writes microhap_singleSNP_metadata.tsv when mode != "none"
    if args.single_snp_output != "none":
        try:
            export_single_snp_outputs(
                catalog_csv=catalog_csv,
                microhap_genotypes_tsv=outdir / "microhap_genotypes.tsv",
                microhap_phased_def_tsv=outdir / "microhap_phased_def.tsv",
                out_first_tsv=single_snp_first_tsv,
                out_highest_maf_tsv=single_snp_highest_maf_tsv,
                out_metadata_tsv=single_snp_metadata_tsv,
                mode=args.single_snp_output,
                exclude_loci=set(sex_markers) | set(parasite_markers),
            )
        except Exception as e:
            print(f"Warning: single-SNP export failed: {e}", file=sys.stderr)

    # Optional reference-anchored VCF export. This uses BWA to place catalog
    # alleles, then projects called microhap allele-code genotypes into single
    # SNP/indel records at genomic coordinates. Failure here does not invalidate
    # the normal microhap outputs unless --vcf-strict is used.
    vcf_results = None
    if args.emit_vcf:
        try:
            if not args.reference_fasta:
                raise RuntimeError("--emit-vcf requires --reference-fasta")
            vcf_exclude = set() if args.vcf_include_assay_markers else (set(sex_markers) | set(parasite_markers))
            vcf_results = export_reference_anchored_vcf(
                catalog_csv=catalog_csv,
                microhap_genotypes_tsv=outdir / "microhap_genotypes.tsv",
                reference_fasta=Path(args.reference_fasta),
                out_vcf=vcf_path,
                out_placements_tsv=vcf_placements_tsv,
                out_variant_map_tsv=vcf_variant_map_tsv,
                out_excluded_loci_tsv=vcf_excluded_loci_tsv,
                bwa_path=args.bwa_path,
                build_bwa_index=args.vcf_build_bwa_index,
                min_mapq=args.vcf_min_mapq,
                min_pid=args.vcf_min_pid,
                min_aln_frac=args.vcf_min_aln_frac,
                workers=args.workers if args.workers and args.workers > 0 else _n_workers(None),
                exclude_loci=vcf_exclude,
            )
            print(
                f"VCF export wrote {vcf_results['n_records']} record(s); "
                f"excluded {vcf_results['n_excluded_loci']} locus/loci.",
                file=sys.stderr,
            )
        except Exception as e:
            if args.vcf_strict:
                raise
            print(f"Warning: VCF export failed: {e}", file=sys.stderr)

    sexid_results = None
    parasite_results = None
    if sex_markers or parasite_markers:
        sample_locus_counts, sample_total_pb, sample_locus_seq_counts = scan_resolved_fastq_counts_and_sequences(resolved_dir)
        if sex_markers:
            sexid_results = summarize_sexid_markers(
                sample_locus_counts=sample_locus_counts,
                sample_total_pb=sample_total_pb,
                sex_markers=sex_markers,
                out_prefix=str(outdir / 'sexid'),
                sex_system=args.sex_system,
                min_total_pb_reads=args.sex_min_total_pb,
                min_marker_reads=args.sex_min_marker_reads,
                min_marker_pot=args.sex_min_marker_pot,
                provisional_min_markers=args.sex_provisional_min_markers,
                strong_marker_pot=args.sex_strong_marker_pot,
                homo_score_max=args.sex_homo_score_max,
                hetero_score_min=args.sex_hetero_score_min,
                sample_locus_seq_counts=sample_locus_seq_counts,
                sex_consensus_min_pid=args.sex_consensus_min_pid,
            )
            if sexid_results is not None:
                harmonize_sexid_outputs(outdir, sex_markers, sexid_results, sex_system=args.sex_system)
                append_assay_rows_to_wide_table(outdir / 'microhap_genotypes.tsv', '__sex_call__', sexid_results['sex_row'])
                append_assay_rows_to_wide_table(outdir / 'microhap_genotypes.tsv', '__gametic_class__', sexid_results['gametic_row'])
                append_assay_rows_to_wide_table(outdir / 'microhap_genotypes.tsv', '__sex_score__', sexid_results['score_row'])
                append_assay_rows_to_wide_table(outdir / 'microhap_phasedVARs.tsv', '__sex_call__', sexid_results['sex_row'])
                append_assay_rows_to_wide_table(outdir / 'microhap_phasedVARs.tsv', '__gametic_class__', sexid_results['gametic_row'])
                append_assay_rows_to_wide_table(outdir / 'microhap_phasedVARs.tsv', '__sex_score__', sexid_results['score_row'])
        if parasite_markers:
            parasite_results = summarize_parasite_markers(
                sample_locus_counts=sample_locus_counts,
                sample_total_pb_reads=sample_total_pb,
                parasite_markers=parasite_markers,
                outdir=str(outdir / 'parasite'),
                min_total_pb_reads=args.parasite_min_total_pb,
                min_reads_detect=args.parasite_min_reads,
                low_pot=args.parasite_low_pot,
                med_pot=args.parasite_med_pot,
                high_pot=args.parasite_high_pot,
                very_high_pot=args.parasite_very_high_pot,
            )

    if plots_dir is not None and HAVE_MPL and sexid_results is not None:
        try:
            padded_stats = _augment_locus_stats_with_zeros(locus_stats, list(primers.keys()))
            loci_sorted_by_share = padded_stats["sorted_loci_by_mean"]
            locus_pct_mean = padded_stats["mean_pct"]
            locus_pct_sd = padded_stats["sd_pct"]
            for marker in sex_markers:
                try:
                    plot_sexid_dashboard(
                        locus=marker,
                        score_rows=sexid_results['score_rows'],
                        call_rows=sexid_results['call_rows'],
                        locus_stats={
                            "loci_sorted_by_share": loci_sorted_by_share,
                            "mean_pct": locus_pct_mean,
                            "sd_pct": locus_pct_sd,
                        },
                        outdir=plots_dir,
                        sex_system=args.sex_system,
                    )
                except Exception as e:
                    try:
                        write_placeholder_png(Path(plots_dir)/f"{marker}_sexid_dashboard.png", marker, f"error: {e}")
                    except Exception:
                        pass
        except Exception as e:
            print(f"Warning: sexID plotting failed: {e}", file=sys.stderr)

    if plots_dir is not None and HAVE_MPL and parasite_results is not None:
        try:
            padded_stats = _augment_locus_stats_with_zeros(locus_stats, list(primers.keys()))
            for marker in parasite_markers:
                try:
                    plot_parasite_dashboard(
                        marker=marker,
                        marker_rows=parasite_results['marker_rows'],
                        locus_stats=padded_stats,
                        outdir=plots_dir,
                        low_pot=0.0005,
                        med_pot=0.0020,
                        high_pot=0.0100,
                        very_high_pot=0.0500,
                    )
                except Exception as e:
                    print(f"Warning: parasite plot failed for {marker}: {e}", file=sys.stderr)
                    try:
                        write_placeholder_png(Path(plots_dir)/f"{marker}_parasite_dashboard.png", marker, f"error: {e}")
                    except Exception as e2:
                        print(f"Warning: placeholder parasite plot failed for {marker}: {e2}", file=sys.stderr)
        except Exception as e:
            print(f"Warning: parasite plotting failed: {e}", file=sys.stderr)

    if plots_dir is not None and HAVE_MPL:
        try:
            plot_library_summary(
                sample_stats_tsv=sample_stats_tsv,
                a2_metrics_tsv=outdir / "microhap_a2_metrics.tsv",
                primers_csv=Path(args.primers),
                out_prefix=str(outdir / "microhap"),
            )
        except Exception as e:
            print(f"Warning: plot_library_summary failed: {e}", file=sys.stderr)

    print("Done.")
    print(f"- Resolved FASTQs: {resolved_dir}")
    print(f"- Catalog used:    {catalog_csv}")
    if not args.use_catalog:
        print(f"- Consensus FASTA: {consensus_fa}")
    print(f"- Metrics:         {outdir / 'microhap_a2_metrics.tsv'}")
    print(f"- Genotypes:       {outdir / 'microhap_genotypes.tsv'}")
    print(f"- Phased defs:     {outdir / 'microhap_phased_def.tsv'}")
    print(f"- Phased calls:    {outdir / 'microhap_phased.tsv'}")
    print(f"- Phased VARs:     {outdir / 'microhap_phasedVARs.tsv'}")

    if args.single_snp_output in ("first", "both"):
        print(f"- SNP calls 1st:   {single_snp_first_tsv}")
    if args.single_snp_output in ("highest_maf", "both"):
        print(f"- SNP calls MAF:   {single_snp_highest_maf_tsv}")
    if args.single_snp_output != "none":
        print(f"- SNP metadata:    {single_snp_metadata_tsv}")

    if args.emit_vcf:
        print(f"- VCF:             {vcf_path}")
        print(f"- VCF placements:  {vcf_placements_tsv}")
        print(f"- VCF variant map: {vcf_variant_map_tsv}")
        print(f"- VCF exclusions:  {vcf_excluded_loci_tsv}")

    if sexid_results is not None:
        print(f"- SexID calls:     {outdir / 'sexid_calls.tsv'}")
        print(f"- SexID p_OT:      {outdir / 'sexid_marker_pot.tsv'}")
        print(f"- SexID scores:    {outdir / 'sexid_sample_marker_scores.tsv'}")
    if parasite_results is not None:
        print(f"- Parasite calls:  {outdir / 'parasite_calls.tsv'}")
        print(f"- Parasite summary:{outdir / 'parasite_summary.tsv'}")
    if plots_dir is not None:
        print(f"- Plots:           {plots_dir}")

if __name__ == "__main__":
    main()
