#!/usr/bin/env python3
"""
gtseq_microhap_catalog_and_call.py

Alignment-free microhaplotype genotyping for GT-seq amplicon sequencing data.

This script implements an analysis pipeline for generating microhaplotype
genotypes directly from paired-end GT-seq FASTQ files without requiring
read alignment or conventional variant calling.

Pipeline overview
-----------------
1. Primer-bounded paired-end read resolution
   - Identify reads beginning with locus-specific forward primers
   - Confirm reverse primer in paired read
   - Merge paired reads to reconstruct full amplicon sequences

2. Allele discovery and catalog construction
   - Identify unique amplicon sequences within each sample
   - Rank sequences by read abundance
   - Retain up to two candidate alleles per locus under a diploid model
   - Aggregate across samples to construct a catalog of haplotypes

3. Genotype inference
   - Perform a second pass through resolved reads
   - Assign reads to catalog alleles using exact sequence matching
   - Infer diploid genotypes based on allele abundance ratios

4. Microhaplotype extraction
   - Align catalog haplotypes
   - Identify SNP and indel variation
   - Generate phased microhaplotype representations

Outputs (written to --outdir)
------------------------------
sample_stats.tsv
microhap_catalog.csv
microhap_consensus.fasta
microhap_a2_metrics.tsv
microhap_genotypes.tsv
microhap_phased.tsv
microhap_phasedSNPs.tsv

Optional outputs (if matplotlib available)
------------------------------------------
plots/  (per-locus dashboards)
microhap_on_target_fraction.png
microhap_raw_reads_hist.png
...

Repository
----------
Source code and documentation:
https://github.com/GTseq/gtseq_microhap

Citation
--------
If you use this software, please cite:

Campbell N. (2026) Alignment-free microhaplotype genotyping for
GT-seq amplicon sequencing data. GTseek LLC.

Version
-------
0.2.1

Author
------
Nathan Campbell
GTseek LLC

"""

__version__ = "0.2.1"

import os, sys, re, gzip, csv, collections, difflib
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
def trim_r1_short(seq: str, fwd: str, rev: str) -> Optional[Tuple[str, int]]:
    s = seq.upper()
    if not s.startswith(fwd):
        return None
    start = len(fwd)
    rev_rc = rc(rev)
    j = s.find(rev_rc, start)
    if j == -1:
        return None
    end = j + len(rev_rc)
    return (s[:end], end)

def r2_starts_with_rev(seq: str, rev: str) -> bool:
    return seq.upper().startswith(rev.upper())

def overlap_merge(a: str, b: str, min_ov: int, max_mismatch_frac: float) -> Optional[str]:
    """Ungapped suffix/prefix overlap merge; allows substitutions."""
    a = a.upper(); b = b.upper()
    max_len = min(len(a), len(b))

    # Stage 1: exact overlap (fast path)
    for ov in range(max_len, min_ov - 1, -1):
        if a[-ov:] == b[:ov]:
            return a + b[ov:]

    # Stage 2: allow mismatches
    for ov in range(max_len, min_ov - 1, -1):
        mism = 0
        limit = max(1, int(max_mismatch_frac * ov))
        aa = a[-ov:]; bb = b[:ov]
        for x, y in zip(aa, bb):
            if x != y:
                mism += 1
                if mism > limit:
                    break
        else:
            return a + b[ov:]

    return None

def _gapped_merge_from_alignment(a: str, b: str, aln_a: str, aln_b: str) -> str:
    merged_cols = []
    for ca, cb in zip(aln_a, aln_b):
        if ca == "-" and cb == "-":
            continue
        if ca == "-":
            merged_cols.append(cb)
        elif cb == "-":
            merged_cols.append(ca)
        else:
            merged_cols.append(ca)
    overlap_seq = "".join(merged_cols)

    a_ungapped = aln_a.replace("-", "")
    b_ungapped = aln_b.replace("-", "")

    idx_a = a.rfind(a_ungapped)
    prefix = a[:idx_a] if idx_a != -1 else a

    if b.startswith(b_ungapped):
        suffix = b[len(b_ungapped):]
    else:
        idx_b = b.find(b_ungapped)
        suffix = b[idx_b + len(b_ungapped):] if idx_b != -1 else ""

    return (prefix + overlap_seq + suffix).upper()

def overlap_merge_gapped(
    a: str,
    b: str,
    min_ov: int,
    max_mismatch_frac: float,
    max_indels: int,
    window: int = 140,
) -> Optional[str]:
    """Optional slower fallback; uses Biopython pairwise2 if available (lazy import)."""
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from Bio import pairwise2  # deprecated but still present in many installs
    except Exception:
        return None

    a = a.upper(); b = b.upper()
    if len(a) < min_ov or len(b) < min_ov:
        return None

    a_w = a[-min(window, len(a)):]
    b_w = b[:min(window, len(b))]

    alns = pairwise2.align.localms(a_w, b_w, 2, -3, -4, -1, one_alignment_only=True)
    if not alns:
        return None
    aln_a, aln_b, _score, _start, _end = alns[0]

    aligned_cols = 0
    mism = 0
    indel_cols = 0
    for ca, cb in zip(aln_a, aln_b):
        aligned_cols += 1
        if ca == "-" or cb == "-":
            indel_cols += 1
        elif ca != cb:
            mism += 1

    if aligned_cols < min_ov:
        return None
    if indel_cols > max_indels:
        return None
    if mism > max(1, int(max_mismatch_frac * aligned_cols)):
        return None

    merged = _gapped_merge_from_alignment(a_w, b_w, aln_a, aln_b)
    a_prefix = a[:-len(a_w)] if len(a) > len(a_w) else ""
    b_suffix = b[len(b_w):] if len(b) > len(b_w) else ""
    return (a_prefix + merged + b_suffix).upper()

def build_merged_amplicon(
    r1_seq: str,
    r2_seq: str,
    fwd: str,
    rev: str,
    min_ov: int,
    max_mismatch_frac: float,
    use_gapped_merge: bool = False,
    max_indels: int = 3,
) -> Optional[Tuple[str,int]]:
    """Return (amplicon_sequence, end1_for_quality)."""
    t = trim_r1_short(r1_seq, fwd, rev)
    if t:
        return t

    if not r2_starts_with_rev(r2_seq, rev):
        return None

    b = rc(r2_seq.upper())
    a = r1_seq.upper()

    merged = overlap_merge(a, b, min_ov, max_mismatch_frac)
    if merged is None and use_gapped_merge:
        merged = overlap_merge_gapped(a, b, min_ov, max_mismatch_frac, max_indels=max_indels)
    if merged is None:
        return None

    rev_rc = rc(rev)

    if not merged.startswith(fwd.upper()):
        return None

    if not merged.endswith(rev_rc):
        start = len(fwd)
        j = merged.find(rev_rc, start)
        if j == -1:
            return None
        merged = merged[: j + len(rev_rc)]

    end1 = min(len(r1_seq), len(merged))
    return (merged, end1)

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
     use_gapped_merge, max_indels, force) = args
    return process_sample_pair(
        r1, r2, sample,
        primers=primers,
        min_amplicon=min_amplicon,
        max_amplicon=max_amplicon,
        min_overlap=min_overlap,
        max_ov_mismatch_frac=max_ov_mismatch_frac,
        use_gapped_merge=use_gapped_merge,
        max_indels=max_indels,
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
    use_gapped_merge: bool = False,
    max_indels: int = 3,
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
         use_gapped_merge, max_indels, force)
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
    use_gapped_merge=False,
    max_indels=3,
    out_dir=Path("resolved_fastqs"),
    force=False,
    **_ignored,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"resolved_{sample_name}.fastq.gz"
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    if out_path.exists() and not force:
        return {"sample": sample_name, "read_pairs": 0, "resolved": 0, "skipped_existing": True}

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
            _h2, s2, _p2, _q2 = rec2

            for locus, (fwd, rev) in primers.items():
                if not s1.upper().startswith(fwd.upper()):
                    continue

                built = build_merged_amplicon(
                    s1, s2, fwd, rev,
                    min_ov=min_overlap,
                    max_mismatch_frac=max_ov_mismatch_frac,
                    use_gapped_merge=use_gapped_merge,
                    max_indels=max_indels,
                )
                if not built:
                    continue

                amplicon, _end1 = built
                L = len(amplicon)
                if L < min_amplicon or L > max_amplicon:
                    break

                resolved += 1
                out_h = h1.split()[0] + f"|locus={locus}\n"

                q = q1.rstrip()
                if len(q) < L:
                    q = q + ("I" * (L - len(q)))

                out_fh.write(out_h)
                out_fh.write(amplicon + "\n")
                out_fh.write("+\n")
                out_fh.write(q[:L] + "\n")
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
                _seq = fh.readline()
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

def build_catalog_from_resolved(
    resolved_dir,
    out_catalog_csv,
    out_consensus_fa,
    hom_min=0.85,
    het_min=0.30,
    min_catalog_depth=10,
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

    with out_catalog_csv.open("w") as out_csv, out_consensus_fa.open("w") as out_fa:
        out_csv.write("locus,allele_code,sequence,support_count,total_support\n")
        for locus in sorted(locus_seq_support.keys()):
            seq_counter = locus_seq_support[locus]
            if not seq_counter:
                continue
            sorted_seqs = sorted(seq_counter.items(), key=lambda kv: kv[1], reverse=True)
            total_support = locus_total_support[locus]
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
):
    """Second pass: call genotypes from resolved FASTQs against the catalog.

    Also emits optional phased hap strings derived from catalog allele sequences:
      - *_phased_def.tsv : variable-site definitions per locus
      - *_phased.tsv     : per-sample hap strings for each called locus
    """
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
                if a2_code is not None and frac_a2 >= 0.25:
                    call = "HET"
                    gt = f"{min(a1_code, a2_code)}{max(a1_code, a2_code)}"
                elif (a1_ct / denom) >= 0.90:
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

    # phased SNP hap matrix (locus x sample), analogous to microhap_genotypes.tsv
    # Each cell is the SNP-only hap string for the called genotype (e.g. "CC|CT").
    # Missing/uncalled genotypes are always reported as ".|.".
    phased_hap_snp: dict[str, dict[str, str]] = defaultdict(dict)
    for r in phased_rows:
        phased_hap_snp[r["locus"]][r["sample"]] = r.get("hap_snp", ".|.")

    with Path(f"{out_prefix}_phasedSNPs.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus"] + all_samples)
        for loc in loci_to_write:
            d = phased_hap_snp.get(loc, {})
            row = []
            for s in all_samples:
                v = d.get(s, ".|.")
                # Normalize any legacy/edge-case encodings to '.|.'
                if not v or v == "." or v == "|" or v == ".|" or v == "|.":
                    v = ".|."
                row.append(v)
            w.writerow([loc] + row)


    # ---- plots (optional) ----
    if plots_dir is not None and HAVE_MPL:
        plots_dir = Path(plots_dir)
        plots_dir.mkdir(parents=True, exist_ok=True)

        loci_sorted = locus_stats.get("sorted_loci_by_mean", [])
        locus_pct_mean = locus_stats.get("mean_pct", {})
        locus_pct_sd = locus_stats.get("sd_pct", {})

        for loc in PROG(loci_to_write, desc="Plotting per-locus dashboards", unit="locus", total=len(loci_to_write)):
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
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)

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
def _nw_align_global(ref: str, seq: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str]:
    """Needleman-Wunsch global alignment (deterministic tie-break).

    Returns (aligned_ref, aligned_seq) with '-' gaps.
    For our short GT-seq amplicons this is fast enough and FAR more reliable than difflib,
    especially in the presence of repeats/indels.
    """
    n = len(ref)
    m = len(seq)
    # DP score matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    # Trace: 0=diag, 1=up (gap in seq), 2=left (gap in ref)
    tr = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
        tr[i][0] = 1
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
        tr[0][j] = 2

    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            sj = seq[j - 1]
            s_diag = dp[i - 1][j - 1] + (match if ri == sj else mismatch)
            s_up = dp[i - 1][j] + gap
            s_left = dp[i][j - 1] + gap

            # Deterministic tie-break: diag > up > left
            best = s_diag
            bt = 0
            if s_up > best:
                best = s_up
                bt = 1
            elif s_up == best and bt != 0:
                pass
            if s_left > best:
                best = s_left
                bt = 2
            elif s_left == best and bt not in (0, 1):
                pass

            dp[i][j] = best
            tr[i][j] = bt

    # Traceback
    i, j = n, m
    aref = []
    aseq = []
    while i > 0 or j > 0:
        bt = tr[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and bt == 0:
            aref.append(ref[i - 1])
            aseq.append(seq[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or bt == 1):
            aref.append(ref[i - 1])
            aseq.append("-")
            i -= 1
        else:
            aref.append("-")
            aseq.append(seq[j - 1])
            j -= 1

    aref.reverse()
    aseq.reverse()
    return "".join(aref), "".join(aseq)


def _gap_profile_before_ref_bases(aligned_ref: str, ref_len: int) -> List[int]:
    """For an aligned_ref string, return gaps_before[i] = number of '-' before ref base i (0..ref_len),
    with i==ref_len representing trailing gaps after the last ref base.
    """
    gaps_before = [0] * (ref_len + 1)
    r_i = 0
    run = 0
    for ch in aligned_ref:
        if ch == "-":
            run += 1
        else:
            if r_i <= ref_len:
                gaps_before[r_i] = run
            run = 0
            r_i += 1
    gaps_before[ref_len] = run
    return gaps_before


def _expand_pair_to_gref(aligned_ref: str, aligned_seq: str, gref: str, ref: str) -> str:
    """Expand a (aligned_ref, aligned_seq) pair to match gref (a gapped version of ref).

    gref is constructed as ref with additional '-' runs inserted before each ref base (and possibly at end).
    This function inserts extra '-' in aligned_seq where gref contains more gaps than aligned_ref.
    """
    ref_len = len(ref)
    # Precompute max gaps-before profile from gref
    max_gaps = _gap_profile_before_ref_bases(gref, ref_len)

    # Parse the pairwise alignment into per-ref-base components
    # For each ref position i: collect insertion-run characters (where aligned_ref == '-') before the base,
    # then the base-aligned character (where aligned_ref is a base).
    ins_runs: List[str] = [""] * (ref_len + 1)  # ins before base i, and trailing at ref_len
    base_chars: List[str] = ["-"] * ref_len

    r_i = 0
    buf = []
    for rch, sch in zip(aligned_ref, aligned_seq):
        if rch == "-":
            # insertion relative to ref
            buf.append(sch)
        else:
            if r_i < ref_len:
                ins_runs[r_i] = "".join(buf)
                buf = []
                base_chars[r_i] = sch
            r_i += 1
    # trailing insertions
    ins_runs[ref_len] = "".join(buf)

    # Now expand to match gref
    out = []
    for i in range(ref_len):
        run = ins_runs[i]
        if len(run) < max_gaps[i]:
            run = run + ("-" * (max_gaps[i] - len(run)))
        out.append(run)
        out.append(base_chars[i])
    # trailing
    run = ins_runs[ref_len]
    if len(run) < max_gaps[ref_len]:
        run = run + ("-" * (max_gaps[ref_len] - len(run)))
    out.append(run)
    return "".join(out)


def _multi_align_to_ref(ref: str, seqs: List[str]) -> Tuple[str, List[str]]:
    """Progressive-to-ref multiple alignment.

    We align each seq to ref (pairwise NW), then build a gapped ref (gref) that is the UNION
    of all gap placements, and finally expand each aligned seq to that gref.

    Returns (gref, gapped_seqs) in the same order as seqs.
    """
    ref_len = len(ref)
    pairwise = [_nw_align_global(ref, s) for s in seqs]

    # Build max gap profile across all pairwise aligned refs
    max_gaps = [0] * (ref_len + 1)
    for aref, _ in pairwise:
        gb = _gap_profile_before_ref_bases(aref, ref_len)
        for i in range(ref_len + 1):
            if gb[i] > max_gaps[i]:
                max_gaps[i] = gb[i]

    # Construct gref
    parts = []
    for i in range(ref_len):
        parts.append("-" * max_gaps[i])
        parts.append(ref[i])
    parts.append("-" * max_gaps[ref_len])
    gref = "".join(parts)

    gapped_seqs = [_expand_pair_to_gref(aref, aseq, gref, ref) for aref, aseq in pairwise]
    return gref, gapped_seqs

def _nw_global_align(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> tuple[str, str]:
    """Simple Needleman–Wunsch global alignment (no affine gaps).

    Good enough here because GTseq alleles are short (typically < 200 bp) and
    differences are usually sparse (few SNPs / small indels).
    """
    a = a.upper()
    b = b.upper()
    n, m = len(a), len(b)
    # score matrix
    S = [[0] * (m + 1) for _ in range(n + 1)]
    # traceback: 0=diag, 1=up, 2=left
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = i * gap
        T[i][0] = 1
    for j in range(1, m + 1):
        S[0][j] = j * gap
        T[0][j] = 2

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            sd = S[i - 1][j - 1] + (match if ai == bj else mismatch)
            su = S[i - 1][j] + gap
            sl = S[i][j - 1] + gap
            if sd >= su and sd >= sl:
                S[i][j] = sd
                T[i][j] = 0
            elif su >= sl:
                S[i][j] = su
                T[i][j] = 1
            else:
                S[i][j] = sl
                T[i][j] = 2

    # traceback
    i, j = n, m
    aa = []
    bb = []
    while i > 0 or j > 0:
        tb = T[i][j]
        if i > 0 and j > 0 and tb == 0:
            aa.append(a[i - 1])
            bb.append(b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or tb == 1):
            aa.append(a[i - 1])
            bb.append("-")
            i -= 1
        else:
            aa.append("-")
            bb.append(b[j - 1])
            j -= 1
    return "".join(reversed(aa)), "".join(reversed(bb))



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
    current_boundary = 0
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
    boundary = 0
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
        last_ref_pos = 0  # for insertions (anchor after previous ref base)
        var_index = 0

        # First determine which columns are variable
        for col in range(aln_len):
            rch = master_ref[col]
            # update ref coordinate
            if rch != '-':
                ref_pos += 1
                last_ref_pos = ref_pos

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
                rch = master_ref[col]
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
def _parent_from_catalog_locus(loc: str) -> str:
    if loc.endswith("a") or loc.endswith("b"):
        return loc[:-1]
    return loc

def call_from_resolved_fastqs_v2(
    resolved_dir: Path,
    catalog_csv: Path,
    out_prefix: str,
    min_depth: int,
    locus_stats: dict,
    plots_dir: Path | None = None,
    a2_lo: float = 10.0,
    a2_hi: float = 25.0,
    locus_order: list[str] | None = None,
):
    catalog_child: dict[str, dict[str,int]] = defaultdict(dict)
    seq_to_child: dict[str, dict[str, tuple[str,int]]] = defaultdict(dict)

    with Path(catalog_csv).open() as fh:
        rdr = csv.DictReader(fh)
        cols = {k.lower(): k for k in rdr.fieldnames}
        for row in rdr:
            child = row[cols["locus"]]
            seq   = row[cols["sequence"]]
            code  = int(row[cols["allele_code"]])
            parent = _parent_from_catalog_locus(child)
            catalog_child[child][seq] = code
            seq_to_child[parent][seq] = (child, code)

    resolved_paths = sorted(Path(resolved_dir).glob("resolved_*.fastq.gz"))

    per_locus_sample_gt: dict[str, dict[str, str]] = defaultdict(dict)
    a2_rows = []
    locus_points_a1a2 = defaultdict(list)
    locus_points_depth_pct = defaultdict(list)

    locus_called_samples = defaultdict(int)
    locus_total_samples  = defaultdict(int)
    locus_background_sum = defaultdict(float)

    loci_to_iter = locus_order[:] if locus_order else sorted(catalog_child.keys())

    for fpath in PROG(resolved_paths, desc="Calling from resolved", unit="sample", total=len(resolved_paths)):
        sample = fpath.stem.replace("resolved_", "")

        sample_locus_cnts: dict[str, Counter] = defaultdict(Counter)
        sample_parent_offcat: dict[str, int] = defaultdict(int)
        sample_parent_total: dict[str, int] = defaultdict(int)

        with gzip.open(fpath, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip()
                _ = fh.readline(); _ = fh.readline()

                if "|locus=" not in h:
                    continue
                parent = h.strip().split("|locus=")[-1]
                if not parent:
                    continue
                sample_parent_total[parent] += 1

                if parent in seq_to_child and seq in seq_to_child[parent]:
                    child, _code = seq_to_child[parent][seq]
                    sample_locus_cnts[child][seq] += 1
                else:
                    sample_parent_offcat[parent] += 1

        for loc in loci_to_iter:
            cnts = sample_locus_cnts.get(loc, Counter())
            total_cat = sum(cnts.values())

            parent = _parent_from_catalog_locus(loc)
            off_cat = sample_parent_offcat.get(parent, 0)
            total_parent_reads = sample_parent_total.get(parent, 0)

            if total_parent_reads > 0:
                locus_total_samples[loc] += 1

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

                if total_parent_reads > 0:
                    locus_background_sum[loc] += (off_cat / total_parent_reads) if total_parent_reads else 0.0
                continue

            top = cnts.most_common(2)
            (a1_seq, a1_ct) = top[0]
            a1_code = catalog_child[loc][a1_seq]
            if len(top) > 1:
                (a2_seq, a2_ct) = top[1]; a2_code = catalog_child[loc][a2_seq]
            else:
                a2_ct, a2_code = 0, None

            denom = a1_ct + a2_ct
            frac_a2 = (a2_ct / denom) if denom > 0 else 0.0
            a2_pct  = frac_a2 * 100.0

            if denom < min_depth:
                call = "LOW"
                gt   = "000000"
            else:
                if a2_code is not None and frac_a2 >= 0.25:
                    call = "HET"
                    gt   = f"{min(a1_code, a2_code)}{max(a1_code, a2_code)}"
                elif (a1_ct / denom) >= 0.90:
                    call = "HOM"
                    gt   = f"{a1_code}{a1_code}"
                else:
                    call = "NC"
                    gt   = "000000"

            if call in ("HOM","HET"):
                locus_called_samples[loc] += 1

            per_locus_sample_gt[loc][sample] = gt
            a2_rows.append({
                "sample": sample,
                "locus": loc,
                "depth": denom,
                "A1_count": a1_ct,
                "A2_count": a2_ct,
                "A1_code": int(a1_code) if a1_code is not None else 0,
                "A2_code": int(a2_code) if a2_code is not None else 0,
                "A2_pct": round(a2_pct, 2),
                "call": call
            })

            locus_points_a1a2[loc].append((a1_ct, a2_ct, call, a1_code if call=="HOM" else None))
            locus_points_depth_pct[loc].append((denom, a2_pct, call, a1_code if call=="HOM" else None))

            other_catalog = total_cat - denom
            if total_parent_reads > 0:
                locus_background_sum[loc] += (off_cat + other_catalog) / total_parent_reads

    with Path(f"{out_prefix}_a2_metrics.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample","locus","depth","A1_count","A2_count","A1_code","A2_code","A2_pct","call"])
        for r in a2_rows:
            w.writerow([
                r["sample"], r["locus"], r["depth"],
                r["A1_count"], r["A2_count"],
                r["A1_code"], r["A2_code"],
                r["A2_pct"], r["call"]
            ])

    all_samples = sorted({s for d in per_locus_sample_gt.values() for s in d})
    loci_to_write = loci_to_iter

    with Path(f"{out_prefix}_genotypes.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus"] + all_samples)
        for loc in loci_to_write:
            d = per_locus_sample_gt.get(loc, {})
            w.writerow([loc] + [d.get(s, "000000") for s in all_samples])

    locus_call_rate = {}
    locus_background_pct = {}
    for loc in loci_to_write:
        tot = locus_total_samples.get(loc, 0)
        called = locus_called_samples.get(loc, 0)
        locus_call_rate[loc] = (called / tot * 100.0) if tot else 0.0
        bg_sum = locus_background_sum.get(loc, 0.0)
        locus_background_pct[loc] = (bg_sum / tot * 100.0) if tot else 0.0

    if plots_dir and HAVE_MPL:
        plots_dir.mkdir(parents=True, exist_ok=True)
        padded_stats = _augment_locus_stats_with_zeros(locus_stats, loci_to_write)
        loci_sorted_by_share = padded_stats["sorted_loci_by_mean"]
        locus_pct_mean = padded_stats["mean_pct"]
        locus_pct_sd   = padded_stats["sd_pct"]

        for loc in PROG(loci_to_write, desc="Plotting loci", unit="locus", total=len(loci_to_write)):
            allele_freq_codes = _allele_freq_from_genos(per_locus_sample_gt.get(loc, {}))
            try:
                plot_locus_dashboard(
                    locus=loc,
                    points_a1a2=locus_points_a1a2.get(loc, []),
                    depths_a2pct=locus_points_depth_pct.get(loc, []),
                    loci_sorted_by_share=loci_sorted_by_share,
                    locus_pct_mean=locus_pct_mean,
                    locus_pct_sd=locus_pct_sd,
                    allele_freq_codes=allele_freq_codes,
                    outdir=plots_dir,
                    extra_info={
                        "call_rate": locus_call_rate.get(loc, 0.0),
                        "background_pct": locus_background_pct.get(loc, 0.0),
                    },
                    a2_lo=a2_lo,
                    a2_hi=a2_hi,
                )
            except Exception as e:
                try:
                    write_placeholder_png(Path(plots_dir)/f"{loc}_dashboard.png", loc, f"error: {e}")
                except Exception:
                    pass

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
# CLI + main
# ---------------------------------------------------------------------
def _expand_locus_order(primers: dict, split_map: dict) -> list[str]:
    order = []
    for loc in primers.keys():
        if loc in split_map:
            order.append(f"{loc}a")
            order.append(f"{loc}b")
        else:
            order.append(loc)
    return order



def main():

    import argparse

    ap = argparse.ArgumentParser(description="GT-seq microhap catalog + call")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    ap.add_argument("--primers", required=True, help="CSV/TSV with locus,fwd_primer,rev_primer")
    ap.add_argument("--indir", required=True, help="Directory of FASTQ(.gz) pairs")
    ap.add_argument("--outdir", default="microhap_out", help="Output directory")
    ap.add_argument("--resolved-dir", default=None, help="Override resolved_fastqs directory (default: <outdir>/resolved_fastqs)")
    ap.add_argument("--min-amplicon", type=int, default=60)
    ap.add_argument("--max-amplicon", type=int, default=250)
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--max-ov-mismatch-frac", type=float, default=0.05)

    ap.add_argument("--gapped-merge", action="store_true", help="Enable slow gapped overlap merge fallback (Biopython)")
    ap.add_argument("--max-indels", type=int, default=3, help="Max indel columns allowed if --gapped-merge is enabled")
    ap.add_argument("--workers", type=int, default=0, help="Parallel workers (0 = auto)")
    ap.add_argument("--force", action="store_true", help="Rebuild resolved FASTQs even if they exist")

    ap.add_argument("--min-catalog-depth", type=int, default=10)
    ap.add_argument("--hom-min", type=float, default=0.85)
    ap.add_argument("--het-min", type=float, default=0.30)
    ap.add_argument("--min-depth", type=int, default=10, help="Min A1+A2 depth for genotype calling")

    ap.add_argument(
        "--use-catalog",
        default=None,
        help="Use a pre-existing microhap catalog CSV for genotype calling instead of building a new one from resolved FASTQs"
    )

    ap.add_argument("--no-plots", action="store_true", help="Disable all plotting even if matplotlib is installed")
    ap.add_argument("--a2-lo", type=float, default=10.0, help="Lower A2%% band for dashboards")
    ap.add_argument("--a2-hi", type=float, default=25.0, help="Upper A2%% band for dashboards")

    args = ap.parse_args()

    print(f"gtseq_microhap_catalog_and_call.py version {__version__}", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    resolved_dir = Path(args.resolved_dir) if args.resolved_dir else (outdir / "resolved_fastqs")
    sample_stats_tsv = outdir / "sample_stats.tsv"
    built_catalog_csv = outdir / "microhap_catalog.csv"
    consensus_fa = outdir / "microhap_consensus.fasta"
    plots_dir = None if (args.no_plots or (not HAVE_MPL)) else (outdir / "plots")

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
        use_gapped_merge=args.gapped_merge,
        max_indels=args.max_indels,
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
        )

    locus_stats = compute_locus_readshare_stats_from_resolved(resolved_dir)

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
    )

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
    print(f"- Phased SNPs:     {outdir / 'microhap_phasedSNPs.tsv'}")
    if plots_dir is not None:
        print(f"- Plots:           {plots_dir}")

if __name__ == "__main__":
    main()
