#!/usr/bin/env python3
"""
GT-seq microhap catalog + call

Pipeline:
1. load primers (locus,fwd_primer,rev_primer)
2. find FASTQ pairs in input dir
3. for each sample:
   - read R1/R2 in lockstep
   - for each read pair, try to make amplicon for ONE locus:
       * short-amplicon trim: R1 starts FWD, find RC(REV) in R1, trim
       * else long-amplicon merge: R2 must start REV, merge R1 with RC(R2)
         with min-overlap and max mismatch frac
       * length-gate amplicon [min_amplicon,max_amplicon]
   - write amplicon to resolved_fastqs/resolved_<sample>.fastq.gz
     as a 4-line FASTQ with header “@...|locus=<locus>”
   - tally per-locus counts for this sample
   - collect per-locus sequences for catalog (raw amplicon sequences)
4. after all samples:
   - build locus read-share stats (mean % + SD)
   - build catalog: for each locus, pull all sample-level counts for that locus,
     and **keep only alleles that were top-2 in at least one sample**
   - write microhap_catalog.csv
   - write microhap_consensus.fasta (consensus = most common allele per locus)
5. SECOND PASS:
   - for each resolved_<sample>.fastq.gz
     * count ONLY sequences that are present in the catalog for that locus
     * call genotype using %A2 thresholds
   - write microhap_a2_metrics.tsv
   - write microhap_genotypes.tsv
   - make per-locus dashboards

This is the “cleaned” version that **scores only catalog sequences** to avoid background.
"""

import os, sys, re, gzip, collections, csv
from math import sqrt
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter
from collections import defaultdict
import signal, gc, traceback, datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

# ---- Optional progress bars (tqdm if available; graceful fallback) ----
try:
    from tqdm import tqdm as _tqdm
    def PROG(iterable, desc="", unit="", total=None):
        return _tqdm(iterable, desc=desc, unit=unit, total=total)
except Exception:
    def PROG(iterable, desc="", unit="", total=None):
        # Minimal fallback: prints every 25 items (and on first/last) with a running counter
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
        # accept either (locus,fwd_primer,rev_primer) OR (locus,forward_primer,reverse_primer)
        # but we’ll normalize to these names:
        name_map = {
            "locus": None,
            "fwd_primer": None,
            "rev_primer": None,
        }
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
    if not s.startswith(fwd): return None
    start = len(fwd)
    rev_rc = rc(rev)
    j = s.find(rev_rc, start)
    if j == -1: return None
    end = j + len(rev_rc)
    return (s[:end], end)

def r2_starts_with_rev(seq: str, rev: str) -> bool:
    return seq.upper().startswith(rev.upper())

def overlap_merge(a: str, b: str, min_ov: int, max_mismatch_frac: float) -> Optional[str]:
    a = a.upper(); b = b.upper()
    max_len = min(len(a), len(b))
    best = None
    for ov in range(max_len, min_ov - 1, -1):
        mism = 0
        limit = max(1, int(max_mismatch_frac * ov))
        for i in range(ov):
            if a[len(a)-ov+i] != b[i]:
                mism += 1
                if mism > limit:
                    break
        else:
            best = ov
            break
    if best is None:
        return None
    return a + b[best:]

def build_merged_amplicon(r1_seq: str, r2_seq: str, fwd: str, rev: str,
                          min_ov: int, max_mismatch_frac: float) -> Optional[Tuple[str,int]]:
    # Try short-amplicon trim first
    t = trim_r1_short(r1_seq, fwd, rev)
    if t:
        return t
    # Else overlap-merge R1 with rc(R2)
    if not r2_starts_with_rev(r2_seq, rev):
        return None
    merged = overlap_merge(r1_seq.upper(), rc(r2_seq.upper()), min_ov, max_mismatch_frac)
    if merged is None: return None
    rev_rc = rc(rev)
    if not merged.startswith(fwd.upper()):
        return None
    if not merged.endswith(rev_rc):
        start = len(fwd)
        j = merged.find(rev_rc, start)
        if j == -1: return None
        merged = merged[: j + len(rev_rc)]
    end1 = min(len(r1_seq), len(merged))  # R1 qualities up to here
    return (merged, end1)
    
def _select_top2_for_catalog(cnts: "collections.Counter[str]", hom_min=0.85, het_min=0.30):
    """
    From a Counter of primer-bounded sequences for one sample×locus, return which
    sequences to count toward the *catalog* for this sample.

    Rules:
      - HOM evidence: top1 >= hom_min and (top2 absent or top2 < het_min) -> [top1]
      - HET evidence: top1 >= het_min and top2 >= het_min -> [top1, top2]
      - otherwise -> []
    """
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

    # Strong homozygote evidence
    if f1 >= hom_min and f2 < het_min:
        return [s1]

    # Heterozygote evidence
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

# --- DROP-IN: parallel resolver helpers ---

def _resolve_one_sample(args):
    """
    Worker wrapper. Returns (sample, ok, err_msg_or_None).
    """
    (sample, r1, r2, primers, resolved_dir,
     min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac) = args
    try:
        # Your process_sample_pair() already writes the resolved_<sample>.fastq.gz
        # and returns stats; we don't care about the stats here.
        process_sample_pair(
            r1, r2,
            primers=primers,
            min_amplicon=min_amplicon,
            max_amplicon=max_amplicon,
            min_overlap=min_overlap,
            max_ov_mismatch_frac=max_ov_mismatch_frac,
            resolved_dir=resolved_dir,
        )
        return (sample, True, None)
    except Exception as e:
        return (sample, False, str(e))

# --- TOP-LEVEL worker for ProcessPoolExecutor (must not be nested) ---
def _worker_resolve(args):
    """
    args = (r1, r2, sample, primers, resolved_dir,
            min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac, force)
    Returns whatever process_sample_pair returns (or raises).
    """
    (r1, r2, sample, primers, resolved_dir,
     min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac, force) = args
    # IMPORTANT: match process_sample_pair's signature exactly
    return process_sample_pair(
        r1, r2, sample,
        primers=primers,
        min_amplicon=min_amplicon,
        max_amplicon=max_amplicon,
        min_overlap=min_overlap,
        max_ov_mismatch_frac=max_ov_mismatch_frac,
        out_dir=resolved_dir,   # your function uses out_dir (not resolved_dir)
        force=force,
    )

def resolve_fastqs_multiprocess(
    paired_samples,               # list of (r1, r2, sample)
    resolved_dir: Path,
    primers: dict,
    min_amplicon: int = 60,
    max_amplicon: int = 250,
    min_overlap: int = 12,
    max_ov_mismatch_frac: float = 0.05,
    force: bool = False,
    workers: int | None = None,
):
    """
    Parallel first-pass resolver. Returns (n_done, errors:list[str])
    paired_samples: [(r1_path, r2_path, sample_name), ...]
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    resolved_dir.mkdir(parents=True, exist_ok=True)

    if workers is None:
        try:
            workers = max(1, (os.cpu_count() or 2) // 2)
        except Exception:
            workers = 1

    if not paired_samples:
        return 0, []

    # Build fully-picklable arg tuples for each job
    jobs = [
        (r1, r2, sample, primers, resolved_dir,
         min_amplicon, max_amplicon, min_overlap, max_ov_mismatch_frac, force)
        for (r1, r2, sample) in paired_samples
    ]

    done = 0
    errors: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_worker_resolve, args): args for args in jobs}
        try:
            # progress bar if tqdm is installed; else plain iterator
            try:
                from tqdm import tqdm as _tqdm
                iterator = _tqdm(as_completed(futs), total=len(futs),
                                 desc="Resolving FASTQs (parallel)", unit="sample")
            except Exception:
                iterator = as_completed(futs)

            for fut in iterator:
                args = futs[fut]
                try:
                    _ = fut.result()   # triggers exception here if worker failed
                    done += 1
                except Exception as e:
                    # args contains (r1, r2, sample, ...)
                    sample = args[2]
                    errors.append(f"{sample}: {e}")
        finally:
            pass

    return done, errors

# ---------------------------------------------------------------------
# 1st pass: process one sample, write resolved fastq, collect candidates
# ---------------------------------------------------------------------
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
    **_ignored,  # swallow any stray kwargs
):
    """
    Read paired FASTQs, extract primer-bounded amplicons for exactly one locus per read pair,
    and write them to resolved_fastqs/resolved_<sample>.fastq.gz as 4-line FASTQ entries.
    Returns simple stats dict.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"resolved_{sample_name}.fastq.gz"
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    # If not forcing and output exists, skip rework.
    if out_path.exists() and not force:
        return {"read_pairs": 0, "resolved": 0, "skipped_existing": True}

    # Open inputs & temp output
    fh1 = open_read(r1_path)
    fh2 = open_read(r2_path)
    import gzip
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

            h1, s1, p1, q1 = rec1
            _h2, s2, _p2, _q2 = rec2

            matched = False
            for locus, (fwd, rev) in primers.items():
                # short-amplicon path requires R1 to start with fwd
                if not s1.upper().startswith(fwd.upper()):
                    continue

                built = build_merged_amplicon(
                    s1, s2, fwd, rev,
                    min_ov=min_overlap,
                    max_mismatch_frac=max_ov_mismatch_frac,
                )
                if not built:
                    continue

                amplicon, _end1 = built
                L = len(amplicon)
                if L < min_amplicon or L > max_amplicon:
                    matched = True
                    break  # this pair matched primers but length-gated out; stop checking other loci

                # write one FASTQ record per resolved amplicon
                matched = True
                resolved += 1
                out_h = h1.split()[0] + f"|locus={locus}\n"
                # Ensure quality string at least as long as sequence (pad with 'I' if needed)
                q = q1.rstrip()
                if len(q) < L:
                    q = q + ("I" * (L - len(q)))

                out_fh.write(out_h)
                out_fh.write(amplicon + "\n")
                out_fh.write("+\n")
                out_fh.write(q[:L] + "\n")
                break  # stop after first locus match

            # if not matched: just move on to next pair

    finally:
        fh1.close()
        fh2.close()
        out_fh.close()

    # Atomically publish the file
    tmp_path.replace(out_path)

    return {"read_pairs": read_pairs, "resolved": resolved}

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
    return {
        "sorted_loci_by_mean": sorted_loci_by_mean,
        "mean_pct": mean_pct,
        "sd_pct": sd_pct,
    }

# ---------------------------------------------------------------------
# write catalog & consensus  (pruned to top-2-seen alleles)
# ---------------------------------------------------------------------
def write_catalog(catalog_candidates: dict, catalog_path: Path, consensus_path: Path):
    # Build rows
    rows = []
    for locus, seqdict in catalog_candidates.items():
        # seqdict: seq -> {samples:set(), reads:int}
        # we want to order alleles by total_reads desc
        items = []
        for seq, info in seqdict.items():
            items.append((seq, len(info["samples"]), info["reads"]))
        items.sort(key=lambda x: x[2], reverse=True)
        # assign 3-digit codes starting at 101
        for i, (seq, nsamp, nreads) in enumerate(items, start=1):
            code = 100 + i
            rows.append({
                "locus": locus,
                "allele_code": str(code),
                "sequence": seq,
                "total_samples": nsamp,
                "total_reads": nreads,
            })
    # write csv
    with catalog_path.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["locus", "allele_code", "sequence", "total_samples", "total_reads"])
        for r in sorted(rows, key=lambda x: (x["locus"], int(x["allele_code"]))):
            w.writerow([r["locus"], r["allele_code"], r["sequence"], r["total_samples"], r["total_reads"]])

    # write consensus fasta = allele 101 if present, else first
    with consensus_path.open("w") as out:
        loci = sorted({r["locus"] for r in rows})
        for locus in loci:
            # pick 101 if exists
            seq_101 = None
            for r in rows:
                if r["locus"] == locus and r["allele_code"] == "101":
                    seq_101 = r["sequence"]
                    break
            if seq_101 is None:
                # pick first row for this locus
                for r in rows:
                    if r["locus"] == locus:
                        seq_101 = r["sequence"]
                        break
            if seq_101 is None:
                seq_101 = "N"
            out.write(f">{locus}\n{seq_101}\n")
            
def build_catalog_from_resolved(
    resolved_dir: Path,
    out_catalog_csv: Path,
    out_consensus_fa: Path,
    hom_min: float = 0.85,   # top1 ≥ 85% -> HOM evidence
    het_min: float = 0.30,   # top1 & top2 ≥ 30% each -> HET evidence
):
    """
    Build the microhap catalog from resolved fastqs *only* using the strongest per-sample evidence:
      - HOM sample evidence: top1 >= hom_min and top2 < het_min -> count top1 only
      - HET sample evidence: top1 >= het_min and top2 >= het_min -> count top1 and top2
      - otherwise ignore this sample for catalog building at that locus
    """
    resolved_paths = sorted(Path(resolved_dir).glob("resolved_*.fastq.gz"))

    from collections import defaultdict, Counter
    catalog_counts: dict[str, Counter] = defaultdict(Counter)

    for fpath in resolved_paths:
        with gzip.open(fpath, "rt") as fh:
            sample_locus_counts: dict[str, Counter] = defaultdict(Counter)
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip()
                _ = fh.readline(); _ = fh.readline()
                loc = h.strip().split("|locus=")[-1] if "|locus=" in h else None
                if not loc:
                    continue
                sample_locus_counts[loc][seq] += 1

        for loc, cnts in sample_locus_counts.items():
            chosen = _select_top2_for_catalog(cnts, hom_min=hom_min, het_min=het_min)
            for s in chosen:
                catalog_counts[loc][s] += cnts[s]

    out_catalog_csv = Path(out_catalog_csv)
    out_consensus_fa = Path(out_consensus_fa)
    with out_catalog_csv.open("w", newline="") as cw, out_consensus_fa.open("w") as fw:
        w = csv.writer(cw)
        w.writerow(["locus", "allele_code", "sequence", "support_count", "total_support"])
        for loc in sorted(catalog_counts.keys()):
            cnts = catalog_counts[loc]
            if not cnts:
                continue
            total_support = sum(cnts.values())
            ranked = cnts.most_common()
            for idx, (seq, c) in enumerate(ranked, start=1):
                allele_code = 100 + idx  # 101, 102, ...
                w.writerow([loc, allele_code, seq, c, total_support])
            fw.write(f">{loc}\n{ranked[0][0]}\n")

# ---------------------------------------------------------------------
# helper: allele freq from genotype strings
# ---------------------------------------------------------------------
def _allele_freq_from_genos(per_sample_gt) -> dict[int, float]:
    """
    Accepts either:
      - dict[sample] -> "101102" or "000000" (skip no-calls)
      - dict[sample] -> {"status": "CALLED"/..., "gt": (code1, code2)}
    Returns allele-code frequencies (0..1) that sum to exactly 1.0 (if any calls), else {}.
    """
    from collections import Counter

    counts = Counter()
    n_alleles = 0

    for rec in (per_sample_gt or {}).values():
        # Case 1: string like "101102" or "000000"
        if isinstance(rec, str):
            s = rec.strip()
            if len(s) == 6 and s != "000000":
                a = int(s[:3]); b = int(s[3:])
                counts[a] += 1; counts[b] += 1
                n_alleles += 2
            continue

        # Case 2: dict with {"status": "...", "gt": (a,b)}
        if isinstance(rec, dict):
            if rec.get("status") != "CALLED":
                continue
            gt = rec.get("gt")
            if not gt or len(gt) != 2:
                continue
            try:
                a = int(gt[0]); b = int(gt[1])
            except Exception:
                # If codes came as strings that aren't ints, skip
                continue
            counts[a] += 1; counts[b] += 1
            n_alleles += 2

    if n_alleles == 0:
        return {}

    freqs = {k: v / float(n_alleles) for k, v in counts.items()}
    # Normalize to sum exactly 1.0 (avoid float drift)
    s = sum(freqs.values())
    if s and abs(s - 1.0) > 1e-12:
        # adjust the largest bin to close the gap
        kmax = max(freqs, key=freqs.get)
        freqs[kmax] += (1.0 - s)
    return freqs

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
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
    """
    Make the 2×2 dashboard for one locus.

    - Top-left: A1 vs A2 (HOMs colored/shaped by allele code; HET orange; NC black X; LOW red ▲)
    - Top-right: Depth vs %A2 (fixed 0–50 y, guide lines at 10 and 25)
    - Middle band: centered "Call rate | Background" text (keeps headroom for logo)
    - Bottom-left: Locus read share among panel (mean ±1 SD, red line = uniform share)
    - Bottom-right: Allele frequency bar chart (codes on x-axis, 0–1 y-axis)
    """
    if not HAVE_MPL:
        return

    from pathlib import Path
    import collections
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    outdir = Path(outdir) if outdir else Path("plots")
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------- figure layout (top headroom for logo) ----------
    fig = plt.figure(figsize=(12, 8.6))
    gs = GridSpec(2, 2, height_ratios=[1.0, 1.4], width_ratios=[1.0, 1.0],
                  hspace=0.35, wspace=0.28)
    fig.subplots_adjust(top=0.88)  # leave space for logo

    # ---------- styles ----------
    het_style = dict(marker="o", linestyle="none", alpha=0.9, color="C1", rasterized=True)
    nc_style  = dict(marker="x", linestyle="none", alpha=0.9, color="k",  rasterized=True)
    low_style = dict(marker="^", linestyle="none", alpha=0.9, color="red", rasterized=True)

    # HOM “lucky charms” mapping
    markers_cycle = ['o','s','D','^','v','P','X','*','<','>','h','H']
    colors_cycle  = ["#2ca02c","#d62728","#9467bd","#1f77b4","#ff7f0e","#17becf",
                     "#8c564b","#e377c2","#7f7f7f","#bcbd22","#aec7e8","#ffbb78"]
    # explicit preferred styles for common codes
    explicit = {101: ('o', "#2ca02c"), 102: ('s', "#d62728"), 103: ('D', "#9467bd")}

    # Which allele codes appear?
    present_codes = sorted(allele_freq_codes.keys()) if allele_freq_codes else sorted(
        {a for (_, _, c, a) in points_a1a2 if c == "HOM" and a is not None} |
        {a for (_, _, c, a) in depths_a2pct  if c == "HOM" and a is not None}
    )

    # Build a finite palette and cycle deterministically (no while loop; cannot hang)
    palette = [(m, c) for m in markers_cycle for c in colors_cycle]  # 12*12 = 144 combos
    code_style = {}
    for i, cd in enumerate(present_codes):
        code_style[cd] = palette[i % len(palette)]

    # Override with explicit styles, if present
    for cd, sty in explicit.items():
        if cd in code_style:
            code_style[cd] = sty

    def plot_hom(ax, xvals, yvals, codes):
        import collections
        groups = collections.defaultdict(lambda: ([], []))
        for x, y, cd in zip(xvals, yvals, codes):
            if cd is None:
                continue
            groups[int(cd)][0].append(x)
            groups[int(cd)][1].append(y)
        for cd in sorted(groups):
            m, c = code_style.get(cd, ('o', "#2ca02c"))
            xs, ys = groups[cd]
            ax.scatter(xs, ys, marker=m, c=c, edgecolors='none', linewidths=0, s=36, alpha=0.95, rasterized=True)

    # ---------- Top-left: A1 vs A2 ----------
    ax1 = fig.add_subplot(gs[0, 0])
    xs = [x for (x, y, call, code) in points_a1a2 if call == "HOM"]
    ys = [y for (x, y, call, code) in points_a1a2 if call == "HOM"]
    cs = [code for (x, y, call, code) in points_a1a2 if call == "HOM"]
    if xs: plot_hom(ax1, xs, ys, cs)
    xs = [x for (x, y, call, _) in points_a1a2 if call == "HET"]
    ys = [y for (x, y, call, _) in points_a1a2 if call == "HET"]
    if xs: ax1.plot(xs, ys, **het_style)
    xs = [x for (x, y, call, _) in points_a1a2 if call == "NC"]
    ys = [y for (x, y, call, _) in points_a1a2 if call == "NC"]
    if xs: ax1.plot(xs, ys, **nc_style)
    xs = [x for (x, y, call, _) in points_a1a2 if call == "LOW"]
    ys = [y for (x, y, call, _) in points_a1a2 if call == "LOW"]
    if xs: ax1.plot(xs, ys, **low_style)
    ax1.set_xlabel("A1 count"); ax1.set_ylabel("A2 count")
    ax1.set_title(f"{locus}: A1 vs A2")

    # ---------- Top-right: depth vs %A2 ----------
    ax2 = fig.add_subplot(gs[0, 1])
    xs = [d for (d, p, call, code) in depths_a2pct if call == "HOM"]
    ys = [p for (d, p, call, code) in depths_a2pct if call == "HOM"]
    cs = [code for (d, p, call, code) in depths_a2pct if call == "HOM"]
    if xs: plot_hom(ax2, xs, ys, cs)
    xs = [d for (d, p, call, _) in depths_a2pct if call == "HET"]
    ys = [p for (d, p, call, _) in depths_a2pct if call == "HET"]
    if xs: ax2.plot(xs, ys, **het_style)
    xs = [d for (d, p, call, _) in depths_a2pct if call == "NC"]
    ys = [p for (d, p, call, _) in depths_a2pct if call == "NC"]
    if xs: ax2.plot(xs, ys, **nc_style)
    xs = [d for (d, p, call, _) in depths_a2pct if call == "LOW"]
    ys = [p for (d, p, call, _) in depths_a2pct if call == "LOW"]
    if xs: ax2.plot(xs, ys, **low_style)
    ax2.axhline(a2_lo, ls="--", lw=1, color="gray")
    ax2.axhline(a2_hi, ls="--", lw=1, color="gray")
    ax2.set_ylim(0, 50)
    ax2.set_xlabel("Total locus depth"); ax2.set_ylabel("% A2")
    ax2.set_title(f"{locus}: depth vs %A2")

    # ---------- Middle band: call rate | background ----------
    if extra_info is None: extra_info = {}
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

    # ---------- Bottom-left: read distribution (mean ± SD) ----------
    ax3 = fig.add_subplot(gs[1, 0])
    means = [locus_pct_mean[L] for L in loci_sorted_by_share]
    sds   = [locus_pct_sd[L]   for L in loci_sorted_by_share]
    xs = list(range(len(loci_sorted_by_share)))
    n_loci = max(1, len(loci_sorted_by_share))

    # Detect units robustly: if the sum is near 1, treat as fractions; else percents
    total = sum(means) if means else 0.0
    if total <= 5.0:
        means = [m * 100.0 for m in means]
        sds   = [sd * 100.0 for sd in sds]
    avg_share = 100.0 / n_loci

    bars = ax3.bar(xs, means, color="#b3b3b3", edgecolor="#b3b3b3", linewidth=0.0, zorder=1)
    if locus in loci_sorted_by_share:
        idx = loci_sorted_by_share.index(locus)
        bars[idx].set_facecolor("tab:red"); bars[idx].set_edgecolor("tab:red"); bars[idx].set_zorder(2)
    if sds:
        ax3.errorbar(xs, means, yerr=sds, fmt="none",
                     ecolor="k", elinewidth=0.8, capthick=0.8, capsize=1, zorder=3)
    ax3.axhline(avg_share, color="tab:red", lw=1, zorder=5)
    ymax = max((m + sd) for m, sd in zip(means or [0.0], sds or [0.0]))
    ax3.set_ylim(0, max(ymax * 1.15, avg_share * 1.8, 2.0))
    ax3.set_xlim(-0.5, len(xs) - 0.5 + 1.0)
    ax3.set_ylabel("Avg % reads per locus (±1 SD)")
    ax3.set_xlabel("Loci (sorted by mean %, least → most)")
    ax3.set_title("Read distribution among loci")
    if len(xs) > 20:
        ax3.set_xticks([])

    # ---------- Bottom-right: allele frequency (0..1, must sum to 1 when present) ----------
    ax4 = fig.add_subplot(gs[1, 1])
    freqs_dict = allele_freq_codes or {}
    codes = sorted(freqs_dict.keys())
    freqs = [freqs_dict[c] for c in codes]
    if codes:
        # final nudge to avoid 0.999999 due to float drift
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
    ax4.set_title(f"{locus}: allele frequency (101→)")

    # ---------- save ----------
    outpath = outdir / f"{locus}_dashboard.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)
    
class _Timeout(Exception):
    pass

def _timeout_handler(signum, frame):
    raise _Timeout("plotting timed out")

def run_with_timeout(func, timeout_s, *args, **kwargs):
    """Run func(*args, **kwargs) with a SIGALRM timeout (Linux/Unix only)."""
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    try:
        signal.alarm(max(1, int(timeout_s)))
        return func(*args, **kwargs)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

def write_placeholder_png(outpath, locus, reason):
    """Write a small placeholder image indicating why this locus was skipped."""
    # use matplotlib since it's already a dependency and Agg backend is set
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 2.2))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    msg = f"{locus}\nPlot skipped: {reason}\nGenerated {ts}"
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11)
    fig.savefig(outpath, dpi=140, bbox_inches="tight")
    plt.close(fig)

# ----------------------------------------------------------------------
# SECOND PASS: re-read resolved_*.fastq.gz but ONLY count catalog alleles
# ----------------------------------------------------------------------
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
    # load catalog locus->seq->allele_code
    catalog: dict[str, dict[str,int]] = defaultdict(dict)
    with Path(catalog_csv).open() as fh:
        rdr = csv.DictReader(fh)
        cols = {k.lower(): k for k in rdr.fieldnames}
        for row in rdr:
            loc  = row[cols["locus"]]
            seq  = row[cols["sequence"]]
            code = int(row[cols["allele_code"]])
            catalog[loc][seq] = code

    resolved_paths = sorted(Path(resolved_dir).glob("resolved_*.fastq.gz"))

    per_locus_sample_gt: dict[str, dict[str, str]] = defaultdict(dict)
    a2_rows = []
    locus_points_a1a2 = defaultdict(list)
    locus_points_depth_pct = defaultdict(list)

    locus_called_samples = defaultdict(int)
    locus_total_samples  = defaultdict(int)
    locus_background_sum = defaultdict(float)

    for fpath in resolved_paths:
        sample = fpath.stem.replace("resolved_", "")
        per_sample_background_frac = defaultdict(list)

        with gzip.open(fpath, "rt") as fh:
            from collections import Counter, defaultdict as dd
            sample_locus_cnts: dict[str, Counter] = dd(Counter)
            sample_locus_offcat: dict[str, int]  = dd(int)

            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().strip()
                _ = fh.readline(); _ = fh.readline()
                loc = h.strip().split("|locus=")[-1] if "|locus=" in h else None
                if not loc:
                    continue
                if loc in catalog and seq in catalog[loc]:
                    sample_locus_cnts[loc][seq] += 1
                else:
                    sample_locus_offcat[loc] += 1

        loci_iter = list(locus_universe) if locus_universe else list(catalog.keys())
        for loc in loci_iter:
            cnts = sample_locus_cnts.get(loc, Counter())
            total_cat = sum(cnts.values())
            off_cat   = sample_locus_offcat.get(loc, 0)

            if total_cat > 0 or off_cat > 0:
                locus_total_samples[loc] += 1

            if total_cat == 0:
                per_locus_sample_gt[loc][sample] = "000000"
                a2_rows.append({"sample": sample, "locus": loc, "depth": 0, "A1_count": 0, "A2_count": 0, "A2_pct": 0.0, "call": "LOW"})
                locus_points_a1a2[loc].append((0, 0, "LOW", None))
                locus_points_depth_pct[loc].append((0, 0.0, "LOW", None))
                if (total_cat + off_cat) > 0:
                    per_sample_background_frac[loc].append(off_cat / (total_cat + off_cat))
                continue

            top = cnts.most_common(2)
            (a1_seq, a1_ct) = top[0]
            a1_code = catalog[loc][a1_seq]
            if len(top) > 1:
                (a2_seq, a2_ct) = top[1]; a2_code = catalog[loc][a2_seq]
            else:
                a2_seq, a2_ct, a2_code = None, 0, None

            denom = a1_ct + a2_ct
            if denom < min_depth:
                call = "LOW"; gt = "000000"; a2_pct = 0.0
            else:
                frac_a2 = a2_ct / denom
                if a2_code is not None and frac_a2 >= 0.25:
                    call = "HET"; gt = f"{min(a1_code, a2_code)}{max(a1_code, a2_code)}"
                elif (a1_ct / denom) >= 0.90:
                    call = "HOM"; gt = f"{a1_code}{a1_code}"
                else:
                    call = "NC";  gt = "000000"
                a2_pct = frac_a2 * 100.0

            if call in ("HOM","HET"):
                locus_called_samples[loc] += 1
            per_locus_sample_gt[loc][sample] = gt
            a2_rows.append({"sample": sample, "locus": loc, "depth": denom, "A1_count": a1_ct, "A2_count": a2_ct, "A2_pct": round(a2_pct,2), "call": call})
            locus_points_a1a2[loc].append((a1_ct, a2_ct, call, a1_code if call=="HOM" else None))
            locus_points_depth_pct[loc].append((denom, a2_pct, call, a1_code if call=="HOM" else None))

            other_catalog = total_cat - denom
            total_reads_for_loc = off_cat + total_cat
            if total_reads_for_loc > 0:
                per_sample_background_frac[loc].append((off_cat + other_catalog) / total_reads_for_loc)

        for loc, fracs in per_sample_background_frac.items():
            if fracs:
                locus_background_sum[loc] += (sum(fracs) / len(fracs))

    # write per-sample A1/A2 metrics
    with Path(f"{out_prefix}_a2_metrics.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample","locus","depth","A1_count","A2_count","A2_pct","call"])
        for r in a2_rows:
            w.writerow([r["sample"], r["locus"], r["depth"], r["A1_count"], r["A2_count"], r["A2_pct"], r["call"]])

    # genotype matrix: write all primer loci, even if never in catalog
    all_samples = sorted({s for d in per_locus_sample_gt.values() for s in d})
    loci_to_write = sorted(locus_universe) if locus_universe else sorted(catalog.keys())
    with Path(f"{out_prefix}_genotypes.tsv").open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["locus"] + all_samples)
        for loc in loci_to_write:
            d = per_locus_sample_gt.get(loc, {})
            w.writerow([loc] + [d.get(s, "000000") for s in all_samples])

    # call-rate / background for headers
    locus_call_rate = {}
    locus_background_pct = {}
    for loc in loci_to_write:
        tot = locus_total_samples.get(loc, 0)
        called = locus_called_samples.get(loc, 0)
        locus_call_rate[loc] = (called / tot * 100.0) if tot else 0.0
        bg_sum = locus_background_sum.get(loc, 0.0)
        locus_background_pct[loc] = (bg_sum / tot * 100.0) if tot else 0.0

    # plots for *all* primer loci
    if plots_dir and HAVE_MPL:
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Pad zero-signal loci so they show up as 0-height bars in the bottom-left plot
        padded_stats = _augment_locus_stats_with_zeros(locus_stats, locus_universe)
        loci_sorted_by_share = padded_stats["sorted_loci_by_mean"]
        locus_pct_mean = padded_stats["mean_pct"]
        locus_pct_sd   = padded_stats["sd_pct"]
        try:
            from tqdm import tqdm as _tqdm
            itr = _tqdm(loci_to_write, desc="Plotting loci")
        except Exception:
            itr = loci_to_write

        for loc in itr:
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
                )
            except Exception as e:
                try:
                    write_placeholder_png(Path(plots_dir)/f"{loc}_dashboard.png", loc, f"error: {e}")
                except Exception:
                    pass
            
def compute_locus_readshare_stats_from_resolved(resolved_dir: Path) -> dict:
    """
    Recompute locus read-share stats (mean % and SD) from already-resolved FASTQs.

    For each resolved_<sample>.fastq.gz:
      - count reads per locus for that sample
      - convert to % of that sample's resolved reads
      - append to locus_depths[locus]
    Returns:
      {
        "sorted_loci_by_mean": [locus,...]  # ascending by mean %
        "mean_pct": {locus: mean_percent},
        "sd_pct":   {locus: sd_percent},
      }
    """
    import gzip
    from collections import defaultdict

    # Reuse the global progress helper if available
    def _iter_with_progress(it, total):
        try:
            from tqdm import tqdm as _tqdm
            return _tqdm(it, desc="Scanning resolved FASTQs", unit="sample", total=total)
        except Exception:
            return it

    # Find resolved FASTQs (skip partials)
    resolved_paths = sorted(
        p for p in resolved_dir.glob("resolved_*.fastq.gz")
        if not str(p).endswith(".part")
    )

    locus_depths: dict[str, list[float]] = defaultdict(list)

    for fpath in _iter_with_progress(resolved_paths, total=len(resolved_paths)):
        per_locus_counts: dict[str, int] = defaultdict(int)

        with gzip.open(fpath, "rt") as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline()   # sequence
                _   = fh.readline()   # '+'
                _   = fh.readline()   # qualities

                # Header is like "...|locus=<locus>\n"
                if "|locus=" not in h:
                    continue
                loc = h.strip().split("|locus=")[-1]
                if not loc:
                    continue
                per_locus_counts[loc] += 1

        total_resolved = sum(per_locus_counts.values())
        if total_resolved > 0:
            for loc, ct in per_locus_counts.items():
                pct = (ct / total_resolved) * 100.0
                locus_depths[loc].append(pct)

    # Delegate to the existing reducer to compute mean/sd + sorted order
    return compute_locus_readshare_stats(locus_depths)

# --- MP helper: call process_sample_pair *positionally* to avoid kwarg-name drift
def _mp_process_sample_pair(args_tuple):
    """
    args_tuple = (sample_id, r1_path, r2_path,
                  primers_csv_path, min_amplicon, max_amplicon,
                  min_overlap, max_ov_mismatch_frac)
    Returns whatever process_sample_pair returns (or raises).
    """
    return process_sample_pair(*args_tuple)
    
def _augment_locus_stats_with_zeros(locus_stats: dict, locus_universe) -> dict:
    """
    Ensure every locus in locus_universe appears in locus_stats with mean=0, sd=0,
    and rebuild the sorted order used by the bottom-left panel.
    """
    mean = dict(locus_stats.get("mean_pct", {}))
    sd   = dict(locus_stats.get("sd_pct", {}))

    if locus_universe:
        for loc in locus_universe:
            if loc not in mean:
                mean[loc] = 0.0
                sd[loc]   = 0.0

    sorted_by_mean = sorted(mean.keys(), key=lambda l: mean[l])
    return {"mean_pct": mean, "sd_pct": sd, "sorted_loci_by_mean": sorted_by_mean}

# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    import argparse
    from pathlib import Path
    import sys
    from multiprocessing import cpu_count

    p = argparse.ArgumentParser(
        description="GT-seq microhap: parallel resolve → strict catalog → genotypes → dashboards"
    )
    p.add_argument("-i", "--input-dir", required=True,
                   help="directory with paired fastqs (R1/R2, .fastq or .fastq.gz)")
    p.add_argument("-p", "--primers", required=True,
                   help="CSV with columns: locus,fwd_primer,rev_primer (extra cols OK)")
    p.add_argument("-o", "--out-prefix", default="microhap",
                   help="output prefix (default: microhap)")
    p.add_argument("--min-amplicon", type=int, default=60)
    p.add_argument("--max-amplicon", type=int, default=250)
    p.add_argument("--min-depth", type=int, default=10,
                   help="min per-sample per-locus depth to call genotype (A1+A2)")
    p.add_argument("--min-overlap", type=int, default=12)
    p.add_argument("--max-ov-mismatch-frac", type=float, default=0.05)
    p.add_argument("--plots", dest="plots_dir", default=None,
                   help="directory for per-locus dashboards")
    p.add_argument("--force-resolve", action="store_true",
                   help="rebuild resolved_fastqs/ even if .fastq.gz already exist")
    p.add_argument("--a2-lo", type=float, default=10.0)
    p.add_argument("--a2-hi", type=float, default=25.0)
    p.add_argument("--threads", type=int, default=None,
                   help="worker processes for resolve pass (default: half of CPUs)")

    args = p.parse_args()

    input_dir   = Path(args.input_dir)
    primers_csv = Path(args.primers)
    out_prefix  = args.out_prefix
    plots_dir   = Path(args.plots_dir) if args.plots_dir else None

    # ---------- load primers / pairs ----------
    primers = load_primers(primers_csv)
    pairs = find_fastq_pairs(input_dir)  # {sample: (r1, r2)}
    print(f"Found {len(pairs)} paired samples.")

    resolved_dir = Path("resolved_fastqs")
    resolved_dir.mkdir(parents=True, exist_ok=True)

    # Determine worker count (default: half the CPUs, at least 1)
    if args.threads is None:
        try:
            half = max(1, (cpu_count() or 2) // 2)
        except Exception:
            half = 1
        workers = half
    else:
        workers = max(1, int(args.threads))

    # Determine if we must (re)resolve
    existing = sorted(
        p for p in resolved_dir.glob("resolved_*.fastq.gz")
        if not str(p).endswith(".part")
    )
    do_first_pass = args.force_resolve or (len(existing) == 0)

    if do_first_pass:
        # ---------- First pass: resolve in parallel ----------
        # NOTE: resolve_fastqs_multiprocess expects tuples (r1, r2, sample)
        paired_samples = [(r1, r2, sample) for sample, (r1, r2) in pairs.items()]
        done, errors = resolve_fastqs_multiprocess(
            paired_samples=paired_samples,
            resolved_dir=resolved_dir,
            min_amplicon=args.min_amplicon,
            max_amplicon=args.max_amplicon,
            min_overlap=args.min_overlap,
            max_ov_mismatch_frac=args.max_ov_mismatch_frac,
            primers=primers,
            force=True,
            workers=workers,
        )
        print(f"Resolved {done} samples into {resolved_dir}/")
        if errors:
            print(f"[WARN] {len(errors)} worker errors (showing first 3):")
            for e in errors[:3]:
                print("  -", e, file=sys.stderr)

        # sanity check: ensure we actually have outputs
        produced = sorted(
            p for p in resolved_dir.glob("resolved_*.fastq.gz")
            if not str(p).endswith(".part")
        )
        if not produced:
            sys.exit("ERROR: parallel resolve produced no files. See warnings above.")

        # Locus read-share stats (bottom-left panel) from resolved files
        locus_stats = compute_locus_readshare_stats_from_resolved(resolved_dir)

        # Strict catalog from resolved FASTQs (top-2 with thresholds)
        catalog_path   = Path(f"{out_prefix}_catalog.csv")
        consensus_path = Path(f"{out_prefix}_consensus.fasta")
        build_catalog_from_resolved(
            resolved_dir=resolved_dir,
            out_catalog_csv=catalog_path,
            out_consensus_fa=consensus_path,
            hom_min=0.85,  # ≥85% = HOM evidence
            het_min=0.30,  # both ≥30% = HET evidence
        )

    else:
        print(f"Found {len(existing)} resolved FASTQs in {resolved_dir}/; "
              f"skipping first pass (use --force-resolve to rebuild).")

        # Need a catalog to proceed if we skipped resolving
        catalog_path = Path(f"{out_prefix}_catalog.csv")
        if not catalog_path.exists():
            sys.exit(f"ERROR: {catalog_path} not found. Re-run with --force-resolve to (re)build the catalog.")

        # Compute locus stats from existing resolved files
        locus_stats = compute_locus_readshare_stats_from_resolved(resolved_dir)

    # ---------- Second pass: catalog-only counting, calls, plots ----------
    # ensure we use the path just built/checked
    catalog_path = Path(f"{out_prefix}_catalog.csv")
    locus_universe = set(primers.keys())  # plot every primer locus, even if zero signal

    # create plots dir if requested
    if plots_dir:
        plots_dir.mkdir(parents=True, exist_ok=True)

    call_from_resolved_fastqs(
        resolved_dir=resolved_dir,
        catalog_csv=catalog_path,
        out_prefix=out_prefix,
        min_depth=args.min_depth,
        locus_stats=locus_stats,
        plots_dir=plots_dir,
        a2_lo=args.a2_lo,
        a2_hi=args.a2_hi,
        locus_universe=locus_universe,
    )

    # ---------- Footer ----------
    print("Wrote:")
    print("  resolved_fastqs/resolved_<sample>.fastq.gz")
    print(f"  {out_prefix}_a2_metrics.tsv")
    print(f"  {out_prefix}_catalog.csv")
    print(f"  {out_prefix}_consensus.fasta")
    print(f"  {out_prefix}_genotypes.tsv")
    if plots_dir:
        print(f"  {plots_dir}/<locus>_dashboard.png")

if __name__ == "__main__":
    main()

