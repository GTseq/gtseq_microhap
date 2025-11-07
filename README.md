# gtseq_microhap

Identifies and genotypes microhaplotypes from GT-seq paired-end FASTQ data using a `primers.csv` file.  
Builds a catalog of unique alleles, calls genotypes across samples, and generates per-locus QC dashboards.  
Alignment-free except for short R1/R2 overlap stitching.

---

## Features

- Parallelized **FASTQ resolving** from raw GT-seq paired-end reads.
- Strict, data-driven **microhaplotype catalog** construction:
  - Uses only strong per-sample evidence (top-2 alleles with configurable thresholds).
  - Filters background / spurious alleles.
- **Catalog-only second pass**:
  - Counts only catalog alleles.
  - Calls genotypes using A1/A2 rules (HOM / HET / NC / LOW).
  - Allele frequencies computed from called genotypes only.
- **Per-locus dashboards** (`plots/<locus>_dashboard.png`):
  - A1 vs A2 read counts.
  - Total depth vs %A2.
  - Call rate and background summary.
  - Read distribution among loci (panel-wide on-target balance).
  - Allele frequency barplot per locus.
- Handles loci with many alleles (true microhaplotypes), plus loci with zero signal.
- Designed for multi-sample GT-seq panels (hundreds of loci, hundreds of individuals).

---

## Quick start

### Inputs

1. **Paired-end FASTQs** in one directory (e.g. `R1_R2/`), named in R1/R2 pairs.
2. **Primers CSV** with at least:
   ```text
   locus,fwd_primer,rev_primer
   NC_XXXXXX_1_ABC123,ACGT...,TGCA...
   ...

### Basic usage

python3 gtseq_microhap_catalog_and_call.py \
  -i ./R1_R2 \
  -p ./primers.csv \
  --plots plots

This will:

1- Resolve reads into resolved_fastqs/resolved_<sample>.fastq.gz

2- Build a strict allele catalog (microhap_catalog.csv, microhap_consensus.fasta)

3- Call genotypes (microhap_genotypes.tsv, microhap_a2_metrics.tsv)

4- Generate dashboards (plots/<locus>_dashboard.png)

### Important options

--min-amplicon / --max-amplicon
Acceptable stitched amplicon length.

--min-depth
Minimum (A1 + A2) reads at a locus in a sample to attempt a call.

--a2-lo / --a2-hi
%A2 thresholds for HOM / NC / HET decision.

--force-resolve
Rebuild resolved_fastqs/ from raw FASTQs.

--threads
Worker processes for the resolve step (default: half your CPUs).

Run with -h to see full options.

### Requirements
Python 3.8+
matplotlib
tqdm

### Citation / Authors
Developed by GTseek for GT-seq-based microhaplotype analysis.
If you use this tool in a publication, please cite the GT-seq method paper (Campbell et al. 2015) and this repository.
