
# gtseq_microhap

Direct microhaplotype cataloging and genotyping for GT-seq paired-end and single-end amplicon data.

`gtseq_microhap_catalog_and_call.py` resolves primer-bounded paired-end reads or probe-validated single-end reads, builds a microhaplotype allele catalog, calls diploid genotypes, exports phased microhaplotypes, optional single-SNP matrices, accessory assay calls, publication-quality QC plots, and reference-anchored VCF output.

**Current release:** **v0.3.4**

---

## Features

- Quality-aware paired-end amplicon reconstruction (default)
- Single-end resolution using forward primers and internal probe sequences (`--single-end`)
- Diploid-aware catalog construction
- Catalog-based genotype calling (HOM/HET/NC/LOW)
- Per-sample, per-locus catalog allele depth and locus-fraction diagnostics
- Phased SNP-only and SNP+indel outputs
- Optional single-SNP exports
- SexID and parasite/presence/absence assays
- Reference-anchored VCF export
- Publication-quality QC dashboards

---

## Citation

If you use **gtseq_microhap** in published work, please cite:

Campbell NR, Campbell AR, Blair SK, Finger AJ. 2026. *Direct Microhaplotype Genotyping for GT-seq (Genotyping-in-Thousands by Sequencing) Using a Diploid Abundance Model.* PLOS Computational Biology. https://doi.org/10.1371/journal.pcbi.1014808

Please also cite:

Campbell NR, Harmon SA, Narum SR. 2015. *Genotyping-in-Thousands by sequencing (GT-seq): A cost effective SNP genotyping method based on custom amplicon sequencing.* Molecular Ecology Resources 15:855–867.

---

## Example dataset and reproducibility resources

The Delta Smelt validation dataset, BWA-generated reference VCF, validation scripts, concordance analyses, PIC analyses, and supporting files are archived on Zenodo:

https://doi.org/10.5281/zenodo.19069550

---

## Requirements

- Python 3.10+
- matplotlib
- tqdm (required for single-end progress reporting; optional for paired-end)
- bwa (required only for `--emit-vcf`)

---

## Installation

```bash
git clone https://github.com/GTseq/gtseq_microhap.git
cd gtseq_microhap
python3 gtseq_microhap_catalog_and_call.py --help
```

---

## Inputs

### Paired-end (default)

Provide paired R1/R2 FASTQ or FASTQ.GZ files in `--indir`, named so that the pipeline can identify matching R1 and R2 files. Provide a primer CSV/TSV containing `locus`, `fwd_primer`, and `rev_primer` (aliases `fwd` and `rev` are accepted). Paired-end reads are merged and trimmed to retain the complete forward-primer-through-reverse-complement-of-reverse-primer amplicon.

### Single-end (`--single-end`)

Provide single-end FASTQ or FASTQ.GZ files in `--indir` and a **CSV** with five columns:

```csv
locus,fwd_primer,rev_primer,probe1,probe2
Locus001,ACGTACGT,TGCATGCA,AGC.TGAC,AGC.TGTC
```

`fwd` and `rev` are also accepted as header names. These example sequences are illustrative, not a validated assay. Each locus requires a forward primer, reverse primer, and at least one nonempty probe.

The SE resolver requires a read to start with the locus's forward primer and contain **at least one** downstream match to `probe1`, `probe2`, or the reverse complement of either probe. In probe patterns, `.` matches exactly one nucleotide (`A`, `C`, `G`, or `T`); the same wildcard rule applies after reverse complementation. Probe sequences are matching patterns, not replacements for observed bases.

If the reverse complement of the reverse primer (`RC-REV`) occurs downstream, the resolver trims **after** the complete RC-REV sequence. Both primer sequences remain in the resolved read. If RC-REV is absent, the read is retained to its sequenced endpoint, subject to the amplicon-length limits. Ambiguous matches to multiple loci are discarded. R2 files and existing `resolved_` FASTQs are excluded from SE input discovery; multiple SE files collapsing to the same sample name must be merged beforehand.

**Important limitation:** Reads from the same biological allele that end at different positions are distinct full-length sequences to the unchanged catalog. For reliable SE genotyping, use read lengths that consistently span the target or evaluate potential terminal-length artifacts and background fractions.

### Optional inputs

SexID and parasite marker lists may also be supplied using the existing command-line options.

---

## Quick start

**Paired-end (default):**

```bash
python3 gtseq_microhap_catalog_and_call.py \
  --primers primers.csv \
  --indir R1_R2 \
  --outdir microhap_PE
```

**Single-end:**

```bash
python3 gtseq_microhap_catalog_and_call.py \
  --single-end \
  --primers targets_with_probes.csv \
  --indir SE_FASTQs \
  --outdir microhap_SE \
  --workers 20
```

Single-end mode displays a combined read-count progress bar with `tqdm`. To estimate total reads and remaining time, it makes an initial counting pass over the input FASTQs. The resolver produces the same `resolved_fastqs/resolved_<sample>.fastq.gz` interface used by the existing catalog and genotype caller. Keep separate output directories for PE and SE runs to avoid mixing resolved FASTQs.

---

## Main outputs

All outputs are written to `--outdir` unless otherwise noted.

| Output | Description |
|---------|-------------|
| `sample_stats.tsv` | Raw read counts (read pairs for PE; single reads for SE) and accepted resolved-read counts for every sample. |
| `resolved_fastqs/resolved_<sample>.fastq.gz` | Resolved amplicon FASTQs used for downstream analysis; SE reads may lack RC-REV if sequencing ended before the reverse primer. |
| `microhap_catalog.csv` | Final allele catalog containing locus, allele code, sequence, support count, and total support. |
| `microhap_catalog_length_pruned.tsv` | Candidate alleles removed by length/divergence filtering. Provides an audit trail even when no alleles are removed. |
| `microhap_consensus.fasta` | Highest-supported consensus allele sequence for each locus. |
| `microhap_a2_metrics.tsv` | Per-sample A1/A2 counts, allele codes, depth, A2 percentage, and genotype call. |
| `microhap_allele_depths.tsv` | Long-format per-sample, per-locus read depths and locus fractions for every observed final-catalog allele. `BACKGROUND` rows report primer-bounded reads that do not match a final catalog allele; allele and background fractions use all primer-bounded reads at that sample × locus as the denominator. |
| `microhap_genotypes.tsv` | Locus-by-sample matrix of allele-coded diploid genotypes. Missing/no-calls are `000000`. |
| `microhap_phased_def.tsv` | Variant definitions used to generate phased haplotypes. |
| `microhap_phased.tsv` | Per-sample phased SNP-only and SNP+indel haplotypes. |
| `microhap_phasedVARs.tsv` | Compact phased SNP+indel haplotype matrix. |
| `microhap_gt_vs_raw.png` | Genotype call rate versus raw read count across samples. |
| `microhap_gt_vs_pb.png` | Genotype call rate versus primer-bounded read count. |
| `microhap_raw_reads_hist.png` | Histogram of raw sequencing reads per sample. |
| `microhap_pb_reads_hist.png` | Histogram of primer-bounded reads per sample. |
| `microhap_on_target_fraction.png` | Distribution of on-target fractions across samples. |
| `plots/` | Per-locus QC dashboards including A1/A2 behavior, depth, allele frequencies, background, and call rate. |

## Optional outputs

| Output | Trigger | Description |
|---------|---------|-------------|
| `microhap_singleSNP_first.tsv` | `--single-snp-output first/both` | First biallelic SNP from each locus. |
| `microhap_singleSNP_highestMAF.tsv` | `--single-snp-output highest_maf/both` | Highest-MAF SNP from each locus. |
| `microhap_singleSNP_metadata.tsv` | Any SNP export | Metadata describing SNP selection. |
| `sexid_calls.tsv` | `--sexid-markers` | Sample-level sex assignments. |
| `sexid_marker_pot.tsv` | `--sexid-markers` | Marker-level pOT summaries. |
| `sexid_sample_marker_scores.tsv` | `--sexid-markers` | Per-sample marker scores. |
| `parasite_calls.tsv` | `--parasite-markers` | Presence/absence calls. |
| `parasite_summary.tsv` | `--parasite-markers` | Marker summary statistics. |
| `parasite_sample_marker_scores.tsv` | `--parasite-markers` | Per-sample parasite scores. |
| `microhap_variants.vcf` | `--emit-vcf` | Reference-anchored VCF. |
| `microhap_vcf_locus_placements.tsv` | `--emit-vcf` | BWA placement of loci. |
| `microhap_vcf_variant_map.tsv` | `--emit-vcf` | Mapping of catalog variants to VCF records. |
| `microhap_vcf_excluded_loci.tsv` | `--emit-vcf` | Excluded loci and reasons. |

---

## Command-line documentation

Run:

```bash
python3 gtseq_microhap_catalog_and_call.py --help
```

The help output documents every option and its default value.

---

## Version history

### v0.3.4

- Added optional `--single-end` mode with five-column locus/primer/probe CSV input.
- SE locus validation requires the forward primer plus either internal probe in either orientation, including `.` single-nucleotide wildcard patterns.
- Preserves observed forward and reverse-primer sequences and trims adapter read-through after RC-REV when present.
- Added multiprocess, read-based SE resolution progress reporting with ETA.
- PE processing remains the default; downstream catalog construction and genotype calling are unchanged.

### v0.3.3

- Added `microhap_allele_depths.tsv` with per-sample, per-locus allele depths, fractions, rank, and off-catalog `BACKGROUND` reads (requested for locus QC).
- Retained the existing diploid A1/A2 genotype inference.

### v0.3.2 (previous documented release)

- Previous paired-end release with catalog length/divergence filtering and diploid genotyping.
- The prior README attributed allele-depth reporting to v0.3.2; it is described under v0.3.3 above for the intended release sequence.

### v0.3.1

- Reference-anchored VCF export
- Phased SNP and SNP+indel outputs
- Catalog length/divergence filtering
- Single-SNP exports
- SexID and parasite assays
- Publication-ready documentation

---

## Author

Nathan R. Campbell  
GTseek LLC
