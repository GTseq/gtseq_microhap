
# gtseq_microhap

Direct microhaplotype cataloging and genotyping for GT-seq paired-end amplicon data.

`gtseq_microhap_catalog_and_call.py` resolves primer-bounded paired-end reads, builds a microhaplotype allele catalog, calls diploid genotypes, exports phased microhaplotypes, optional single-SNP matrices, accessory assay calls, publication-quality QC plots, and reference-anchored VCF output.

**Current release:** **v0.3.1**

---

## Features

- Quality-aware paired-end amplicon reconstruction
- Diploid-aware catalog construction
- Catalog-based genotype calling (HOM/HET/NC/LOW)
- Phased SNP-only and SNP+indel outputs
- Optional single-SNP exports
- SexID and parasite/presence/absence assays
- Reference-anchored VCF export
- Publication-quality QC dashboards

---

## Citation

If you use **gtseq_microhap** in published work, please cite the GT-seq microhap preprint:

**Preprint DOI:** https://doi.org/10.64898/2026.04.01.715880

Please also cite:

Campbell NR, Harmon SA, Narum SR. 2015. *Genotyping-in-Thousands by sequencing (GT-seq): A cost effective SNP genotyping method based on custom amplicon sequencing.* Molecular Ecology Resources 15:855–867.

---

## Example dataset and reproducibility resources

The Delta Smelt validation dataset, BWA-generated reference VCF, validation scripts, concordance analyses, PIC analyses, and supporting files are archived on Zenodo:

https://doi.org/10.5281/zenodo.15780159

---

## Requirements

- Python 3.8+
- matplotlib
- tqdm (optional)
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

* Paired-end FASTQ/FASTQ.GZ files
* Primer CSV/TSV containing `locus`, `fwd_primer`, and `rev_primer`
* Optional SexID and parasite marker lists

---

## Quick start

```bash
python3 gtseq_microhap_catalog_and_call.py   --primers primers.csv   --indir R1_R2   --outdir microhap_out
```

---

## Main outputs

All outputs are written to `--outdir` unless otherwise noted.

| Output | Description |
|---------|-------------|
| `sample_stats.tsv` | Raw read pair counts and primer-bounded read counts for every sample. |
| `resolved_fastqs/resolved_<sample>.fastq.gz` | Primer-bounded resolved amplicon FASTQs used for downstream analysis. |
| `microhap_catalog.csv` | Final allele catalog containing locus, allele code, sequence, support count, and total support. |
| `microhap_catalog_length_pruned.tsv` | Candidate alleles removed by length/divergence filtering. Provides an audit trail even when no alleles are removed. |
| `microhap_consensus.fasta` | Highest-supported consensus allele sequence for each locus. |
| `microhap_a2_metrics.tsv` | Per-sample A1/A2 counts, allele codes, depth, A2 percentage, and genotype call. |
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
