# Outputs

This document describes the output directory structure produced by the pipeline. All paths are relative to the run's output directory (e.g. `ChinaTown/` in the examples below).

## `fastqc/`

FastQC reports generated on raw input reads, before trimming.

| Path | Description |
|---|---|
| `fastqc/html/<sample>_fastqc.html` | Per-sample FastQC HTML report |
| `fastqc/zip/<sample>_fastqc.zip` | Per-sample FastQC data archive |

---

## `multiqc/`

Aggregated QC report combining FastQC (raw and trimmed), Cutadapt, and STAR alignment logs.

| Path | Description |
|---|---|
| `multiqc/multiqc_report.html` | Combined QC summary across all samples |

---

## `star/`

### `star/index/`

| Path | Description |
|---|---|
| `star/index/star` | STAR genome index directory |

### `star/align/log/`

Per-sample STAR alignment logs.

| Path | Description |
|---|---|
| `star/align/log/<sample>.Log.final.out` | Final alignment summary statistics |
| `star/align/log/<sample>.Log.out` | Full STAR run log |
| `star/align/log/<sample>.Log.progress.out` | Alignment progress log |

### `star/align/sortedbam/`

| Path | Description |
|---|---|
| `star/align/sortedbam/<sample>.Aligned.sortedByCoord.out.bam` | Coordinate-sorted aligned reads, before deduplication |

---

## `dedupped/`

| Path | Description |
|---|---|
| `dedupped/<sample>_converted.sam` | Deduplicated reads (via UMI-tools if `--umi true`, otherwise `samtools markdup`), converted to SAM for downstream site-counting |

---

## `temple/`

Outputs from the editing-site calling pipeline (the `temple`/`sam_to_map`/`find_edit_sites`/`combine_edit_sites` toolchain). Named after the internal Rust crate.

### `temple/map/`

| Path | Description |
|---|---|
| `temple/map/<sample>_map.bin` | Per-sample per-position base counts (`sam_to_map` output). Binary format: one record per covered position, tracking A->G and T->C mismatch/match counts separately, to support strand-aware editing detection downstream. |

### `temple/pooled_controls/`

| Path | Description |
|---|---|
| `temple/pooled_controls/<condition>.bin` | Control condition's per-replicate maps pooled into a single binary file (via `--control_combination`), used as the reference/background for site-calling. Only produced for conditions used as a `control` in `comparisons_csv`. |

### `temple/filtered_replicate_bins/`

| Path | Description |
|---|---|
| `temple/filtered_replicate_bins/<treatment>_vs_<control>.bin` | Per-replicate candidate editing sites that passed all `find_edit_sites` thresholds, one file per replicate, before cross-replicate combination. |

### `temple/combined_tsvs/<combination_mode>/`

| Path | Description |
|---|---|
| `temple/combined_tsvs/<combination_mode>/<treatment>_vs_<control>.tsv` | Per-replicate filtered bins merged across replicates per `--combination_mode` (`and`/`or`). Columns: `chr, start, end, rna_G, rna_total, ctrl_G, ctrl_total, rna_edit_frac`. This is the editing-site table after background-subtraction (if a `background` condition was specified in `comparisons_csv`) but before gene annotation. |

### `temple/annotated_edit_sites/`

| Path | Description |
|---|---|
| `temple/annotated_edit_sites/<treatment>_vs_<control>_annotated.tsv` | Combined editing sites intersected against the synthetic gene-span GTF (`bedtools intersect -wa -wb`). One row per (site, overlapping gene) pair , a site inside multiple overlapping genes produces multiple rows. Columns are the site columns above, followed by the standard 9 GTF columns (`gtf_chr, gtf_source, gtf_feature, gtf_start, gtf_end, gtf_score, gtf_strand, gtf_frame, gtf_attributes`). |

### `temple/summarized_edit_sites/<edit_threshold>/`

Final gene and transcript-level summary tables, filtered to sites with `rna_edit_frac >= --edit_threshold`. The directory name matches the `--edit_threshold` value used (e.g. `0.00/`, `0.05/`).

| Path | Description |
|---|---|
| `<comparison>_summary.tsv` | Transcript-level summary. One row per `transcript_id`: gene name, chromosome, strand, transcript length, number of edit sites, total G/read counts, mean editing percentage, and a comma-separated list of `chr:pos:frac` for every contributing site. |
| `<comparison>_summary_gene.tsv` | Gene-level summary. One row per `gene_name`, deduplicated at the genomic-site level (so a site counted under multiple transcripts of the same gene is not double-counted): number of transcripts, chromosome, strand, number of edit sites, total G/read counts, mean editing percentage, and a comma-separated list of `chr:pos` for every contributing site. 

---

## Notes

- `<treatment>_vs_<control>` and `<comparison>` refer to the `treatment`/`control` pairing defined in `--comparisons_csv` (e.g. `HyperTRIBE_vs_wtRNA`). If a `background` condition is set for that comparison, its own `<background>_vs_<control>` files are generated as an intermediate (visible in `filtered_replicate_bins/` and `combined_tsvs/`) and consumed by `SUBTRACTBKG` before annotation, no separate top-level output is published for the background comparison itself beyond those intermediates.
- Site coordinates throughout are 0-based half-open BED-style intervals (`start`, `end`), except within the `edit_sites` list columns in the summary tables, which report the 1-based `end` position only, for compactness.
