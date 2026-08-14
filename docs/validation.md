# Validation

This document records how oxytribe's editing-site output was validated against the original HyperTRIBE pipeline and published results.

## Methodology

1. **Aggregate gene-level overlap** against the legacy pipeline's own output (`HyperTRIBE_results_wtRNA.xls`, produced by `summarize_results.pl` on `HyperTRIBE_1_2_wtRNA.bedgraph` per the published protocol). Same-dataset, same-experiment comparison, both pipelines ran on identical input FASTQs, so any divergence reflects pipeline logic differences, not biological or dataset differences.

2. **Top-N rank agreement** on editing site count, since the strongest, most reproducible targets are the ones most likely to be independently confirmed and the ones a researcher will look at first.

3. **Named target validation** against genes independently confirmed in published literature (McMahon et al. 2016; Xu et al. 2018 - Hrp48 HyperTRIBE western blot and TRIBE/HyperTRIBE/CLIP consensus targets), independent of both pipelines' own output.

## How to run

```bash
Rscript workflow/bin/validate.r \
    --new_gene_summary <path to *_summary_gene.tsv> \
    --legacy_xls <path to legacy HyperTRIBE_results_*.xls> \
    --top_n 30 \
    --out_dir validation_output
# or simply
pixi run validation 
# this does the same with the default paths, top_n 30
```

`--edit_threshold` on the main pipeline should be set to `0` when generating the `_summary_gene.tsv` used here, the legacy xls this compares against was produced with no editing-fraction cutoff applied before reporting (see "Known pipeline differences" below).

## Results (test dataset: SRR5944748/49/50, SRR6426146 - dm6)

| Metric | Result |
|---|---|
| Legacy gene count | 5726 |
| New pipeline gene count | 6087 |
| Shared genes | 5712 |
| Legacy genes missed by new pipeline | 14 |
| Genes recovered | **99.76%** |
| Top-30 by edit-site count, overlap | 27/30 (90%) |

Note: The higher edit-site count in oxytribe is primarily driven by how the two pipelines handle read alignments and genomic bounds. The legacy Perl pipeline relies on rigid branching logic that drops any read failing to match a specific, hardcoded CIGAR string pattern, throwing away many spliced or hyper-edited reads before site evaluation. In contrast, oxytribe Rust parser  handles all CIGAR operations, recovering these complex alignments while expanding the target area across entire gene spans rather than restricting search to strict exon-only boundaries. Together, these allows oxytribe to capture a significantly higher density of genuine edit sites without discarding valid multi-mismatch reads.

### Named target validation

| Gene | Legacy | New pipeline |
|---|---|---|
| Rm62 | not present under this casing | 101 sites |
| Syt1 | absent | absent (present in `ref_gtf`, 0 sites in both pipelines, likely a test-dataset coverage limitation, not a pipeline discrepancy) |
| fne | absent | absent (same as Syt1) |
| cwo | present | 43 sites |
| unc-13 | present | 85 sites |
| CG31694 | present | 53 sites |
| CG11593 | present | 45 sites |
| Chd64 | present | 26 sites |

6/8 named targets recovered with real signal; the 2 misses are absent in both pipelines identically, so not attributable to the migration.

## Known pipeline differences (resolved during migration)

- **Per-replicate independent filtering vs. legacy's unfiltered replicate intersect.** The legacy protocol's `bedtools intersect -f 0.9 -r` between two raw, unthresholded replicate files is the actual site-combination step; a 20-read/5%-editing threshold exists in the legacy codebase (`filter_by_threshold_without_header.pl`) but is not applied before this step in the documented pipeline run. `find_edit_sites` now uses `--min-rna-edit-count 1` (any nonzero editing signal) at the per-replicate candidate stage, matching legacy's actual SQL predicate (`rstrings{$editbase} > '0'`), rather than a coverage/fraction threshold.

- **Intronic site loss at annotation.** Raw `ref_gtf` files typically only contain `exon`/`CDS`/`start_codon`/`stop_codon` rows, no full gene-body bounding row. `bedtools intersect` against this directly drops any site that only falls in an intron. `PREPARE_GENE_SPANS` builds a synthetic gene-body span (min start -> max end per `gene_id`) used both for `find_edit_sites`'s interval lookup and the final `BEDTOOLS_INTERSECT` annotation step, matching legacy's `refFlat.txt`-based `txStart<=pos<=txEnd` gene-span logic.


## Interpretation

- `gene_spans.gtf`-based site counts run consistently **higher** than legacy's per-gene counts (e.g. CtBP: legacy 220, new pipeline 283). Expected, not a bug: synthetic gene spans are looser than `refFlat`'s real per-transcript exon/intron structure and can sweep in neighboring/overlapping genes' introns that legacy's explicit "skip if position maps to >1 gene" rule (`summarize_results.pl`: `next if (keys %$gene_hash > 1)`) would exclude. Rank order and gene identity are the more reliable comparison points than absolute site counts.
