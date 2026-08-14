# Parameters

## Input/output

### `--input_csv`
Path to the sample sheet csv file. Each row represents one sample.

| Column | Required | Description |
|---|---|---|
| `sample` | Yes | Unique sample identifier |
| `fastq_1` | Yes | Path to r1 fastq file (or single-end fastq) |
| `fastq_2` | No | Path to r2 fastq file. Leave empty for single-end |
| `condition` | Yes | Condition label used to group replicates |

**Example:**
```csv
sample,fastq_1,fastq_2,condition
rep1,/data/rep1_R1.fastq,/data/rep1_R2.fastq,treatment
rep2,/data/rep2_R1.fastq,/data/rep2_R2.fastq,treatment
ctrl,/data/ctrl_R1.fastq,/data/ctrl_R2.fastq,control
```
**Default:** `${baseDir}/contract/input.csv`

---

### `--comparisons_csv`
Path to the comparisons CSV file. Each row defines one treatment vs control comparison.

| Column | Required | Description |
|---|---|---|
| `treatment` | Yes | Condition name to use as the treatment group |
| `control` | Yes | Condition name to use as the control group |
| `background` | No | Condition whose editing sites will be subtracted from the results |

The `background` column is used to remove non-specific editing events, typically the HyperADARcd-alone condition. If omitted, no background subtraction is performed.

**Example:**
```csv
treatment,control,background
HyperTRIBE,wtRNA,HyperADARcd
```
**Default:** `${baseDir}/contract/comparisons.csv`

---

### `--chr_list`
Path to a plain text file listing chromosome names to include in the analysis, one per line. Reads mapping to chromosomes not in this list are silently discarded. Names must exactly match the chromosome naming convention of your reference genome (e.g. UCSC uses `chr1`, Ensembl uses `1`).

**Default:** `${baseDir}/contract/chromosomes.txt`

---

## Reference genome

### `--ref_fasta`
Path to the reference genome fasta file. Used for STAR genome generation and `samtools faidx`.

**Default:** `${baseDir}/contract/ref_genome/mm39.fa`

---

### `--ref_gtf`
Path to the genome annotation GTF file. Used for STAR genome generation and as the source for `PREPARE_GENE_SPANS`, which derives one synthetic gene-body interval (min start -> max end across all exons/transcripts) per `gene_id`, keyed by `(gene_id, chromosome)` so genes sharing a symbol across multiple loci each retain their own span rather than being collapsed into one. This synthetic-span GTF, not the raw `ref_gtf` is what `find_edit_sites` uses to resolve editing strand (A->G vs T->C) and restrict candidates to positions inside a gene, and what `BEDTOOLS_INTERSECT` uses for final gene annotation. 

**Default:** `${baseDir}/contract/ref_genome/mm39.gtf`

---

### `--star_ignore_sjdbgtf`
If `true`, STAR will not use the GTF file for splice junction annotation during alignment. Useful when aligning to a genome without a splice junction database or when using a custom index.

| Type | Default |
|---|---|
| `Boolean` | `false` |

---

## Alignment and preprocessing

### `--umi`
Set to `true` if your library was prepared with UMI barcodes. When enabled, UMI-based deduplication (`UMI-tools dedup`) is used instead of coordinate-based deduplication (`samtools markdup`).

| Type | Default |
|---|---|
| `Boolean` | `true` |

---

### `--skip_umiextract`
Set to `true` to skip the UMI extraction step (`UMI-tools extract`). Use this when UMI barcodes have already been moved into the read names upstream. Only relevant when `--umi true`.

| Type | Default |
|---|---|
| `Boolean` | `true` |

---

### `--skip_trimming`
Set to `true` to skip adapter trimming with Cutadapt. Use this when reads have already been trimmed upstream.

| Type | Default |
|---|---|
| `Boolean` | `true` |

---

## Replicate combination

### `--combination_mode`
The logic used to combine independently-filtered treatment editing evidence across biological replicates within a condition (via `combine_edit_sites`).

| Value | Description |
|---|---|
| `and` | A position must be present in all replicates to be retained. Counts are summed across replicates. More strict, reduces false positives. |
| `or` | A position present in any replicate is retained. Counts are summed across contributing replicates. More sensitive, reduces false negatives. |

The original HyperTRIBE pipeline's `bedtools intersect -f 0.9 -r` between two replicates is approximately equivalent to `and` mode.

| Type | Default |
|---|---|
| `String` | `"or"` |

---

### `--control_combination`
The logic used to pool raw per-position counts across control replicates (via `COMBINE_CONTROLS`) into a single control map, before any editing-site thresholds are applied. Unlike `combination_mode`, this operates on raw, unfiltered counts (`sam_to_map` output) rather than already-filtered treatment sites, pooling maximizes control coverage rather than restricting to only positions every replicate happened to cover.

| Type | Default |
|---|---|
| `String` | `"and"` |

---

## Editing site thresholds

These thresholds are applied by `find_edit_sites` when comparing the treatment condition against the control. A candidate position must fall inside a gene span (from `PREPARE_GENE_SPANS`) to be considered at all; the gene's strand at that position determines whether A->G or T->C mismatches are used as the editing signal.

### `--min_control_coverage`
Minimum number of reads required at a position in the control sample. Positions with fewer reads in the control are excluded to avoid calling editing sites where the reference base is uncertain.

| Type | Default |
|---|---|
| `Integer` | `10` |

---

### `--max_control_edit_frac`
Maximum allowed edit fraction at a position in the control sample (G-fraction for + strand genes, C-fraction for - strand genes). Positions exceeding this threshold are excluded as likely SNPs or endogenous editing events rather than RBP-specific editing.

| Type | Default |
|---|---|
| `BigDecimal` | `0.005` (0.5%) |

---

### `--min_control_non_g_frac`
Minimum fraction of unedited-base reads required at a position in the control sample (A-fraction for + strand genes, T-fraction for - strand genes). Ensures the control position is predominantly the unedited reference base.

| Type | Default |
|---|---|
| `BigDecimal` | `0.8` (80%) |

---

### `--min_rna_edit_count`
Minimum number of editing-signal reads (raw count) required at a position in the treatment sample for it to be considered a candidate site. Matches the legacy pipeline's site-calling criterion (`> 0` editing reads), filtering on editing strength happens later, at the `edit_threshold` summarization step.

| Type | Default |
|---|---|
| `Integer` | `1` |

---

## Downstream analysis

### `--edit_threshold`
Minimum editing fraction used by `summarize_edit_sites.r` when producing the gene-level and transcript-level summary tables. Sites below this threshold are excluded from the summary output, and the value is also used as the output subdirectory name (`temple/summarized_edit_sites/<edit_threshold>/`). 

| Type | Default |
|---|---|
| `BigDecimal` | `0.05` (5%) |

