# Miscellaneous configuration notes

This document covers pipeline configuration beyond the core `params` (see `params.md`): the Nextflow `process` scope and profile-based config layering, plus known gotchas worth knowing before running on new data.

## `process` config (`process.config` or the `process {}` block)

Nextflow's `process` scope sets resource allocation and tool-specific arguments per process, separately from the biological/pipeline `params`. Two mechanisms are used throughout this pipeline:

### Resource allocation via labels

Processes are tagged with a `label` (`process_single`, `process_low`, `process_medium`, `process_high`) in their module definition, and `process.config` assigns CPU/memory per label:

```groovy
process {
    cpus = 4
    memory = '16 GB'

    withLabel: 'process_low' {
        cpus = 4
        memory = '16 GB'
    }
    withLabel: 'process_medium' {
        cpus = 16
        memory = '64 GB'
    }
    withLabel: 'process_high' {
        cpus = 32
        memory = '128 GB'
    }
}
```

The top-level `cpus`/`memory` (outside any `withLabel`) are the default for any process without a matching label. Adjust these based on your cluster/machine's actual resources, `STAR_GENOMEGENERATE` and `STAR_ALIGN` in particular are memory-requiring on large genomes and are given `process_high`-equivalent resources explicitly rather than relying on label defaults.

### Tool arguments via `withName`

Per-process command-line arguments (adapter sequences, quality trimming thresholds, alignment filters, output formatting) are set via `ext.args`, keyed by process name:

```groovy
withName: 'CUTADAPT' {
    ext.args = [
        '--trim-n',
        '-q 25,25',
        ...
    ].join(' ')
}
```

`ext.prefix` and `ext.suffix` control output file naming per process where relevant (e.g. `SAMTOOLS_VIEW`'s `ext.prefix = { "${meta.id}_view" }`).

**Important:** these values are dataset- and protocol-specific. Adapter sequences, UMI barcode patterns (`UMITOOLS_EXTRACT`'s `--bc-pattern`), quality thresholds, and STAR filter stringency (`--outFilterMismatchNoverLmax`, `--outFilterMultimapNmax`, etc.) used for the bundled test profile are tuned for that specific published dataset and will not be correct defaults for arbitrary new data. Review and override every `withName` block relevant to your library prep and organism before running on a new dataset, don't assume the test profile's tool arguments generalize.

**Note on STAR_ALIGN:** The first rust crate `SAM_TO_MAP` depends entirly on the existance of the MD tag in the aligned reads, so it cannot be excluded when running alignment.

## Profiles (`profiles.config` or `-profile` on the CLI)

Profiles layer additional config on top of the base `nextflow.config`/`process.config`, selected via `-profile <name>[,<name>...]` on the command line. This pipeline uses profiles for two independent concerns that get combined:

- **Container engine**: `docker`, `singularity`, selects how process containers are run. Exactly one of these should be specified.
- **Dataset/test config**: `test`, points `params` at the bundled small validation dataset (`workflow/contract/test/`), used for testing and the validation work documented in `validation.md`.

Example: `-profile test,singularity` runs the bundled test dataset under Singularity. For a real dataset, you supply your own `params.config` (or CLI `--param value` overrides) instead of `test`, combined with whichever container engine profile matches your environment.

Profiles are additive don't assume selecting `test` alone is sufficient without also specifying a container engine, and don't assume a profile you didn't specify is implicitly applied.

---

## Known gotchas

### `--skip_trimming true` requires pre-gzipped input fastqs

`STAR_ALIGN` is configured in test profile with:

```groovy
withName: 'STAR_ALIGN' {
    ext.args = [
        ...
        '--readFilesCommand zcat',
    ].join(' ')
}
```

`--readFilesCommand zcat` unconditionally decompresses every input fastq through `zcat` before STAR reads it. This is correct when reads have passed through `CUTADAPT` first (its output is gzip-compressed), but **if `params.skip_trimming = true`, raw input fastqs go directly to STAR**, and if those raw files are plain, uncompressed text (not `.gz`), `zcat` fails on them.

**This failure is silent.** `zcat` on a non-gzip file produces empty output, not an error; STAR receives an empty read stream, reports `Number of input reads: 0` in its own log, and exits with status 0, no crash, no error message anywhere in the pipeline. Every downstream step (dedup, `SAMTOOLS_VIEW`, `sam_to_map`, `find_edit_sites`, etc.) then correctly processes an empty file, silently propagating the emptiness all the way to zero-byte final output with no indication of where the actual problem originated. Tracing this back to its source without knowing to check `STAR_ALIGN`'s own log first can cost significant debugging time.

**Before running with `--skip_trimming true`:** confirm your input fastqs are actually gzip-compressed (`file <path>` should report `gzip compressed data`, not `ASCII text`). If they're plain text, either:
- gzip them and point `input_csv` at the `.gz` versions, or
- remove `'--readFilesCommand zcat'` from `STAR_ALIGN`'s `ext.args` in your run's `process.config` if your fastqs will always be uncompressed for that run.

**How to check if this has already happened on a run:** for any sample producing empty/zero-byte output downstream, check that sample's own STAR log directly before investigating anything else:

```bash
find <work dir> -iname "*<sample_id>*Log.final.out" | xargs grep "Number of input reads"
```

`Number of input reads: 0` confirms this issue at the earliest possible stage, saving the need to trace through every downstream process individually.
