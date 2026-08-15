# OXYTRIBE
OXYTRIBE is a Nextflow pipeline for HyperTRIBE RNA editing identification and gene-level target prioritization. It is a full rewrite of the original [HyperTRIBE](https://github.com/rosbashlab/HyperTRIBE) computational pipeline, replacing the Perl/Python/MySQL stack with Rust binaries and Nextflow, removing the database dependency entirely, and adding support for multi-replicate comparison with configurable AND/OR logic.

Note: this is not a typical rewrite where the same logic is retained. OXYTRIBE uses a substantially different approach across all main processing steps.

## Documentation

- [Parameters](docs/parameters.md) : configuration and command-line parameters
- [Output](docs/output.md) : output files and their structure
- [Validation](docs/validation.md) : validation and comparisons with the HyperTRIBE pipeline

## TODO

- [x] Converting `sam_to_matrix.pl` to Rust
- [x] Using a simpler approach without Database
- [x] All the processes until the final TSV file are either in Nextflow or Rust
- [x] Wrap the workflow with Nextflow
- [x] Delete the modules that already exist in nf-core
- [x] Annotate the editing sites in the TSV file
- [x] Add an R script for top edited genes with threshold like 5% and 1% and some stats
- [x] Add UMI extract
- [x] USE UMI dedup
- [x] OW, use cordinate based dedupping
- [x] Trimming process
- [x] Add containers for R analysis script
- [x] Allow for subtracting the editing sites from controlVscontrol
- [x] Add containers for the rust binaries to support multiple platforms
- [x] Testing
- [x] Validation using the HyperTRIBE sample data
- [x] Update the documentation for installation and test profile
- [x] Docs for parameters
- [x] Docs for validation
- [x] Docs for outputs

## Architecture

The pipeline is implemented using Nextflow (DSL2). Read processing, alignment, and QC are handled using standard bioinformatics modules from nf-core, followed by three custom Rust binaries for candidate site identification. Downstream processing uses R and Bedtools modules for background subtraction, GTF annotation, and target summarization.

All modules support Docker containers and Apptainer/Singularity. Custom binaries use SHA digests to ensure reproducibility.

![](./assets/arc.png)

## Installation

OXYTRIBE uses [Pixi](https://pixi.sh) to manage dependencies including Nextflow, nf-core tools, R, and the SRA toolkit.

```bash
# Install Pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | sh

# Clone the repository
git clone https://github.com/fragilefort/oxytribe.git
cd oxytribe

# Install dependencies
pixi install
```

Containers (Docker or Singularity/Apptainer) are required to run the pipeline because nf-core modules and the custom Rust binaries are distributed as container images. A compatible container runtime must be available on your system.

## Usage

OXYTRIBE is run through Nextflow. The main workflow is:

```bash
pixi run nextflow run workflow/main.nf
```

The workflow uses CSV files to define the input samples and treatment/control comparisons.

### Basic usage

Edit the paths, file names in `input.csv`, tweak the parameters and fire it up:

```bash
pixi run nextflow run workflow/main.nf
```

By default, OXYTRIBE expects:

- `${baseDir}/contract/input.csv` : sample sheet
- `${baseDir}/contract/comparisons.csv` : treatment/control comparisons
- `${baseDir}/contract/chromosomes.txt` : chromosomes to include
- `${baseDir}/contract/ref_genome/mm39.fa` : reference genome
- `${baseDir}/contract/ref_genome/mm39.gtf` : genome annotation

For a custom analysis, these can be overridden with the corresponding parameters, or simply edit the config files:

```bash
pixi run nextflow run workflow/main.nf \
    --input_csv /path/to/input.csv \
    --comparisons_csv /path/to/comparisons.csv \
    --ref_fasta /path/to/genome.fa \
    --ref_gtf /path/to/annotation.gtf \
    --chr_list /path/to/chromosomes.txt
```

### Replicate combination

OXYTRIBE supports configurable replicate-combination logic.

To require an editing site to be present in all treatment replicates:

```bash
pixi run nextflow run workflow/main.nf \
    --input_csv /path/to/input.csv \
    --comparisons_csv /path/to/comparisons.csv \
    --combination_mode and
```

The default is `or`, where a position present in any replicate is retained.

Control replicates can be combined independently using `--control_combination`.

### UMI processing

UMI-based deduplication is enabled by default:

```bash
pixi run nextflow run workflow/main.nf \
    --input_csv /path/to/input.csv \
    --comparisons_csv /path/to/comparisons.csv \
    --umi true
```

If UMIs have already been extracted into the read names upstream, UMI extraction can be skipped:

```bash
pixi run nextflow run workflow/main.nf \
    --umi true \
    --skip_umiextract true
```

Adapter trimming can similarly be skipped with:

```bash
pixi run nextflow run workflow/main.nf \
    --skip_trimming true
```

### Editing-site thresholds

Editing-site thresholds can be adjusted using the corresponding parameters. For example:

```bash
pixi run nextflow run workflow/main.nf \
    --min_control_coverage 10 \
    --max_control_edit_frac 0.005 \
    --min_control_non_g_frac 0.8 \
    --min_rna_edit_count 1 \
    --edit_threshold 0.05
```

See [Parameters](docs/parameters.md) for the complete list of parameters, defaults, input formats, and descriptions.

## Quick Start

The repository includes a test workflow using the validation dataset from the original HyperTRIBE paper (test dataset: SRR5944748/49/50, SRR6426146 - dm6).

Run the test using Docker:

```bash
pixi run test-docker
```

Or using Singularity/Apptainer:

```bash
pixi run test-singularity
```

The test workflow:

1. Downloads four FASTQ files from SRA:
   - Hrp48 HyperTRIBE replicate 1
   - Hrp48 HyperTRIBE replicate 2
   - wtRNA control
   - HyperADARcd background
2. Downloads the dm6 reference genome from UCSC.
3. Runs the full OXYTRIBE pipeline end-to-end using the test profile.

Expected output is written to `ChinaTown/`, including editing sites for the `HyperTRIBE_vs_wtRNA` comparison.

The above test designed for server runs, to run it locally, we put a cap on available resources:
```bash
pixi run download-test-data
pixi run nextflow run workflow/main.nf -profile local,test,docker

```

See [Validation](docs/validation.md) for details about the test dataset, validation procedure, and comparison with the original HyperTRIBE pipeline.

## Output

The main output files contain editing sites and gene/transcript-level summaries generated from the treatment/control comparisons.

See [Output](docs/output.md) for a description of the output directory structure and individual output files.

## Contributions

Contributions are welcome via issues and pull requests.

# HyperTRIBE
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. This is an improved version of a previously developed technique called TRIBE (Targets of RNA-binding proteins Identified By Editing).

Please get the code from GitHub or download from here: [![DOI](https://zenodo.org/badge/114820120.svg)](https://zenodo.org/badge/latestdoi/114820120)

Detailed doumentation is available at: http://hypertribe.readthedocs.io/en/latest/

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.
