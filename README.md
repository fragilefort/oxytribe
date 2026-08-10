# OXYTRIBE
OXYTRIBE is a Nextflow pipeline for HyperTRIBE RNA editing identification and gene-level target prioritization. It is a full rewrite of the original [HyperTRIBE](https://github.com/rosbashlab/HyperTRIBE) computational pipeline, replacing the Perl/Python/MySQL stack with Rust binaries and Nextflow, removing the database dependency entirely, and adding support for multi-replicate comparison with configurable AND/OR logic. Note, this is not a typical rewrite where the same logic is retained, OXYTRIBE uses a very different approach across all main steps from HyperTRIBE.

## TODO

- [x] Converting sam_to_matrix.pl to Rust
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
- [ ] Update the documentation for installation and usage


## Architecture
he pipeline is handled  using Nextflow (DSL2). Read processing, alignment, and QC are handled via standard bioinformatics modules (nf-core), followed by a custom 3 Rust binaries for candidate site identification, and downstream R / Bedtools modules for background subtraction, GTF annotation, and target summarization.All modules here support docker containers and apptainers, custom binaries uses sha digests to insure reproducablity.

![](./assets/arc.png)
 
## Installation
To install OXYTRIBE, you only need [Pixi](https://pixi.sh) to manage reproducible Conda and R environments. It also have all the software to run the pipeline like nexflow, nf-core. It was tested with `pixi 0.72.2` but should be fine with other versions.

```bash
git clone https://github.com/fragilefort/oxytribe.git
cd oxytibe
pixi install
```

## Quick start
You can run the test profile using:
```bash
pixi run test
```
This does the following:
1. Download test data from NCBI, 4 fastq files used in the HyperTRIBE paper along with downloading the fasta ref genome. 
2. Run end to end pipelines on these samples using the parameters in `workflow/conf/test.config`. This also uses addtional files which are locates in `workflow/contract/test`, these files include chromosme names, input sample sheet and comparsions sample sheet.


# HyperTRIBE
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. This is an improved version of a previously developed technique called TRIBE (Targets of RNA-binding proteins Identified By Editing).

Please get the code from GitHub or download from here: [![DOI](https://zenodo.org/badge/114820120.svg)](https://zenodo.org/badge/latestdoi/114820120)

Detailed doumentation is available at: http://hypertribe.readthedocs.io/en/latest/

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.
