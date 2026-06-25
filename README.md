# oxytribe
This is an attempt to improve the HyperTRIBE software by moving some parts from Perl to Rust and Python.
## TODO
[x] Converting sam_to_matrix.pl to Rust
[x] Using a simpler approach without Database
[x] All the processes until the final tsv file is either in nf or in Rust
[x] Wrap the workflow with a nextflow
[x] Delete the modules that already exists in nf-core
[] Annotate the editing sites in the tsv file
[] Containerize the pixi enviroment
[] Update the documentation for installation and usage
[] Testing and benchmarking

# HyperTRIBE
HyperTRIBE is a technique used for the identification of the targets of RNA binding proteins (RBP) in vivo. This is an improved version of a previously developed technique called TRIBE (Targets of RNA-binding proteins Identified By Editing).

Please get the code from GitHub or download from here: [![DOI](https://zenodo.org/badge/114820120.svg)](https://zenodo.org/badge/latestdoi/114820120)

Detailed doumentation is available at: http://hypertribe.readthedocs.io/en/latest/

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.
