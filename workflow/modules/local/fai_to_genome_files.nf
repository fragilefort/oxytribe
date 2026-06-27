#!/usr/bin/env nextflow

process FAI_TO_GENOME_FILES {
    input:
    tuple val(meta), path(fai)

    output:
    path "chromosomes.txt", emit: chr_list
    path "chrom.sizes", emit: chrom_sizes

    script:
    """
    cut -f1 ${fai} > chromosomes.txt
    cut -f1,2 ${fai} > chrom.sizes
    """
}
