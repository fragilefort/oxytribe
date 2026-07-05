#!/usr/bin/env nextflow

process SUMMARIZE_EDIT_SITES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(tsv)
    val threshold

    output:
    tuple val(meta), path("${meta.id}_summary.tsv"), path("${meta.id}_summary_gene.tsv")

    script:
    """
    summarize_edit_sites.r \
        --input ${tsv} \
        --output ${meta.id}_summary.tsv \
        --threshold ${threshold}
    """
}
