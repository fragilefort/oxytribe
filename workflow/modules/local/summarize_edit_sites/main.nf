#!/usr/bin/env nextflow

process SUMMARIZE_EDIT_SITES {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'ghcr.io/fragilefort/oxytribe-r@sha256:4c51bdd260ab0d96702f91f050ccd8b576dc6f81c136c6f39a305a82f670eec6'
        : 'ghcr.io/fragilefort/oxytribe-r@sha256:4c51bdd260ab0d96702f91f050ccd8b576dc6f81c136c6f39a305a82f670eec6'}"

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
