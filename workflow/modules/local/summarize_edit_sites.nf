#!/usr/bin/env nextflow

process SUMMARIZE_EDIT_SITES {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'oras://ghcr.io/fragilefort/oxytribe-r@sha256:7fd2a3a37fc9401e89748115e9a9ed8baa093138445506802ce69d7d73c8ff2a'
        : 'ghcr.io/fragilefort/oxytribe-r@sha256:bc6d4009679fb3a3615b4fb095e272f7c5697e4ec0b15487de83c2430270b087'}"

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
