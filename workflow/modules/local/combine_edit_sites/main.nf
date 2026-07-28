#!/usr/bin/env nextflow

process COMBINE_EDIT_SITES {
    tag "${condition}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:2dc8b045a9a9fecd186209712ec083a94bda15ce8305f438bccdd1ed6feeea12'
        : 'ghcr.io/fragilefort/oxytribe@sha256:2dc8b045a9a9fecd186209712ec083a94bda15ce8305f438bccdd1ed6feeea12'}"

    input:
    tuple val(condition), path(bins)
    val mode

    output:
    tuple val(condition), path("${condition}.bin")

    script:
    """
    combine_edit_sites ${bins} \
        --output ${condition}.bin \
        --mode ${mode}
    """
}
