#!/usr/bin/env nextflow

process COMBINE_EDIT_SITES {
    tag "${condition}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'
        : 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'}"

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
