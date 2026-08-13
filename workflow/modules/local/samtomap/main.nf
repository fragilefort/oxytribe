#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'
        : 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'}"

    input:
    tuple val(meta), path(sam)
    path chr_list

    output:
    tuple val(meta), path("${meta.id}_map.bin"), emit: map

    script:

    """
    sam_to_map ${sam} ${meta.id}_map.bin ${chr_list}
    """
}
