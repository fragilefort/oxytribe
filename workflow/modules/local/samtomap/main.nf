#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:4a90665e21e712f3d82cfcf36a0a120bc400dfe12da5a2afc58fbc64c60c2b71'
        : 'ghcr.io/fragilefort/oxytribe@sha256:4a90665e21e712f3d82cfcf36a0a120bc400dfe12da5a2afc58fbc64c60c2b71'}"

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
