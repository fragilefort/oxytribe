#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:a6741c336433d07ea4475f56dd18a3ef5d4d0100ee54efe639016618c8953ef9'
        : 'ghcr.io/fragilefort/oxytribe@sha256:a6741c336433d07ea4475f56dd18a3ef5d4d0100ee54efe639016618c8953ef9'}"

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
