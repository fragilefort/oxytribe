#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:3b8c8f06cccacc07705b3be579e796df166abf6eb4397cd9bb61722745cccb27'
        : 'ghcr.io/fragilefort/oxytribe@sha256:3b8c8f06cccacc07705b3be579e796df166abf6eb4397cd9bb61722745cccb27'}"

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
