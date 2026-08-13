#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'
        : 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'}"

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
