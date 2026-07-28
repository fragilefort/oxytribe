#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:2dc8b045a9a9fecd186209712ec083a94bda15ce8305f438bccdd1ed6feeea12'
        : 'ghcr.io/fragilefort/oxytribe@sha256:2dc8b045a9a9fecd186209712ec083a94bda15ce8305f438bccdd1ed6feeea12'}"

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
