#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:7fd7c1daa02017efd5433cebba48f2ebd11f174b2e2407f0ce9acafcb7f38d1b'
        : 'ghcr.io/fragilefort/oxytribe@sha256:7fd7c1daa02017efd5433cebba48f2ebd11f174b2e2407f0ce9acafcb7f38d1b'}"

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
