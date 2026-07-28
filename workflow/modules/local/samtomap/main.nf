#!/usr/bin/env nextflow

process SAM_TO_MAP {
    tag "${meta.id}"

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
