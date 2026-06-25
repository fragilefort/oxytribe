#!/usr/bin/env nextflow

process COMBINE_EDIT_SITES {
    tag "${condition}"

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
