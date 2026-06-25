#!/usr/bin/env nextflow

process FIND_EDIT_SITES {
    input:
    tuple val(meta), path(rna_bin), path(ctrl_bin)
    path chr_list

    output:
    tuple val(meta), path("${meta.id}.tsv")

    script:
    """
    find_edit_sites \
        ${rna_bin} \
        ${ctrl_bin} \
        ${meta.id}.tsv \
        ${chr_list} \
        --min-control-coverage ${params.min_control_coverage} \
        --max-control-edit-frac ${params.max_control_edit_frac} \
        --min-control-non-g-frac ${params.min_control_non_g_frac} \
        --min-rna-coverage ${params.min_rna_coverage} \
        --min-rna-edit-frac ${params.min_rna_edit_frac}
    """
}
