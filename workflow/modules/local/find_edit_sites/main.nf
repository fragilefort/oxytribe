#!/usr/bin/env nextflow

process FIND_EDIT_SITES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(rna_bin), path(ctrl_bin)
    path chr_list
    val min_control_coverage
    val max_control_edit_frac
    val min_control_non_g_frac
    val min_rna_coverage
    val min_rna_edit_frac

    output:
    tuple val(meta), path("${meta.id}.tsv")

    script:
    """
    find_edit_sites \
        ${rna_bin} \
        ${ctrl_bin} \
        ${meta.id}.tsv \
        ${chr_list} \
        --min-control-coverage ${min_control_coverage} \
        --max-control-edit-frac ${max_control_edit_frac} \
        --min-control-non-g-frac ${min_control_non_g_frac} \
        --min-rna-coverage ${min_rna_coverage} \
        --min-rna-edit-frac ${min_rna_edit_frac}
    """
}
