process FIND_EDIT_SITES {
    tag "${meta.id}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'
        : 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'}"

    input:
    tuple val(meta), path(rna_bin), path(ctrl_bin)
    path gene_spans
    path chr_list
    val min_control_coverage
    val max_control_edit_frac
    val min_control_non_g_frac
    val min_rna_edit_count

    output:
    tuple val(meta), path("${meta.id}.bin"), emit: bin

    script:
    """
    find_edit_sites \\
        ${rna_bin} \\
        ${ctrl_bin} \\
        ${meta.id}.bin \\
        ${gene_spans} \\
        ${chr_list} \\
        --min-control-coverage ${min_control_coverage} \\
        --max-control-edit-frac ${max_control_edit_frac} \\
        --min-control-non-g-frac ${min_control_non_g_frac} \\
        --min-rna-edit-count ${min_rna_edit_count}
    """
}
