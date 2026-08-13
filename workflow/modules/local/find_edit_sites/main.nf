process FIND_EDIT_SITES {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:5ff5cc23a935ce9ad2ab90836fa4ef2f4ae1e4b0a46e00d82e7e1e3ceaf101fc'
        : 'ghcr.io/fragilefort/oxytribe@sha256:5ff5cc23a935ce9ad2ab90836fa4ef2f4ae1e4b0a46e00d82e7e1e3ceaf101fc'}"

    input:
    tuple val(meta), path(rna_bin), path(ctrl_bin)
    val min_control_coverage
    val max_control_edit_frac
    val min_control_non_g_frac
    val min_rna_coverage
    val min_rna_edit_frac

    output:
    tuple val(meta), path("${meta.id}.bin"), emit: bin

    script:
    """
    find_edit_sites \\
        ${rna_bin} \\
        ${ctrl_bin} \\
        ${meta.id}.bin \\
        --min-control-coverage ${min_control_coverage} \\
        --max-control-edit-frac ${max_control_edit_frac} \\
        --min-control-non-g-frac ${min_control_non_g_frac} \\
        --min-rna-coverage ${min_rna_coverage} \\
        --min-rna-edit-frac ${min_rna_edit_frac}
    """
}
