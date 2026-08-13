process FIND_EDIT_SITES {
    tag "${meta.id}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'
        : 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'}"

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
