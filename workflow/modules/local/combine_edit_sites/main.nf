process COMBINE_EDIT_SITES {
    tag "${condition}"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'
        : 'ghcr.io/fragilefort/oxytribe@sha256:edfe94b517ab63e69ebf032188616b19de594e2269f9f81b24f27573194774ac'}"

    input:
    tuple val(condition), path(bins)
    val mode
    path chr_list

    output:
    tuple val(condition), path("${condition}.tsv"), emit: tsv

    script:
    """
    combine_edit_sites ${bins} \\
        --output ${condition}.tsv \\
        --mode ${mode} \\
        --chr-list ${chr_list}
    """
}
