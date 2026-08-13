process COMBINE_EDIT_SITES {
    tag "${condition}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:5ff5cc23a935ce9ad2ab90836fa4ef2f4ae1e4b0a46e00d82e7e1e3ceaf101fc'
        : 'ghcr.io/fragilefort/oxytribe@sha256:5ff5cc23a935ce9ad2ab90836fa4ef2f4ae1e4b0a46e00d82e7e1e3ceaf101fc'}"

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
