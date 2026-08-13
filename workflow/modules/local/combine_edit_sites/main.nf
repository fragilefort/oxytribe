process COMBINE_EDIT_SITES {
    tag "${condition}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:2a88e5229d866880a9c667af255a38fe385aaf9608c9b86d1891a038a5e5a479'
        : 'ghcr.io/fragilefort/oxytribe@sha256:2a88e5229d866880a9c667af255a38fe385aaf9608c9b86d1891a038a5e5a479'}"

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
