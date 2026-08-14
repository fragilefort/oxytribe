process COMBINE_CONTROLS {
    tag "${condition}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:6cace988efe7bd1fd4961b8997f6c12a6ade940264b5f77e55d46df7b30a166d'
        : 'ghcr.io/fragilefort/oxytribe@sha256:6cace988efe7bd1fd4961b8997f6c12a6ade940264b5f77e55d46df7b30a166d'}"

    input:
    tuple val(condition), path(bins)
    val mode

    output:
    tuple val(condition), path("${condition}.bin"), emit: bin

    script:
    """
    combine_edit_sites ${bins} \\
        --output ${condition}.bin \\
        --mode ${mode} \\
        --raw
    """
}
