process COMBINE_CONTROLS {
    tag "${condition}"
    label 'process_single'

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'
        : 'ghcr.io/fragilefort/oxytribe@sha256:994fabbe29b2e59aab1708a457d70b52c5d03d1f41737f3c31df164674c8d974'}"

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
