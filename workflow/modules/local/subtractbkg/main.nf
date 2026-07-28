process SUBTRACTBKG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'oras://ghcr.io/fragilefort/oxytribe-r@sha256:7fd2a3a37fc9401e89748115e9a9ed8baa093138445506802ce69d7d73c8ff2a'
        : 'ghcr.io/fragilefort/oxytribe-r@sha256:bc6d4009679fb3a3615b4fb095e272f7c5697e4ec0b15487de83c2430270b087'}"

    input:
    tuple val(meta), path(target_tsv), path(bg_tsv)

    output:
    tuple val("${task.process}"), val('subtractbkg'), eval("subtractbkg --version"), topic: versions, emit: versions_subtractbkg

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bg_arg = bg_tsv ? "--bg ${bg_tsv}" : ""
    """
    subtract_background.r \
        --input ${target_tsv} \
        --output ${prefix}_subtracted.tsv \
        ${bg_arg}
    """
}
