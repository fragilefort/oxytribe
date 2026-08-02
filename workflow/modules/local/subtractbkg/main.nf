process SUBTRACTBKG {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'ghcr.io/fragilefort/oxytribe-r@sha256:4c51bdd260ab0d96702f91f050ccd8b576dc6f81c136c6f39a305a82f670eec6'
        : 'ghcr.io/fragilefort/oxytribe-r@sha256:4c51bdd260ab0d96702f91f050ccd8b576dc6f81c136c6f39a305a82f670eec6'}"

    input:
    tuple val(meta), path(target_tsv), path(bg_tsv)

    output:
    tuple val(meta), path("${meta.id}_subtracted.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bg_arg = bg_tsv instanceof List && bg_tsv.isEmpty() ? "" : "--bg ${bg_tsv}"
    """
    subtract_background.r \
        --input ${target_tsv} \
        --output ${prefix}_subtracted.tsv \
        ${bg_arg}
    """
}
