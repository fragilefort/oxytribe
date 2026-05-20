process STAR_INDEX {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(fasta)
    path gtf

    output:
    tuple val(meta), path("star"), emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir star

    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir star \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        ${args}

    cat <<EOF > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    EOF
    """

    stub:
    """
    mkdir star
    touch star/Genome
    touch star/SA
    touch versions.yml
    """
}
