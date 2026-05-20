process STAR_INDEX {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(genome)

    output:
    path ("*.tab"), emit: tab
    path ("*.txt"), emit: txt
    path ("SA*"), emit: SA
    path ("Genome"), emit: genome
    path ("Log.out"), emit: log
    tuple val("${task.process}"), val('STAR'), eval('STAR --version '), emit: versions_star, topic: versions

    script:
    """
    // The number of threads should not be defined here
    STAR --runThreadN 10 \
        --runMode genomeGenerate \
        // How to define the output path?
        --genomeDir results/star/index \
        // Again hardcoded
        --genomeFastaFiles 03_data/genome_files/assembly.fasta \
        --sjdbGTFfile data/genome_files/assembly.gtf \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 11\
        // We also should allow user for extra config
    """

    stub:
    """
    STAR --version
    
    """
}
