#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params {
    input_csv: Path
    ref_fasta: Path
    ref_gtf: Path
}

include { FASTQC } from './modules/nf-core/fastqc/main.nf'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main.nf'

workflow {

    main:
    reads_ch = channel.fromPath(params.input_path)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, single_end: row.fastq2 == '']
            def files = meta.single_end
                ? [file(row.fastq1)]
                : [file(row.fastq_1), file(row.fastq_2)]
            [meta, files]
        }
    FASTQC(reads_ch)

    reffasta_ch = channel.fromPath(params.ref_fasta)
        .map { fasta -> [[id: 'genome'], fasta] }
    refgtf_ch = channel.fromPath(params.ref_gtf)
        .map { gtf -> [[id: 'genome'], gtf] }
    STAR_GENOMEGENERATE(reffasta_ch, refgtf_ch)

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    star_index = STAR_GENOMEGENERATE.out.index
}

output {
    fastqc_html {
        path "fastqc/html"
    }
    fastqc_zip {
        path "fastqc/zip"
    }
    star_index {
        path "star/index"
    }
}
