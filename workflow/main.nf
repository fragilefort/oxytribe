#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params {
    input_path: Path
}

include { FASTQC } from './modules/nf-core/fastqc/main.nf'

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

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
}

output {
    fastqc_html {
        path "fastqc/html"
    }
    fastqc_zip {
        path "fastqc/zip"
    }
}
