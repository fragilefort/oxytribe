#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params {
    input_csv: Path
    ref_fasta: Path
    ref_gtf: Path
    star_ignore_sjdbgtf: Boolean
    chr_list: Path
    combination_mode: String
    comparisons_csv: Path
    umi: Boolean
    outdir: Path
    min_control_coverage: Integer
    max_control_edit_frac: Double
    min_control_non_g_frac: Double
    min_rna_coverage: Integer
    min_rna_edit_frac: Double
}

include { FASTQC } from './modules/nf-core/fastqc/main.nf'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main.nf'
include { STAR_ALIGN } from './modules/nf-core/star/align/main.nf'
include { SAM_TO_MAP } from './modules/local/sam_to_map.nf'
include { COMBINE_EDIT_SITES } from './modules/local/combine_edit_sites.nf'

workflow {

    main:
    reads_ch = channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, single_end: row.fastq_2 == '', condition: row.condition]
            def files = meta.single_end
                ? [file(row.fastq_1)]
                : [file(row.fastq_1), file(row.fastq_2)]
            [meta, files]
        }
    FASTQC(reads_ch)

    reffasta_ch = channel.fromPath(params.ref_fasta)
        .map { fasta -> [[id: 'genome'], fasta] }
    refgtf_ch = channel.fromPath(params.ref_gtf)
        .map { gtf -> [[id: 'genome'], gtf] }
    STAR_GENOMEGENERATE(reffasta_ch, refgtf_ch)

    STAR_ALIGN(
        reads_ch,
        STAR_GENOMEGENERATE.out.index.first(),
        refgtf_ch.first(),
        params.star_ignore_sjdbgtf ?: false,
    )
    star_all_logs = STAR_ALIGN.out.log_final.mix(
        STAR_ALIGN.out.log_out,
        STAR_ALIGN.out.log_progress,
    )

    SAM_TO_MAP(
        STAR_ALIGN.out.sam,
        params.chr_list,
    )

    grouped_bins = SAM_TO_MAP.out.map
        .map { meta, bin -> [meta.condition, bin] }
        .groupTuple()
    COMBINE_EDIT_SITES(
        grouped_bins,
        params.combination_mode,
    )

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    star_index = STAR_GENOMEGENERATE.out.index
    star_align_log = star_all_logs
    sam_files = STAR_ALIGN.out.sam
    sam_maps = SAM_TO_MAP.out.map
    combined_maps = COMBINE_EDIT_SITES.out
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
    star_align_log {
        path "star/align/log"
    }
    sam_files {
        path "star/align/sam"
    }
    sam_maps {
        path "temple/map"
    }
    combined_maps {
        path "temple/combined_maps/${params.combination_mode}"
    }
}
