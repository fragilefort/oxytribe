#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/fastqc/main.nf'
include { MULTIQC } from './modules/nf-core/multiqc/main.nf'
include { UMITOOLS_EXTRACT } from './modules/nf-core/umitools/extract/main.nf'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main.nf'
include { STAR_ALIGN } from './modules/nf-core/star/align/main.nf'
include { SAM_TO_MAP } from './modules/local/samtomap/main.nf'
include { COMBINE_CONTROLS } from './modules/local/combinecontrols/main.nf'
include { COMBINE_EDIT_SITES } from './modules/local/combine_edit_sites/main.nf'
include { FIND_EDIT_SITES } from './modules/local/find_edit_sites/main.nf'
include { PREPARE_GENE_SPANS } from './modules/local/makespan/main.nf'
include { BEDTOOLS_INTERSECT } from './modules/nf-core/bedtools/intersect/main.nf'
include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main.nf'
include { SUMMARIZE_EDIT_SITES } from './modules/local/summarize_edit_sites/main.nf'
include { CUTADAPT } from './modules/nf-core/cutadapt/main.nf'
include { UMITOOLS_DEDUP } from './modules/nf-core/umitools/dedup/main.nf'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_NAMESORT } from './modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FIXMATE } from './modules/nf-core/samtools/fixmate/main.nf'
include { SAMTOOLS_MARKDUP } from './modules/nf-core/samtools/markdup/main.nf'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW } from './modules/nf-core/samtools/view/main.nf'
include { SUBTRACTBKG } from './modules/local/subtractbkg/main.nf'

params {
    input_csv: Path
    ref_fasta: Path
    ref_gtf: Path
    star_ignore_sjdbgtf: Boolean
    chr_list: Path
    combination_mode: String
    control_combination: String
    comparisons_csv: Path
    min_control_coverage: Integer
    max_control_edit_frac: BigDecimal
    min_control_non_g_frac: BigDecimal
    min_rna_edit_count: Integer
    edit_threshold: BigDecimal
    umi: Boolean
    skip_umiextract: Boolean
    skip_trimming: Boolean
}

workflow {

    main:
    reads_ch = channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                single_end: row.fastq_2 == '',
                condition: row.condition,
            ]
            def files = meta.single_end
                ? [file(row.fastq_1)]
                : [file(row.fastq_1), file(row.fastq_2)]
            [meta, files]
        }

    FASTQC(reads_ch)

    reads_ch
        .branch {
            extract: params.umi && !params.skip_umiextract
            bypass: true
        }
        .set { routed_reads_ch }

    UMITOOLS_EXTRACT(routed_reads_ch.extract)
    post_umi_ch = UMITOOLS_EXTRACT.out.reads.mix(routed_reads_ch.bypass)

    post_umi_ch
        .branch {
            trim: !params.skip_trimming
            bypass: true
        }
        .set { routed_trim_ch }

    CUTADAPT(routed_trim_ch.trim)
    align_input_ch = CUTADAPT.out.reads.mix(routed_trim_ch.bypass)
    FASTQC_TRIMMED(CUTADAPT.out.reads)

    reffasta_ch = channel.fromPath(params.ref_fasta).map { fasta -> [[id: 'genome'], fasta] }
    refgtf_ch = channel.fromPath(params.ref_gtf).map { gtf -> [[id: 'genome'], gtf] }
    PREPARE_GENE_SPANS(
        channel.fromPath(params.ref_gtf)
    )

    STAR_GENOMEGENERATE(reffasta_ch, refgtf_ch)

    STAR_ALIGN(
        align_input_ch,
        STAR_GENOMEGENERATE.out.index.first(),
        refgtf_ch.first(),
        params.star_ignore_sjdbgtf ?: false,
    )

    star_all_logs = STAR_ALIGN.out.log_final.mix(
        STAR_ALIGN.out.log_out,
        STAR_ALIGN.out.log_progress,
    )

    multiqc_files = FASTQC.out.zip
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(CUTADAPT.out.log)
        .mix(star_all_logs)

    multiqc_input_ch = multiqc_files
        .collect { _meta, files -> files }
        .map { files -> [[id: 'multiqc'], files, [], [], [], []] }

    MULTIQC(multiqc_input_ch)

    STAR_ALIGN.out.bam_sorted_aligned
        .branch {
            umi: params.umi
            markdup: true
        }
        .set { routed_bam }

    SAMTOOLS_INDEX(routed_bam.umi)
    UMITOOLS_DEDUP(routed_bam.umi.join(SAMTOOLS_INDEX.out.index), false)

    SAMTOOLS_NAMESORT(routed_bam.markdup, [[], [], []], [])
    SAMTOOLS_FIXMATE(SAMTOOLS_NAMESORT.out.bam)
    SAMTOOLS_SORT(SAMTOOLS_FIXMATE.out.bam, [[], [], []], [])
    SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam, [[], [], []])

    dedup_bam_ch = UMITOOLS_DEDUP.out.bam.mix(SAMTOOLS_MARKDUP.out.bam)

    SAMTOOLS_VIEW(
        dedup_bam_ch.map { meta, bam -> [meta, bam, []] },
        [[], [], []],
        [[], []],
        [[], []],
        [],
    )

    SAM_TO_MAP(SAMTOOLS_VIEW.out.sam, params.chr_list)

    channel.fromPath(params.comparisons_csv)
        .splitCsv(header: true)
        .map { row -> row.control }
        .unique()
        .set { control_conditions_ch }

    SAM_TO_MAP.out.map
        .map { meta, bin -> [meta.condition, bin] }
        .combine(control_conditions_ch)
        .filter { cond, _bin, ctrl_cond -> cond == ctrl_cond }
        .map { cond, bin, _ctrl_cond -> [cond, bin] }
        .groupTuple()
        .set { raw_controls_grouped_ch }

    COMBINE_CONTROLS(
        raw_controls_grouped_ch,
        params.control_combination,
    )

    channel.fromPath(params.comparisons_csv)
        .splitCsv(header: true)
        .flatMap { row ->
            def pairs = [[row.treatment, row.control]]
            if (row.containsKey('background') && row.background) {
                pairs << [row.background, row.control]
            }
            return pairs
        }
        .unique()
        .set { comparison_pairs_ch }

    comparison_pairs_ch
        .combine(
            SAM_TO_MAP.out.map.map { meta, bin -> [meta.condition, meta, bin] }
        )
        .filter { treatment_cond, _control_cond, sample_cond, _meta, _bin ->
            treatment_cond == sample_cond
        }
        .map { _treatment_cond, control_cond, _sample_cond, meta, bin ->
            [control_cond, meta, bin]
        }
        .combine(COMBINE_CONTROLS.out.bin)
        .filter { control_cond, _meta, _trt_bin, ctrl_cond_out, _ctrl_bin ->
            control_cond == ctrl_cond_out
        }
        .map { control_cond, meta, trt_bin, _ctrl_cond_out, ctrl_bin ->
            def pair_meta = meta + [
                id: "${meta.id}_vs_${control_cond}",
                comparison_id: "${meta.condition}_vs_${control_cond}",
                control_cond: control_cond,
            ]
            [pair_meta, trt_bin, ctrl_bin]
        }
        .set { find_edit_input_ch }

    FIND_EDIT_SITES(
        find_edit_input_ch,
        PREPARE_GENE_SPANS.out.gtf.first(),
        params.chr_list,
        params.min_control_coverage,
        params.max_control_edit_frac,
        params.min_control_non_g_frac,
        params.min_rna_edit_count,
    )

    FIND_EDIT_SITES.out.bin
        .map { meta, bin -> [meta.comparison_id, bin] }
        .groupTuple()
        .set { combine_treatments_grouped_ch }

    COMBINE_EDIT_SITES(
        combine_treatments_grouped_ch,
        params.combination_mode,
        params.chr_list,
    )

    channel.fromPath(params.comparisons_csv)
        .splitCsv(header: true)
        .map { row ->
            def target_id = "${row.treatment}_vs_${row.control}"
            def bg_id = row.containsKey('background') && row.background
                ? "${row.background}_vs_${row.control}"
                : "NO_BG"
            [target_id, bg_id]
        }
        .set { bg_relations_ch }

    COMBINE_EDIT_SITES.out.tsv
        .map { condition, tsv -> [condition, [id: condition], tsv] }
        .set { combined_tsvs }

    bg_relations_ch
        .join(combined_tsvs, by: 0)
        .branch { _target_id, bg_id, _meta, _target_tsv ->
            no_bg: bg_id == "NO_BG"
            needs_bg: bg_id != "NO_BG"
        }
        .set { routed_bg }

    ch_no_bg = routed_bg.no_bg.map { _target_id, _bg_id, meta, target_tsv ->
        [meta, target_tsv, []]
    }

    ch_with_bg = routed_bg.needs_bg
        .map { _target_id, bg_id, meta, target_tsv -> [bg_id, meta, target_tsv] }
        .join(combined_tsvs.map { bg_id_key, _meta, tsv -> [bg_id_key, tsv] }, by: 0)
        .map { _bg_id, meta, target_tsv, bg_tsv -> [meta, target_tsv, bg_tsv] }

    SUBTRACTBKG(ch_no_bg.mix(ch_with_bg))

    SAMTOOLS_FAIDX(
        reffasta_ch.map { meta, fasta -> [meta, fasta, []] },
        true,
    )

    BEDTOOLS_INTERSECT(
        SUBTRACTBKG.out.tsv.combine(PREPARE_GENE_SPANS.out.gtf),
        SAMTOOLS_FAIDX.out.sizes.first(),
    )

    SUMMARIZE_EDIT_SITES(
        BEDTOOLS_INTERSECT.out.intersect,
        params.edit_threshold,
    )

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    multiqc_report = MULTIQC.out.report
    star_index = STAR_GENOMEGENERATE.out.index
    star_align_log = star_all_logs
    starsorted_bam = STAR_ALIGN.out.bam_sorted_aligned
    dedupped_sam = SAMTOOLS_VIEW.out.sam
    sam_maps = SAM_TO_MAP.out.map
    pooled_controls = COMBINE_CONTROLS.out.bin
    filtered_bins = FIND_EDIT_SITES.out.bin
    combined_tsvs = COMBINE_EDIT_SITES.out.tsv
    annotated_tsv = BEDTOOLS_INTERSECT.out.intersect
    summarized_tsv = SUMMARIZE_EDIT_SITES.out
}

output {
    fastqc_html {
        path "fastqc/html"
    }
    fastqc_zip {
        path "fastqc/zip"
    }
    multiqc_report {
        path "multiqc"
    }
    star_index {
        path "star/index"
    }
    star_align_log {
        path "star/align/log"
    }
    starsorted_bam {
        path "star/align/sortedbam"
    }
    dedupped_sam {
        path "dedupped"
    }
    sam_maps {
        path "temple/map"
    }
    pooled_controls {
        path "temple/pooled_controls"
    }
    filtered_bins {
        path "temple/filtered_replicate_bins"
    }
    combined_tsvs {
        path "temple/combined_tsvs/${params.combination_mode}"
    }
    annotated_tsv {
        path "temple/annotated_edit_sites"
    }
    summarized_tsv {
        path "temple/summarized_edit_sites/${params.edit_threshold}"
    }
}
