#!/usr/bin/env Rscript

# summarize_editing.R
# Summarizes A-to-G editing sites at transcript level from annotated TSV
# produced by bedtools intersect (find_edit_sites output + GTF annotation)
#
# Usage:
#   Rscript summarize_editing.R --input <tsv> --output <tsv> --threshold <float>
#

suppressPackageStartupMessages({
    library(data.table)
    library(purrr)
    library(tibble)
})


args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    keys <- args[seq(1, length(args), 2)]
    values <- args[seq(2, length(args), 2)]
    setNames(as.list(values), gsub("^--", "", keys))
}

opts <- parse_args(args)

if (!all(c("input", "output", "threshold") %in% names(opts))) {
    stop("Usage: Rscript summarize_editing.R --input <tsv> --output <tsv> --threshold <float>")
}

input_path <- opts$input
output_path <- opts$output
edit_threshold <- as.numeric(opts$threshold)

# TSV columns from find_edit_sites (cols 1-8) + GTF columns from bedtools (cols 9+)
# bedtools -wa -wb appends the full GTF line after the query columns
# GTF columns: chr, source, feature, start, end, score, strand, frame, attributes

col_names <- c(
    "chr", "site_start", "site_end",
    "rna_G", "rna_total", "ctrl_G", "ctrl_total", "rna_edit_frac",
    "gtf_chr", "gtf_source", "gtf_feature", "gtf_start", "gtf_end",
    "gtf_score", "gtf_strand", "gtf_frame", "gtf_attributes"
)

dt <- fread(
    input_path,
    sep        = "\t",
    header     = FALSE,
    col.names  = col_names,
    skip       = 1 # skip the # header line from find_edit_sites
)

cat("Total rows:", nrow(dt), "\n")

# We exclude introns by keeping only rows where the GTF feature is "exon".
# Using transcript-level GTF rows would double-count sites and inflate lengths.
# CDS rows are also excluded since we want full exonic context.

dt <- dt[gtf_feature == "exon"]

cat("Rows after exon filter:", nrow(dt), "\n")

# attributes column looks like:
# gene_id "Gsta3"; transcript_id "NM_001420195.1"; gene_name "Gsta3"; ...

parse_attr <- function(attrs, key) {
    pattern <- paste0(key, ' "([^"]+)"')
    m <- regmatches(attrs, regexpr(pattern, attrs))
    ifelse(length(m) == 0 || nchar(m) == 0, NA_character_,
        sub(pattern, "\\1", m)
    )
}

dt[, gene_name := map_chr(gtf_attributes, \(a) parse_attr(a, "gene_name"))]
dt[, transcript_id := map_chr(gtf_attributes, \(a) parse_attr(a, "transcript_id"))]

# Each unique (transcript_id, exon) interval contributes its length.
# We deduplicate exon intervals first since the same exon can appear
# multiple times if a site overlaps multiple transcripts.

exon_lengths <- dt[
    !is.na(transcript_id),
    .(exon_length = unique(gtf_end - gtf_start)),
    by = .(transcript_id, gtf_start, gtf_end)
][, .(transcript_length = sum(exon_length)), by = transcript_id]

# A single editing site can match multiple exon rows of the same transcript
# (e.g. overlapping exon annotations). Keep one row per (site, transcript).

dt_dedup <- unique(dt[, .(
    chr, site_start, site_end,
    rna_G, rna_total, ctrl_G, ctrl_total, rna_edit_frac,
    gtf_strand, gene_name, transcript_id
)])


summarize_transcript <- function(d) {
    # editing percentage: sum(G) / sum(G + other) across all sites
    total_G <- sum(d$rna_G)
    total_reads <- sum(d$rna_total)
    edit_pct <- if (total_reads > 0) round(total_G / total_reads, 4) else 0

    # edit sites as "chr:pos:edit_frac" comma separated for easy parsing
    sites <- paste(
        sprintf("%s:%d:%.4f", d$chr, d$site_start, d$rna_edit_frac),
        collapse = ","
    )

    tibble(
        gene_name         = d$gene_name[1],
        chr               = d$chr[1],
        strand            = d$gtf_strand[1],
        n_edit_sites      = nrow(d),
        total_rna_G       = total_G,
        total_rna_reads   = total_reads,
        editing_pct       = edit_pct,
        edit_sites        = sites
    )
}

result <- dt_dedup[
    !is.na(transcript_id),
    summarize_transcript(.SD),
    by = transcript_id
]

result <- merge(
    result,
    exon_lengths,
    by = "transcript_id",
    all.x = TRUE
)

setcolorder(result, c(
    "transcript_id", "gene_name", "chr", "strand",
    "transcript_length", "n_edit_sites",
    "total_rna_G", "total_rna_reads", "editing_pct",
    "edit_sites"
))

# sort by editing percentage descending
result <- result[order(gene_name, -editing_pct, transcript_id)]
result <- result[editing_pct >= edit_threshold]
fwrite(result, output_path, sep = "\t", quote = FALSE)

gene_result <- result[, .(
    n_transcripts   = .N,
    chr             = chr[1],
    strand          = strand[1],
    n_edit_sites    = sum(n_edit_sites),
    total_rna_G     = sum(total_rna_G),
    total_rna_reads = sum(total_rna_reads),
    editing_pct     = round(sum(total_rna_G) / sum(total_rna_reads), 4),
    edit_sites      = paste(unique(unlist(strsplit(edit_sites, ","))), collapse = ",")
), by = gene_name][order(gene_name, -editing_pct)]

gene_output <- sub("\\.tsv$", "_gene.tsv", output_path)
fwrite(gene_result, gene_output, sep = "\t", quote = FALSE)
