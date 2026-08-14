#!/usr/bin/env Rscript

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

gene_output <- sub("\\.tsv$", "_gene.tsv", output_path)

col_names <- c(
    "chr", "site_start", "site_end",
    "rna_G", "rna_total", "ctrl_G", "ctrl_total", "rna_edit_frac",
    "gtf_chr", "gtf_source", "gtf_feature", "gtf_start", "gtf_end",
    "gtf_score", "gtf_strand", "gtf_frame", "gtf_attributes"
)

write_empty_outputs <- function(tx_out, gene_out) {
    empty_tx <- data.table(
        transcript_id = character(), gene_name = character(), chr = character(),
        strand = character(), transcript_length = integer(), n_edit_sites = integer(),
        total_rna_G = integer(), total_rna_reads = integer(), editing_pct = numeric(),
        edit_sites = character()
    )
    empty_gene <- data.table(
        gene_name = character(), n_transcripts = integer(), chr = character(),
        strand = character(), n_edit_sites = integer(), total_rna_G = integer(),
        total_rna_reads = integer(), editing_pct = numeric(), edit_sites = character()
    )
    fwrite(empty_tx, tx_out, sep = "\t", quote = FALSE)
    fwrite(empty_gene, gene_out, sep = "\t", quote = FALSE)
}

# Check if file is empty
if (!file.exists(input_path) || file.info(input_path)$size == 0) {
    cat("Input file is empty. Writing empty outputs...\n")
    write_empty_outputs(output_path, gene_output)
    quit(status = 0)
}

# Read lines and filter out any '#' header lines dynamically
raw_lines <- readLines(input_path)
clean_lines <- raw_lines[!grepl("^#", raw_lines)]

if (length(clean_lines) == 0) {
    cat("No data rows found after stripping comments. Writing empty outputs...\n")
    write_empty_outputs(output_path, gene_output)
    quit(status = 0)
}

# Parse clean lines into data.table
dt <- fread(text = clean_lines, sep = "\t", header = FALSE, col.names = col_names)

cat("Total rows loaded:", nrow(dt), "\n")

dt[, rna_edit_frac := as.numeric(rna_edit_frac)]
dt <- dt[rna_edit_frac >= edit_threshold]
cat("Rows after site-level threshold filter (>= ", edit_threshold, "): ", nrow(dt), "\n", sep = "")

if (nrow(dt) == 0) {
    cat("No matching edit sites found. Writing empty outputs...\n")
    write_empty_outputs(output_path, gene_output)
    quit(status = 0)
}

parse_attr <- function(attrs, key) {
    pattern <- paste0(".*?", key, ' "([^"]+)".*')
    has_key <- grepl(paste0(key, ' "([^"]+)"'), attrs)
    res <- sub(pattern, "\\1", attrs)
    res[!has_key] <- NA_character_
    res
}

dt[, gene_name := parse_attr(gtf_attributes, "gene_name")]
dt[, transcript_id := parse_attr(gtf_attributes, "transcript_id")]

# 1-based inclusive feature length calculation (+1)
exon_lengths <- dt[
    !is.na(transcript_id),
    .(exon_length = unique(gtf_end - gtf_start + 1)),
    by = .(transcript_id, gtf_start, gtf_end)
][, .(transcript_length = sum(exon_length)), by = transcript_id]

dt_dedup <- unique(dt[, .(
    chr, site_start, site_end,
    rna_G, rna_total, ctrl_G, ctrl_total, rna_edit_frac,
    gtf_strand, gene_name, transcript_id
)])

summarize_transcript <- function(d) {
    total_G <- sum(d$rna_G)
    total_reads <- sum(d$rna_total)
    edit_pct <- round(mean(d$rna_edit_frac), 4)

    sites <- paste(
        sprintf("%s:%d:%.4f", d$chr, d$site_end, d$rna_edit_frac),
        collapse = ","
    )

    tibble(
        gene_name        = d$gene_name[1],
        chr              = d$chr[1],
        strand           = d$gtf_strand[1],
        n_edit_sites     = nrow(d),
        total_rna_G      = total_G,
        total_rna_reads  = total_reads,
        editing_pct      = edit_pct,
        edit_sites       = sites
    )
}

result <- dt_dedup[
    !is.na(transcript_id),
    summarize_transcript(.SD),
    by = transcript_id
]

result <- merge(result, exon_lengths, by = "transcript_id", all.x = TRUE)

setcolorder(result, c(
    "transcript_id", "gene_name", "chr", "strand",
    "transcript_length", "n_edit_sites",
    "total_rna_G", "total_rna_reads", "editing_pct",
    "edit_sites"
))

result <- result[order(gene_name, -n_edit_sites, -editing_pct)]
result <- result[editing_pct >= edit_threshold]
fwrite(result, output_path, sep = "\t", quote = FALSE)

# Gene-level summary deduplicated at genomic site level
gene_sites <- unique(dt[!is.na(gene_name), .(
    chr, site_start, site_end, rna_G, rna_total, rna_edit_frac, gtf_strand, gene_name
)])

gene_result <- gene_sites[, .(
    n_transcripts   = uniqueN(dt[gene_name == .BY$gene_name & !is.na(transcript_id), transcript_id]),
    chr             = chr[1],
    strand          = gtf_strand[1],
    n_edit_sites    = .N,
    total_rna_G     = sum(rna_G),
    total_rna_reads = sum(rna_total),
    editing_pct     = round(mean(rna_edit_frac), 4),
    edit_sites      = paste(sprintf("%s:%d", chr, site_end), collapse = ",")
), by = gene_name][order(-n_edit_sites, -editing_pct)]

gene_result <- gene_result[editing_pct >= edit_threshold]
fwrite(gene_result, gene_output, sep = "\t", quote = FALSE)
