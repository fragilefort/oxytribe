#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
    keys <- args[seq(1, length(args), 2)]
    values <- args[seq(2, length(args), 2)]
    setNames(as.list(values), gsub("^--", "", keys))
}
opts <- parse_args(args)

if (!"input" %in% names(opts) || !"output" %in% names(opts)) {
    stop("Usage: Rscript subtract_background.R --input <target.tsv> --output <filtered.tsv> [--bg <bg.tsv>]")
}

# Safe reader returns an empty data.table with correct columns if file is missing/empty
safe_read <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) {
        return(data.table(V1 = character(), V2 = integer(), V3 = integer()))
    }
    fread(path, sep = "\t", header = FALSE, skip = "#")
}

dt <- safe_read(opts$input)

if ("bg" %in% names(opts) && nrow(dt) > 0) {
    dt_bg <- safe_read(opts$bg)

    if (nrow(dt_bg) > 0) {
        # Anti-join: Keep rows in dt that do NOT match chr(V1), start(V2), end(V3) in dt_bg
        dt <- dt[!dt_bg, on = c("V1", "V2", "V3")]
    }
}

fwrite(dt, opts$output, sep = "\t", col.names = FALSE, quote = FALSE)
