#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tools)
    library(here)
})

new_path <- here("ChinaTown/temple/filtered_edit_sites/HyperTRIBE_vs_wtRNA.tsv")
leg_path <- here("workflow/contract/test/validation_data/present_both_wtRNA.bedgraph")

if (!file.exists(new_path)) {
    stop(sprintf("New filtered sites file not found at: %s", new_path))
}
if (!file.exists(leg_path)) {
    stop(sprintf("Legacy bedgraph file not found at: %s", leg_path))
}

new_dt <- read.delim(new_path, header = TRUE, stringsAsFactors = FALSE)
# Expecting columns: #chr, start, end, ... rna_edit_frac
colnames(new_dt)[1] <- "chr"
colnames(new_dt)[2] <- "start"
colnames(new_dt)[3] <- "end"

# Find editing fraction column dynamically
frac_col <- grep("edit_frac|fraction|pct", colnames(new_dt), ignore.case = TRUE)
if (length(frac_col) > 0) {
    new_dt$rna_edit_frac <- as.numeric(new_dt[, frac_col[1]])
} else {
    new_dt$rna_edit_frac <- as.numeric(new_dt[, 8]) # fallback
}

new_filtered <- new_dt[!is.na(new_dt$rna_edit_frac) & new_dt$rna_edit_frac >= 0.05, ]
leg_dt <- read.delim(leg_path, header = FALSE, stringsAsFactors = FALSE)

new_pos <- paste(new_filtered$chr, new_filtered$end, sep = "_")
leg_pos <- paste(leg_dt$V1, leg_dt$V2, sep = "_")

common <- intersect(new_pos, leg_pos)
pct_new_in_leg <- 100 * length(common) / length(new_pos)
pct_leg_in_new <- 100 * length(common) / length(leg_pos)

cat(sprintf("New filtered sites (>= 5%%) : %d\n", length(new_pos)))
cat(sprintf("Legacy total sites         : %d\n", length(leg_pos)))
cat(sprintf("Overlapping sites          : %d\n", length(common)))
cat(sprintf("%% of new found in legacy   : %.2f%%\n", pct_new_in_leg))
cat(sprintf("%% of legacy found in new   : %.2f%%\n\n", pct_leg_in_new))
