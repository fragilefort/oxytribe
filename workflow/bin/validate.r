#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(data.table)
    library(here)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    if (length(args) == 0) {
        return(list())
    }

    keys <- args[seq(1, length(args), 2)]
    values <- args[seq(2, length(args), 2)]

    setNames(as.list(values), gsub("^--", "", keys))
}

opts <- parse_args(args)

new_gene_path <- if (!is.null(opts$new_gene_summary)) {
    opts$new_gene_summary
} else {
    here("ChinaTown/temple/summarized_edit_sites/0.00/HyperTRIBE_vs_wtRNA_summary_gene.tsv")
}
legacy_xls_path <- if (!is.null(opts$legacy_xls)) {
    opts$legacy_xls
} else {
    here("workflow/contract/test/validation_data/")
}
top_n <- if (!is.null(opts$top_n)) as.integer(opts$top_n) else 30
out_dir <- if (!is.null(opts$out_dir)) opts$out_dir else "validation_output"

if (!file.exists(new_gene_path)) {
    stop(sprintf("New gene summary not found at: %s", new_gene_path))
}
if (!file.exists(legacy_xls_path)) {
    stop(sprintf("Legacy results file not found at: %s", legacy_xls_path))
}

# Named targets confirmed in published literature, independent of both
# pipelines' own output (McMahon et al. 2016, Xu et al. 2018 - Hrp48
# HyperTRIBE western blot / TRIBE-HyperTRIBE-CLIP consensus targets).
# Fixed list, not user-supplied: these are specific, citable claims, not
# a general-purpose input.
KNOWN_TARGETS <- c(
    "Rm62", "Syt1", "fne", "cwo", "unc-13",
    "CG31694", "CG11593", "Chd64"
)

legacy <- fread(legacy_xls_path, sep = "\t", header = TRUE, quote = "")
setnames(legacy, old = names(legacy)[1:2], new = c("gene_name", "n_edit_sites"))
legacy[, gene_name := trimws(gene_name)]
legacy <- legacy[gene_name != ""]

new_dt <- fread(new_gene_path, sep = "\t", header = TRUE, quote = "")
if (!"gene_name" %in% names(new_dt)) setnames(new_dt, old = names(new_dt)[1], new = "gene_name")
new_dt[, gene_name := trimws(gene_name)]
new_dt <- new_dt[gene_name != ""]

legacy_genes <- unique(legacy$gene_name)
new_genes <- unique(new_dt$gene_name)

shared <- intersect(legacy_genes, new_genes)
legacy_only <- setdiff(legacy_genes, new_genes)
new_only <- setdiff(new_genes, legacy_genes)

pct_recovered <- 100 * length(shared) / length(legacy_genes)

cat("==== Gene-level overlap ====\n")
cat(sprintf("Legacy genes         : %d\n", length(legacy_genes)))
cat(sprintf("New pipeline genes   : %d\n", length(new_genes)))
cat(sprintf("Shared               : %d\n", length(shared)))
cat(sprintf("Legacy-only (missed) : %d\n", length(legacy_only)))
cat(sprintf("New-only (novel)     : %d\n", length(new_only)))
cat(sprintf("%% of legacy recovered: %.2f%%\n\n", pct_recovered))

if (length(legacy_only) > 0) {
    cat("Legacy genes missing from new pipeline:\n")
    cat(paste(" -", sort(legacy_only)), sep = "\n")
    cat("\n")
}

legacy_top <- legacy[order(-n_edit_sites)][1:min(top_n, .N), gene_name]
new_top <- new_dt[order(-n_edit_sites)][1:min(top_n, .N), gene_name]
top_shared <- intersect(legacy_top, new_top)

cat(sprintf("==== Top-%d by editing site count ====\n", top_n))
cat(sprintf("Overlap: %d / %d\n\n", length(top_shared), top_n))

comparison_tbl <- data.table(
    rank = seq_len(top_n),
    legacy_gene = legacy_top[seq_len(top_n)],
    new_gene = new_top[seq_len(top_n)]
)
print(comparison_tbl)
cat("\n")

cat("==== Named target validation (literature-confirmed targets) ====\n")
for (g in KNOWN_TARGETS) {
    in_legacy <- g %in% legacy_genes
    row <- new_dt[gene_name == g]
    if (nrow(row) > 0) {
        cat(sprintf(
            "%-10s legacy:%-5s new: FOUND (%d sites)\n",
            g, in_legacy, row$n_edit_sites[1]
        ))
    } else {
        cat(sprintf("%-10s legacy:%-5s new: NOT FOUND\n", g, in_legacy))
    }
}
cat("\n")

# --- Write outputs ---
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(sort(shared), file.path(out_dir, "shared_genes.txt"))
writeLines(sort(legacy_only), file.path(out_dir, "legacy_only_genes.txt"))
writeLines(sort(new_only), file.path(out_dir, "new_only_genes.txt"))
fwrite(comparison_tbl, file.path(out_dir, "top_n_comparison.tsv"), sep = "\t")

cat(sprintf("Gene lists and top-%d comparison written to: %s/\n", top_n, out_dir))
