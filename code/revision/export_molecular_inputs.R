#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript export_molecular_inputs.R EVIDENCE_ROOT OUTPUT_DIR")
base <- file.path(args[1], "processed/reanalysis")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (dataset in c("GSE103278", "GSE128048")) {
  obj <- readRDS(file.path(base, paste0(dataset, ".rds")))
  fstat <- topTableF(obj$fit, number = Inf, adjust.method = "BH", sort.by = "F")
  fdr <- fstat[rownames(obj$fit$coefficients), "adj.P.Val"]
  values <- data.frame(probe = obj$expression$genes$ProbeName,
                       mean_logFC = rowMeans(obj$fit$coefficients), omnibus_FDR = fdr)
  write.table(values, file.path(out, paste0(dataset, "_volcano_values.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  con <- gzfile(file.path(out, paste0(dataset, "_normalized_expression.tsv.gz")), "wt")
  write.table(data.frame(probe = obj$expression$genes$ProbeName, obj$expression$E,
                         check.names = FALSE), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}
capture.output(sessionInfo(), file = file.path(out, "molecular_export_sessionInfo.txt"))
