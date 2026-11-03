#!/usr/bin/env Rscript
# =============================================================================
# Leave-one-sample-out robustness of the module-drought correlations.
#
# With n = 24 a single influential sample can create a strong correlation, so
# every module-drought r is recomputed 24 times, each time leaving one sample
# out.  A result that only survives with all samples present is not credible.
#
# Writes: leave_one_out_GSE103278.tsv
# =============================================================================
suppressPackageStartupMessages(library(WGCNA))
OUT <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis/wgcna")
r <- readRDS(file.path(OUT, "wgcna_results.rds"))
x <- r$res103
expr <- x$exprGeneMat          # genes x samples
colors <- x$colors
tr <- x$traits$drought

MEof <- function(genes, keep) {
  ## use WGCNA's own eigengene definition so r_full reproduces module_trait_*.tsv
  sub <- t(expr[genes, keep, drop = FALSE])          # samples x genes
  suppressWarnings(
    as.numeric(moduleEigengenes(sub, rep("M", length(genes)))$eigengenes[, 1]))
}
rows <- list()
for (m in names(x$moduleSizes)) {
  genes <- rownames(expr)[colors == m]   # colours are positional (no names)
  if (length(genes) < 10) next
  full <- as.numeric(cor(MEof(genes, seq_along(tr)), tr))   # plain numeric, not 1x1 matrix
  loo <- vapply(seq_along(tr), function(i) {
    keep <- setdiff(seq_along(tr), i)
    cc <- as.numeric(cor(MEof(genes, keep), tr[keep]))
    if (sign(cc) != sign(full)) -abs(cc) else cc   # keep sign relative to full
  }, numeric(1))
  rows[[m]] <- data.frame(
    module = m, n_genes = length(genes), r_full = as.numeric(full),
    r_loo_min = min(loo), r_loo_max = max(loo), r_loo_median = median(loo),
    n_sign_flips = sum(sign(loo) != sign(full)),
    worst_sample = colnames(expr)[which.min(abs(loo))], row.names = NULL)
}
res <- do.call(rbind, rows)
res <- res[order(-abs(res$r_full)), ]
write.table(res, file.path(OUT, "leave_one_out_GSE103278.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\n[LOO] modules with any sign flip: %d / %d\n",
            sum(res$n_sign_flips > 0), nrow(res)))
