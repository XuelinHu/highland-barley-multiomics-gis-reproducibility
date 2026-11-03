#!/usr/bin/env Rscript
# =============================================================================
# Are the hub genes of each module enriched for drought consensus / top100 genes?
#
# The module-level test (module_consensus_overlap_*.tsv) asks whether a whole
# module is enriched.  This script asks the sharper question: are the top-20
# hub genes (highest kME) of each module enriched?  n = 20 drawn from the
# 8000-gene network, which contains 1197 (15.0%) consensus and 100 top100 genes.
#
# Writes: hub_gene_enrichment_GSE103278.tsv / _GSE128048.tsv
# =============================================================================
OUT <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis/wgcna")
r <- readRDS(file.path(OUT, "wgcna_results.rds"))

for (nm in c("GSE103278", "GSE128048")) {
  x  <- if (nm == "GSE103278") r$res103 else r$res128
  gt <- x$geneTab; hub <- x$hub
  N  <- nrow(gt)
  Kc <- sum(gt$is_consensus); Kt <- sum(gt$is_top100)

  rows <- do.call(rbind, lapply(sort(unique(hub$module)), function(m) {
    h <- hub$gene[hub$module == m]; n <- length(h)
    kc <- sum(h %in% gt$gene[gt$is_consensus])
    kt <- sum(h %in% gt$gene[gt$is_top100])
    data.frame(module = m, n_hub = n,
               hub_consensus = kc, expected_consensus = n * Kc / N,
               fold_consensus = kc / (n * Kc / N),
               p_consensus = phyper(kc - 1, Kc, N - Kc, n, lower.tail = FALSE),
               hub_top100 = kt, expected_top100 = n * Kt / N,
               fold_top100 = kt / (n * Kt / N),
               p_top100 = phyper(kt - 1, Kt, N - Kt, n, lower.tail = FALSE),
               row.names = NULL)
  }))
  rows$fdr_consensus <- p.adjust(rows$p_consensus, method = "BH")
  rows$fdr_top100    <- p.adjust(rows$p_top100, method = "BH")
  rows <- rows[order(rows$p_consensus), ]
  write.table(rows, file.path(OUT, paste0("hub_gene_enrichment_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  cat(sprintf("\n#### %s : network N = %d | consensus in network = %d (%.1f%%) | top100 = %d\n",
              nm, N, Kc, 100 * Kc / N, Kt))
  print(head(rows, 8), row.names = FALSE, digits = 3)
  cat(sprintf("ALL %d hub genes: %d consensus (%.1f%%, p = %.2e) ; %d top100 (p = %.2e)\n",
              nrow(hub), sum(hub$is_consensus), 100 * sum(hub$is_consensus) / nrow(hub),
              phyper(sum(hub$is_consensus) - 1, Kc, N - Kc, nrow(hub), lower.tail = FALSE),
              sum(hub$is_top100),
              phyper(sum(hub$is_top100) - 1, Kt, N - Kt, nrow(hub), lower.tail = FALSE)))
}
cat("\n[done] hub_gene_enrichment_*.tsv\n")
