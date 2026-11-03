#!/usr/bin/env Rscript
# Print the key numbers needed for WGCNA_REPORT.md straight from the result files.
OUT <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis/wgcna")
r <- readRDS(file.path(OUT, "wgcna_results.rds"))
res103 <- r$res103; res128 <- r$res128; repl <- r$repl

for (nm in c("GSE103278", "GSE128048")) {
  x <- if (nm == "GSE103278") res103 else res128
  cat("\n########", nm, "########\n")
  cat("samples:", x$nSamples, " genes in network:", x$nGenes, "\n")
  cat("soft power:", x$power, " |", x$power_note, "\n")
  sft <- x$sft
  cat("SFT fit at chosen power: R2=", round(sft$SFT.R.sq[sft$Power == x$power], 3),
      " slope=", round(sft$slope[sft$Power == x$power], 3),
      " meanK=", round(sft$mean.k.[sft$Power == x$power], 1), "\n", sep = "")
  cat("max R2 over all powers:", round(max(sft$SFT.R.sq), 3),
      " at power", sft$Power[which.max(sft$SFT.R.sq)], "\n")
  cat("R2 at power 12 / 18 / 20:", round(sft$SFT.R.sq[sft$Power %in% c(12, 18, 20)], 3), "\n")
  cat("modules (non-grey):", length(x$moduleSizes), "\n")
  cat("module sizes:\n"); print(x$moduleSizes)
  cat("grey (unassigned) genes:", sum(x$colors == "grey"), "\n")
  cat("\nmodule-trait (top 12 by p):\n")
  mt <- x$modTrait[order(x$modTrait$p), ]
  print(head(mt, 12), row.names = FALSE)
  cat("\nmodules with FDR<0.05 for drought:\n")
  print(x$modTrait[x$modTrait$trait == "drought" & x$modTrait$fdr < 0.05, ], row.names = FALSE)
  cat("\nconsensus overlap (FDR<0.05):\n")
  print(x$cross[x$cross$fdr < 0.05, ], row.names = FALSE)
  cat("\nn modules with any consensus gene overlap:\n")
  cc <- x$cross[x$cross$gene_set == "drought_consensus_1197", ]
  print(cc[order(cc$p_hypergeom), ][1:min(10, nrow(cc)), ], row.names = FALSE)
  cat("\ngenes in network that are consensus:", sum(x$geneTab$is_consensus),
      " top100:", sum(x$geneTab$is_top100), "\n")
  cat("hub genes: total", nrow(x$hub), " | hub consensus:", sum(x$hub$is_consensus),
      " | hub top100:", sum(x$hub$is_top100), "\n")
  cat("\nhub genes that are in consensus set:\n")
  h <- x$hub[x$hub$is_consensus, c("gene", "module", "kME_module", "GS_drought", "is_top100")]
  print(h[order(h$module, -h$kME_module), ], row.names = FALSE)
}

cat("\n######## REPLICATION ########\n")
print(repl, row.names = FALSE)
cat("\nconcordant:", sum(repl$direction_concordant), "/", nrow(repl), "\n")
cat("significant in GSE128048 (p<0.05):", sum(repl$p_GSE128048 < 0.05),
    " of which concordant:", sum(repl$p_GSE128048 < 0.05 & repl$direction_concordant), "\n")

if (file.exists(file.path(OUT, "cold_crosstalk.rds"))) {
  cat("\n######## COLD ########\n")
  cc <- readRDS(file.path(OUT, "cold_crosstalk.rds"))
  print(head(cc$res[order(cc$res$p_hypergeom),
                    c("gene_set", "n_in_module", "n_in_universe", "overlap",
                      "fold_enrichment", "p_hypergeom", "fdr")], 20), row.names = FALSE)
}
