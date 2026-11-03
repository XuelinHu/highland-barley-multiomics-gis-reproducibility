#!/usr/bin/env Rscript
# =============================================================================
# Low-temperature cross-talk check
#
# GSE41513 (barley shoot/root cold time course, GPL16130) is only 18 samples and
# is provided as log2 ratios vs 0H, so it is NOT used to build a second WGCNA
# network.  Instead the drought co-expression modules found in GSE103278 are
# tested for enrichment among cold-responsive genes using the cold limma
# contrasts that already exist under processed/reanalysis/cold/.
#
# Reads : processed/reanalysis/cold/GSE41513_{shoot,root}_{3,6}H_limma.tsv
#         processed/mapping/GPL16130_morex_mapping.tsv
#         processed/reanalysis/wgcna/wgcna_results.rds
# Writes: processed/reanalysis/wgcna/cold_crosstalk.tsv
# =============================================================================
BASE <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed")
RE   <- file.path(BASE, "reanalysis")
OUT  <- file.path(RE, "wgcna")

## ---- probe -> Morex gene for GPL16130 --------------------------------------
mp <- read.delim(file.path(BASE, "mapping/GPL16130_morex_mapping.tsv"),
                 stringsAsFactors = FALSE)
mp <- mp[mp$status == "unique_gene", c("probe_id", "gene_ids")]
names(mp) <- c("probe", "gene")

## ---- cold-responsive gene sets ---------------------------------------------
contrasts <- c("shoot_3H", "shoot_6H", "root_3H", "root_6H")
cold_sets <- list()
resp_all  <- character(0)
for (cn in contrasts) {
  f <- file.path(RE, "cold", sprintf("GSE41513_%s_limma.tsv", cn))
  x <- read.delim(f, stringsAsFactors = FALSE, row.names = 1)
  x$probe <- rownames(x)
  x$gene  <- mp$gene[match(x$probe, mp$probe)]
  x <- x[!is.na(x$gene), ]
  ## one probe per gene: keep the strongest |logFC|
  x <- x[order(x$gene, -abs(x$logFC)), ]
  x <- x[!duplicated(x$gene), ]
  sig <- x$gene[x$adj.P.Val < 0.05 & abs(x$logFC) >= 1]
  cold_sets[[cn]] <- sig
  resp_all <- union(resp_all, sig)
  cat(sprintf("[cold] %-9s : %d genes tested, %d cold-responsive (FDR<0.05, |log2FC|>=1)\n",
              cn, nrow(x), length(sig)))
}
cat("[cold] union of cold-responsive genes:", length(resp_all), "\n")

## ---- module genes from the drought network ---------------------------------
rr <- readRDS(file.path(OUT, "wgcna_results.rds"))
gt <- rr$res103$geneTab
network_genes <- gt$gene
cons_genes <- gt$gene[gt$is_consensus]
t100_genes <- gt$gene[gt$is_top100]

## background = genes present in BOTH the drought network and the cold platform
cold_universe <- unique(mp$gene)
bg <- intersect(network_genes, cold_universe)
cat("[cold] background (drought network AND cold platform):", length(bg), "\n")

hyper <- function(set, label, universe = bg, module_genes = NULL) {
  if (is.null(module_genes)) module_genes <- universe
  module_genes <- intersect(module_genes, universe)
  S <- intersect(set, universe)
  N <- length(universe)                       # universe size
  n <- length(module_genes)                   # genes drawn (module, in universe)
  K <- length(S)                              # successes in universe
  k <- sum(module_genes %in% S)               # successes drawn
  p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  data.frame(gene_set = label, n_in_module = n, n_in_universe = K,
             overlap = k, expected = n * K / N,
             fold_enrichment = ifelse(n * K > 0, k / (n * K / N), NA_real_),
             p_hypergeom = p, row.names = NULL)
}

## Each row asks: are cold-responsive genes over-represented in this gene set?
## (the reverse framing - are drought genes enriched among cold genes - is the
##  same test, so the consensus / top100 sets are passed as the "module".)
rows <- list(
  hyper(resp_all, "whole_drought_network",   module_genes = network_genes),
  hyper(resp_all, "drought_consensus_1197",  module_genes = cons_genes),
  hyper(resp_all, "top100_drought_candidates", module_genes = t100_genes),
  hyper(resp_all, "consensus_or_top100",     module_genes = union(cons_genes, t100_genes))
)
for (m in names(rr$res103$moduleSizes)) {
  g <- gt$gene[gt$module == m]
  r <- hyper(resp_all, paste0("module_", m), module_genes = g)
  rows[[length(rows) + 1]] <- r
}
res <- do.call(rbind, rows)
res$fdr <- p.adjust(res$p_hypergeom, method = "BH")
res <- res[order(res$p_hypergeom), ]

write.table(res, file.path(OUT, "cold_crosstalk.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
cat("\n[cold] top enrichment results:\n")
print(head(res[, c("gene_set", "n_in_module", "n_in_universe", "overlap",
                   "fold_enrichment", "p_hypergeom", "fdr")], 15))

## ---- how many drought consensus genes respond to cold at all? --------------
cat(sprintf("\n[cold] drought consensus genes on cold platform: %d ; cold-responsive: %d (%.1f%%)\n",
            length(intersect(cons_genes, cold_universe)),
            length(intersect(cons_genes, resp_all)),
            100 * length(intersect(cons_genes, resp_all)) / max(1, length(intersect(cons_genes, cold_universe)))))
cat(sprintf("[cold] background cold-responsive rate: %.1f%%\n",
            100 * length(intersect(bg, resp_all)) / max(1, length(bg))))

## also save the per-module numbers for the report
saveRDS(list(res = res, cold_sets = cold_sets, resp_all = resp_all, bg = bg),
        file.path(OUT, "cold_crosstalk.rds"))
cat("[cold] done\n")
