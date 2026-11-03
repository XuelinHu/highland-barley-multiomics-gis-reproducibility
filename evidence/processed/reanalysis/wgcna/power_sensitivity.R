#!/usr/bin/env Rscript
# =============================================================================
# Soft-threshold power sensitivity check.
#
# The scale-free-topology criterion picks power = 26 for both GSE103278 and
# GSE128048, which is unusually high.  With only 24 samples the R^2 curve does
# not plateau (it keeps climbing to 0.91 at power 30), so the "first power above
# 0.85" rule is not sharply defined.  The WGCNA FAQ recommends power 18 for
# signed networks with 20-30 samples, so this script rebuilds the GSE103278
# network at power 12 / 18 / 26 and asks whether the biological readout -
# the drought-associated modules and their drought-consensus-gene enrichment -
# is stable across that choice.
#
# Writes: power_sensitivity.tsv  (one row per power)
#         power_sensitivity_modules.tsv (per-module drought stats at each power)
# =============================================================================
suppressPackageStartupMessages({ library(WGCNA); library(limma) })
options(stringsAsFactors = FALSE)
disableWGCNAThreads()

BASE <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis")
OUT  <- file.path(BASE, "wgcna")
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25
POWERS <- c(26, 12, 18)   # power 26 first: it is the reference solution

mp <- read.delim(file.path(BASE, "../mapping/GPL14927_morex_mapping.tsv"))
mp <- mp[mp$status == "unique_gene", c("probe_id", "gene_ids")]
names(mp) <- c("probe", "gene")
cons <- unique(read.delim(file.path(BASE, "drought_consensus_genes.tsv"))$gene_id)
t100 <- unique(read.delim(file.path(BASE, "top100_drought_candidates.tsv"))$gene_id)

obj <- readRDS(file.path(BASE, "GSE103278.rds"))
el <- obj$expression; E <- el$E; pn <- el$genes$ProbeName
v <- apply(E, 1, var); ord <- order(pn, -v); keep <- ord[!duplicated(pn[ord])]
E <- E[keep, , drop = FALSE]; pn <- pn[keep]
d <- obj$design; grp <- colnames(d)[apply(d, 1, which.max)]
info <- data.frame(sample = colnames(E),
                   genotype = sub("^(CamB|Maresi)_.*$", "\\1", grp),
                   tissue   = sub("^[^_]+_([a-z]+)_.*$", "\\1", grp),
                   treatment= sub("^.*_([a-z]+)$", "\\1", grp),
                   row.names = colnames(E))
sel <- pn %in% mp$probe; Es <- E[sel, , drop = FALSE]
g <- mp$gene[match(pn[sel], mp$probe)]; vs <- apply(Es, 1, var)
o <- order(g, -vs); Es <- Es[o, ]; g <- g[o]; k <- !duplicated(g)
Eg <- Es[k, , drop = FALSE]; rownames(Eg) <- g[k]

vr <- apply(Eg, 1, var)
top <- names(sort(vr, decreasing = TRUE))[1:8000]
datExpr <- t(Eg[top, , drop = FALSE])
traits <- data.frame(drought  = as.integer(info$treatment == "drought"),
                     genotype = as.integer(info$genotype == "Maresi"),
                     tissue   = as.integer(info$tissue == "root"),
                     row.names = rownames(datExpr))[rownames(datExpr), , drop = FALSE]

res <- list(); modrows <- list()
ref_up <- NULL; ref_dn <- NULL
for (p in POWERS) {
  set.seed(42)
  net <- blockwiseModules(datExpr, power = p, networkType = "signed",
                          TOMType = "signed", minModuleSize = MIN_MODULE_SIZE,
                          mergeCutHeight = MERGE_CUT_HEIGHT, deepSplit = 2,
                          pamRespectsDendro = FALSE, reassignThreshold = 1e-6,
                          numericLabels = TRUE, maxBlockSize = ncol(datExpr) + 100,
                          verbose = 0, saveTOMFileBase = NA)
  colors <- labels2colors(net$colors)
  mods <- sort(unique(colors[colors != "grey"]))
  MEs <- orderMEs(moduleEigengenes(datExpr, colors)$eigengenes)
  mt <- do.call(rbind, lapply(c("drought", "genotype", "tissue"), function(tr) {
    ct <- cor(MEs, traits[[tr]], use = "p")
    data.frame(module = sub("^ME", "", names(MEs)), trait = tr,
               cor = as.numeric(ct),
               p = as.numeric(corPvalueStudent(ct, nrow(traits))), row.names = NULL)
  }))
  mt$fdr <- p.adjust(mt$p, method = "BH")
  dr <- mt[mt$trait == "drought", ]; dr <- dr[order(-abs(dr$cor)), ]

  ## The two poles of the drought response.  NOTE: WGCNA module colours are
  ## assigned by size rank, so the SAME colour means different gene sets at
  ## different powers - colours are NOT stable identifiers.  We therefore track
  ## the up- and the down-in-drought module by their correlation sign.
  N <- ncol(datExpr); K <- sum(colnames(datExpr) %in% cons)
  enr <- function(genes) {
    n <- length(genes); kk <- sum(genes %in% cons)
    c(overlap = kk, expected = n * K / N, fold = kk / (n * K / N),
      p = phyper(kk - 1, K, N - K, n, lower.tail = FALSE))
  }
  up <- dr$module[which.max(dr$cor)]
  dn <- dr$module[which.min(dr$cor)]
  gup <- colnames(datExpr)[colors == up]
  gdn <- colnames(datExpr)[colors == dn]
  eup <- enr(gup); edn <- enr(gdn)

  jac <- function(a, b) if (is.null(a)) NA_real_ else
    length(intersect(a, b)) / length(union(a, b))
  if (is.null(ref_up)) { ref_up <- gup; ref_dn <- gdn; j_up <- j_dn <- NA_real_ }
  else                 { j_up <- jac(ref_up, gup); j_dn <- jac(ref_dn, gdn) }

  cat(sprintf(paste0("[power %2d] %2d modules, grey=%d | UP module %s r=%+.3f (n=%d, ",
                     "consensus %d vs %.1f exp, %.2fx, p=%.1e) | DOWN module %s r=%+.3f ",
                     "(n=%d, consensus %d vs %.1f exp, %.2fx, p=%.1e)\n"),
              p, length(mods), sum(colors == "grey"), up, dr$cor[which.max(dr$cor)],
              length(gup), eup["overlap"], eup["expected"], eup["fold"], eup["p"],
              dn, dr$cor[which.min(dr$cor)], length(gdn), edn["overlap"],
              edn["expected"], edn["fold"], edn["p"]))
  res[[as.character(p)]] <- data.frame(
    power = p, n_modules = length(mods), n_grey = sum(colors == "grey"),
    n_drought_fdr05 = sum(dr$fdr < 0.05),
    up_module = up, up_r = dr$cor[which.max(dr$cor)], up_size = length(gup),
    up_consensus_overlap = eup["overlap"], up_consensus_fold = eup["fold"],
    up_consensus_p = eup["p"],
    down_module = dn, down_r = dr$cor[which.min(dr$cor)], down_size = length(gdn),
    down_consensus_overlap = edn["overlap"], down_consensus_fold = edn["fold"],
    down_consensus_p = edn["p"],
    jaccard_up_vs_power26 = j_up, jaccard_down_vs_power26 = j_dn, row.names = NULL)
  modrows[[as.character(p)]] <- cbind(power = p, dr)
}
write.table(do.call(rbind, res), file.path(OUT, "power_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(do.call(rbind, modrows), file.path(OUT, "power_sensitivity_modules.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("\n[done] power_sensitivity.tsv\n")
