#!/usr/bin/env Rscript
# =============================================================================
# WGCNA co-expression network analysis of barley drought / low-temperature
# expression data (GSE103278 discovery, GSE128048 validation)
#
# Inputs  (read-only):
#   processed/reanalysis/GSE103278.rds, GSE128048.rds   (limma list w/ EList)
#   processed/mapping/GPL14927_morex_mapping.tsv         (probe -> Morex gene)
#   processed/reanalysis/drought_consensus_genes.tsv     (1197 consensus genes)
#   processed/reanalysis/top100_drought_candidates.tsv   (100 top candidates)
# Outputs (new dir): processed/reanalysis/wgcna/
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)

BASE <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis")
OUT  <- file.path(BASE, "wgcna")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_TOP_GENES     <- 8000
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25
R2_CUT          <- 0.85
R2_CUT_FALLBACK <- 0.80

logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
logmsg("WGCNA version: ", as.character(packageVersion("WGCNA")))

## ---------------------------------------------------------------- load data
mp <- read.delim(file.path(BASE, "../mapping/GPL14927_morex_mapping.tsv"))
mp <- mp[mp$status == "unique_gene", c("probe_id", "gene_ids")]
names(mp) <- c("probe", "gene")
logmsg("unique_gene probes in mapping: ", nrow(mp),
       " ; unique genes: ", length(unique(mp$gene)))

cons  <- read.delim(file.path(BASE, "drought_consensus_genes.tsv"))
top100 <- read.delim(file.path(BASE, "top100_drought_candidates.tsv"))
cons_genes   <- unique(cons$gene_id)
top100_genes <- unique(top100$gene_id)
logmsg("consensus genes: ", length(cons_genes), " ; top100: ", length(top100_genes))

#' Build gene-level expression matrix for one GEO series
#' @return list(mat = genes x samples, info = data.frame sample phenotype)
prep_dataset <- function(prefix) {
  obj <- readRDS(file.path(BASE, paste0(prefix, ".rds")))
  el  <- obj$expression          # limma EList
  E   <- el$E                    # 43803 probes x n samples, rownames are NULL
  pn  <- el$genes$ProbeName      # probe IDs in matrix row order
  stopifnot(length(pn) == nrow(E))

  ## recover probe IDs, drop duplicated probe rows (keep the higher-variance one)
  v <- apply(E, 1, var)
  ord <- order(pn, -v)
  keeprows <- ord[!duplicated(pn[ord])]
  E <- E[keeprows, , drop = FALSE]
  pn <- pn[keeprows]

  ## sample phenotype comes from the limma design matrix (one-hot, row = sample)
  d <- obj$design
  grp <- colnames(d)[apply(d, 1, which.max)]
  gt  <- sub("^(CamB|Maresi)_.*$", "\\1", grp)
  ts  <- sub("^[^_]+_([a-z]+)_.*$", "\\1", grp)
  tr  <- sub("^.*_([a-z]+)$", "\\1", grp)
  info <- data.frame(sample = colnames(E), group = grp,
                     genotype = gt, tissue = ts, treatment = tr)
  rownames(info) <- info$sample

  ## collapse probes -> genes, keeping the max-variance probe per gene
  sel <- pn %in% mp$probe
  Esub <- E[sel, , drop = FALSE]
  g    <- mp$gene[match(pn[sel], mp$probe)]
  vs   <- apply(Esub, 1, var)
  o    <- order(g, -vs)
  Esub <- Esub[o, , drop = FALSE]
  g    <- g[o]
  k    <- !duplicated(g)
  Eg   <- Esub[k, , drop = FALSE]
  rownames(Eg) <- g[k]
  logmsg(prefix, ": probes=", nrow(E), " -> mapped probes=", sum(sel),
         " -> genes=", nrow(Eg))
  list(mat = Eg, info = info)
}

## ------------------------------------------------------------ WGCNA pipeline
#' Run the full WGCNA workflow on a gene-level matrix
#' @param Eg genes x samples matrix
#' @param info sample annotation aligned to columns of Eg
#' @param tag dataset label used for filenames / plots
run_wgcna <- function(Eg, info, tag) {

  ## ---- 1. variance selection -------------------------------------------
  v <- apply(Eg, 1, var)
  nsel <- min(N_TOP_GENES, length(v))
  top  <- names(sort(v, decreasing = TRUE))[1:nsel]
  datExpr <- t(Eg[top, , drop = FALSE])          # samples x genes (WGCNA layout)
  logmsg(tag, ": selected ", ncol(datExpr), " genes by variance")

  ## ---- 2. sample outliers ----------------------------------------------
  gsg  <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    logmsg(tag, ": removed ", sum(!gsg$goodGenes), " bad genes / ",
           sum(!gsg$goodSamples), " bad samples")
  } else logmsg(tag, ": goodSamplesGenes OK")

  ## ---- 3. soft-thresholding power --------------------------------------
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30),
                           networkType = "signed", RsquaredCut = R2_CUT, verbose = 0)
  tab <- sft$fitIndices
  ok  <- tab$SFT.R.sq >= R2_CUT
  if (any(ok)) {
    power <- tab$Power[which(ok)[1]]
    power_note <- sprintf("first power with scale-free R^2 >= %.2f", R2_CUT)
  } else {
    ok2 <- tab$SFT.R.sq >= R2_CUT_FALLBACK
    if (any(ok2)) {
      power <- tab$Power[which(ok2)[1]]
      power_note <- sprintf("no power reached R^2 >= %.2f; used R^2 >= %.2f fallback",
                            R2_CUT, R2_CUT_FALLBACK)
    } else {
      power <- tab$Power[which.max(tab$SFT.R.sq)]
      power_note <- sprintf("no power reached R^2 >= %.2f; used power maximising R^2 (%.2f)",
                            R2_CUT_FALLBACK, max(tab$SFT.R.sq))
    }
  }
  logmsg(tag, ": soft power = ", power, " (", power_note, ")")

  ## ---- 4. network construction + module detection ----------------------
  net <- blockwiseModules(
    datExpr, power = power,
    networkType = "signed", TOMType = "signed",
    minModuleSize = MIN_MODULE_SIZE, mergeCutHeight = MERGE_CUT_HEIGHT,
    deepSplit = 2, pamRespectsDendro = FALSE, reassignThreshold = 1e-6,
    numericLabels = TRUE, maxBlockSize = ncol(datExpr) + 100,
    verbose = 0, saveTOMFileBase = NA)

  colors <- labels2colors(net$colors)
  mods   <- sort(unique(colors[colors != "grey"]))
  logmsg(tag, ": ", length(mods), " modules (non-grey); sizes: ",
         paste(mods, as.integer(table(colors[colors != "grey"])[mods]), sep = "=",
               collapse = ", "))

  ## ---- 5. module eigengenes & module-trait correlation -----------------
  MEs0 <- moduleEigengenes(datExpr, colors)$eigengenes
  MEs  <- orderMEs(MEs0)

  traits <- data.frame(
    drought  = as.integer(info$treatment == "drought"),
    genotype = as.integer(info$genotype  == "Maresi"),
    tissue   = as.integer(info$tissue    == "root"),
    row.names = info$sample)
  traits <- traits[rownames(datExpr), , drop = FALSE]

  modTrait <- function(MEs, traits) {
    out <- do.call(rbind, lapply(names(traits), function(tr) {
      ct <- cor(MEs, traits[[tr]], use = "p")
      cp <- corPvalueStudent(ct, nrow(traits))
      data.frame(module = sub("^ME", "", names(MEs)), trait = tr,
                 cor = as.numeric(ct), p = as.numeric(cp), row.names = NULL)
    }))
    out$fdr <- p.adjust(out$p, method = "BH")
    out
  }
  mt <- modTrait(MEs, traits)

  ## ---- 6. module membership / gene significance / hub genes ------------
  kME <- as.data.frame(signedKME(datExpr, MEs, corFnc = "cor", corOptions = "use = 'p'"))
  names(kME) <- sub("^kME", "", names(kME))

  gene_tab <- data.frame(gene = colnames(datExpr), module = colors, row.names = NULL)
  gene_tab$kME_module <- sapply(seq_len(nrow(gene_tab)), function(i) {
    m <- gene_tab$module[i]
    if (m == "grey") return(NA_real_)
    kME[i, m]
  })
  gene_tab$variance <- apply(datExpr, 2, var)
  gene_tab$GS_drought  <- cor(datExpr, traits$drought,  use = "p")
  gene_tab$GS_genotype <- cor(datExpr, traits$genotype, use = "p")
  gene_tab$GS_tissue   <- cor(datExpr, traits$tissue,   use = "p")
  gene_tab$is_consensus <- gene_tab$gene %in% cons_genes
  gene_tab$is_top100    <- gene_tab$gene %in% top100_genes

  hub <- do.call(rbind, lapply(mods, function(m) {
    sub <- gene_tab[gene_tab$module == m, ]
    sub <- sub[order(-sub$kME_module), ]
    head(sub, 20)
  }))
  hub$kME_rank_in_module <- ave(seq_len(nrow(hub)), hub$module, FUN = seq_along)

  ## ---- 7. intersection with consensus / top100 + hypergeometric test ---
  N <- nrow(gene_tab)
  cross <- do.call(rbind, lapply(mods, function(m) {
    inmod <- gene_tab$gene[gene_tab$module == m]
    n <- length(inmod)
    do <- function(set, label) {
      K <- sum(gene_tab$gene %in% set)            # set members in network
      k <- sum(inmod %in% set)                    # set members in module
      p <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
      data.frame(module = m, module_size = n, gene_set = label,
                 set_in_network = K, overlap = k,
                 expected = n * K / N,
                 fold_enrichment = ifelse(n * K / N > 0, k / (n * K / N), NA_real_),
                 p_hypergeom = p, row.names = NULL)
    }
    rbind(do(cons_genes, "drought_consensus_1197"),
          do(top100_genes, "top100_drought_candidates"))
  }))
  cross$fdr <- p.adjust(cross$p_hypergeom, method = "BH")

  ## ---- 8. plots --------------------------------------------------------
  ## 8a. module-trait heatmap
  pretty <- c(drought = "Drought vs Control", genotype = "Genotype (Maresi vs CamB)",
              tissue = "Tissue (Root vs Leaf)")
  M <- matrix(NA_real_, nrow = ncol(MEs), ncol = ncol(traits),
              dimnames = list(sub("^ME", "", names(MEs)), pretty[colnames(traits)]))
  P <- M
  for (i in seq_len(nrow(mt))) {
    M[mt$module[i], pretty[mt$trait[i]]] <- mt$cor[i]
    P[mt$module[i], pretty[mt$trait[i]]] <- mt$p[i]
  }
  txt <- matrix(sprintf("%.2f\n(%.1e)", M, P), nrow = nrow(M), dimnames = dimnames(M))

  png(file.path(OUT, paste0("module_trait_heatmap_", tag, ".png")),
      width = 1500, height = max(600, 110 * nrow(M) + 300), res = 130)
  par(mar = c(7, 11, 4, 2))
  labeledHeatmap(M, xLabels = colnames(M), yLabels = rownames(M),
                 ySymbols = rownames(M), colorLabels = FALSE,
                 colors = blueWhiteRed(50), textMatrix = txt,
                 setStdMargins = FALSE, cex.text = 0.75, zlim = c(-1, 1),
                 main = paste0("Module-trait relationships (", tag, ")\n",
                               "cell: Pearson r and p-value; n = ", nrow(datExpr), " samples"))
  dev.off()

  ## 8b. soft-threshold diagnostic
  png(file.path(OUT, paste0("soft_threshold_", tag, ".png")), width = 1300, height = 600, res = 130)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  plot(tab$Power, tab$SFT.R.sq, xlab = "Soft-threshold power", ylab = "Scale-free topology R^2",
       main = paste0(tag, ": scale independence"), pch = 19, col = "steelblue")
  abline(h = R2_CUT, col = "red", lty = 2); abline(h = R2_CUT_FALLBACK, col = "orange", lty = 3)
  abline(v = power, col = "green4", lty = 2)
  plot(tab$Power, tab$mean.k., xlab = "Soft-threshold power", ylab = "Mean connectivity",
       main = paste0(tag, ": mean connectivity"), pch = 19, col = "darkorange")
  abline(v = power, col = "green4", lty = 2)
  dev.off()

  ## 8c. module dendrogram with colours
  png(file.path(OUT, paste0("module_dendrogram_", tag, ".png")), width = 1500, height = 800, res = 130)
  plotDendroAndColors(net$dendrograms[[1]],
                      colors[net$blockGenes[[1]]], "Module colours",
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                      guideHang = 0.05, main = paste0(tag, " gene dendrogram (power = ", power, ")"))
  dev.off()

  ## 8d. module size barplot
  sz <- as.integer(table(colors[colors != "grey"])[mods])
  png(file.path(OUT, paste0("module_sizes_", tag, ".png")), width = 1300, height = 700, res = 130)
  par(mar = c(7, 4.5, 3, 1))
  bp <- barplot(sz, names.arg = mods, las = 2, col = mods,
                main = paste0("Module sizes (", tag, ")"),
                ylab = "Number of genes", xlab = "")
  text(bp, sz, labels = sz, pos = 3, cex = 0.7, xpd = NA)
  dev.off()

  ## 8e. sample clustering QC (n = 24 is small, so sample structure matters)
  png(file.path(OUT, paste0("sample_dendrogram_", tag, ".png")), width = 1200, height = 700, res = 130)
  par(mar = c(3, 5, 3, 1))
  sTree <- hclust(as.dist(1 - cor(t(datExpr), use = "p")), method = "average")
  lab <- paste0(sTree$labels, "\n", info$group[match(sTree$labels, info$sample)])
  plot(sTree, labels = lab, cex = 0.6, main = paste0(tag, ": sample clustering"),
       xlab = "", sub = "")
  dev.off()

  list(datExpr = datExpr, colors = colors, MEs = MEs, traits = traits,
       modTrait = mt, geneTab = gene_tab, hub = hub, cross = cross,
       power = power, power_note = power_note, sft = tab,
       moduleSizes = setNames(sz, mods), tag = tag,
       nSamples = nrow(datExpr), nGenes = ncol(datExpr),
       exprGeneMat = Eg[top, , drop = FALSE])
}

## =========================================================== discovery set
d103 <- prep_dataset("GSE103278")
res103 <- run_wgcna(d103$mat, d103$info, "GSE103278")

## =========================================================== validation set
d128 <- prep_dataset("GSE128048")
res128 <- run_wgcna(d128$mat, d128$info, "GSE128048")

## =========================================== cross-dataset module replication
## For each GSE103278 module, score its genes in GSE128048 (mean of z-scored
## expression over the genes shared by both networks) and correlate with the
## drought trait there. Sign concordance of the two correlations = replication.
repl <- do.call(rbind, lapply(names(res103$moduleSizes), function(m) {
  g  <- res103$geneTab$gene[res103$geneTab$module == m]
  g  <- intersect(g, rownames(d128$mat))
  if (length(g) < 10) return(NULL)
  z  <- t(scale(t(d128$mat[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  sc <- colMeans(z)
  tr <- as.integer(d128$info$treatment[match(names(sc), d128$info$sample)] == "drought")
  ct <- cor(sc, tr, use = "p")
  cp <- corPvalueStudent(ct, length(tr))
  r1 <- res103$modTrait[res103$modTrait$module == m & res103$modTrait$trait == "drought", ]
  data.frame(module = m, n_genes_tested = length(g),
             r_GSE103278 = r1$cor, p_GSE103278 = r1$p,
             r_GSE128048 = as.numeric(ct), p_GSE128048 = as.numeric(cp),
             direction_concordant = sign(r1$cor) == sign(ct), row.names = NULL)
}))
repl$fdr_GSE128048 <- p.adjust(repl$p_GSE128048, method = "BH")
logmsg("replication: ", sum(repl$direction_concordant), " / ", nrow(repl),
       " GSE103278 modules keep the drought direction in GSE128048")

## ================================================================== outputs
saveRDS(list(res103 = res103, res128 = res128, repl = repl),
        file.path(OUT, "wgcna_results.rds"))

write.tsv <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t",
                                        quote = FALSE, row.names = FALSE, na = "NA")

write.tsv(res103$geneTab,  "module_genes_GSE103278.tsv")
write.tsv(res128$geneTab,  "module_genes_GSE128048.tsv")
write.tsv(res103$hub,      "hub_genes_GSE103278.tsv")
write.tsv(res128$hub,      "hub_genes_GSE128048.tsv")
write.tsv(res103$modTrait, "module_trait_GSE103278.tsv")
write.tsv(res128$modTrait, "module_trait_GSE128048.tsv")
write.tsv(res103$cross,    "module_consensus_overlap_GSE103278.tsv")
write.tsv(res128$cross,    "module_consensus_overlap_GSE128048.tsv")
write.tsv(res103$sft,      "soft_threshold_fit_GSE103278.tsv")
write.tsv(res128$sft,      "soft_threshold_fit_GSE128048.tsv")
write.tsv(repl,            "module_replication_GSE128048.tsv")

## compact key-numbers summary for the report
summ <- list(
  GSE103278 = list(power = res103$power, power_note = res103$power_note,
                   n_samples = res103$nSamples, n_genes = res103$nGenes,
                   n_modules = length(res103$moduleSizes),
                   module_sizes = res103$moduleSizes,
                   top_module_trait = res103$modTrait[order(res103$modTrait$p)[1:6], ]),
  GSE128048 = list(power = res128$power, power_note = res128$power_note,
                   n_samples = res128$nSamples, n_genes = res128$nGenes,
                   n_modules = length(res128$moduleSizes),
                   module_sizes = res128$moduleSizes,
                   top_module_trait = res128$modTrait[order(res128$modTrait$p)[1:6], ]),
  replication = repl
)
saveRDS(summ, file.path(OUT, "key_numbers.rds"))
logmsg("done; outputs in ", OUT)
