#!/usr/bin/env Rscript
# =============================================================================
# FALLBACK co-expression analysis (used only if the WGCNA package cannot be
# installed).  Same data, same trait model, but:
#   * similarity  : 1 - Pearson correlation on the top-variance genes
#   * clustering  : average-linkage hclust
#   * module cut  : dynamicTreeCut::cutreeDynamic (hybrid) when available,
#                   otherwise a fixed-height cut
#   * eigengene   : first principal component of each module (sign-aligned to
#                   the module mean), instead of WGCNA's moduleEigengenes
# Writes files with the "_hclust" suffix so it never overwrites WGCNA output.
# =============================================================================
suppressPackageStartupMessages(library(limma))
options(stringsAsFactors = FALSE)

BASE <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "processed/reanalysis")
OUT  <- file.path(BASE, "wgcna")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_TOP_GENES <- 8000
MIN_MODULE_SIZE <- 30
R2_CUT <- 0.85

has_dtc <- requireNamespace("dynamicTreeCut", quietly = TRUE)
cat("[fallback] dynamicTreeCut available:", has_dtc, "\n")

mp <- read.delim(file.path(BASE, "../mapping/GPL14927_morex_mapping.tsv"))
mp <- mp[mp$status == "unique_gene", c("probe_id", "gene_ids")]
names(mp) <- c("probe", "gene")

cons <- read.delim(file.path(BASE, "drought_consensus_genes.tsv"))
t100 <- read.delim(file.path(BASE, "top100_drought_candidates.tsv"))
cons_genes <- unique(cons$gene_id); top100_genes <- unique(t100$gene_id)

prep_dataset <- function(prefix) {
  obj <- readRDS(file.path(BASE, paste0(prefix, ".rds")))
  el <- obj$expression; E <- el$E; pn <- el$genes$ProbeName
  v <- apply(E, 1, var); ord <- order(pn, -v); keeprows <- ord[!duplicated(pn[ord])]
  E <- E[keeprows, , drop = FALSE]; pn <- pn[keeprows]
  d <- obj$design; grp <- colnames(d)[apply(d, 1, which.max)]
  info <- data.frame(sample = colnames(E), group = grp,
                     genotype = sub("^(CamB|Maresi)_.*$", "\\1", grp),
                     tissue   = sub("^[^_]+_([a-z]+)_.*$", "\\1", grp),
                     treatment= sub("^.*_([a-z]+)$", "\\1", grp),
                     row.names = colnames(E))
  sel <- pn %in% mp$probe; Esub <- E[sel, , drop = FALSE]
  g <- mp$gene[match(pn[sel], mp$probe)]; vs <- apply(Esub, 1, var)
  o <- order(g, -vs); Esub <- Esub[o, , drop = FALSE]; g <- g[o]
  k <- !duplicated(g); Eg <- Esub[k, , drop = FALSE]; rownames(Eg) <- g[k]
  list(mat = Eg, info = info)
}

run_fallback <- function(Eg, info, tag) {
  v <- apply(Eg, 1, var)
  top <- names(sort(v, decreasing = TRUE))[1:min(N_TOP_GENES, length(v))]
  datExpr <- t(Eg[top, , drop = FALSE])

  ## ---- soft-thresholding power (same SFT criterion, done by hand) ----
  sft_tab <- do.call(rbind, lapply(c(1:10, 12, 14, 16, 18, 20), function(p) {
    adj <- abs(cor(datExpr, use = "p"))^p
    k <- rowSums(adj) - 1
    d <- data.frame(Power = p, mean.k. = mean(k), SFT.R.sq = NA_real_, slope = NA_real_)
    # scale-free fit on the connectivity distribution (standard approach)
    bins <- cut(k, breaks = 10, labels = FALSE)
    m <- tapply(k, bins, mean); fr <- as.numeric(table(bins)) / length(k)
    ok <- fr > 0 & m > 0
    d$SFT.R.sq <- if (sum(ok) > 2) summary(lm(log10(fr[ok]) ~ log10(m[ok])))$r.squared else NA
    d$slope <- if (sum(ok) > 2) coef(lm(log10(fr[ok]) ~ log10(m[ok])))[2] else NA
    d
  }))
  ok <- !is.na(sft_tab$SFT.R.sq) & sft_tab$SFT.R.sq >= R2_CUT
  if (any(ok)) { power <- sft_tab$Power[which(ok)[1]]
    note <- sprintf("first power with scale-free R^2 >= %.2f", R2_CUT)
  } else {
    power <- sft_tab$Power[which.max(sft_tab$SFT.R.sq)]
    note <- sprintf("no power reached R^2 >= %.2f; used power maximising R^2 (%.3f)",
                    R2_CUT, max(sft_tab$SFT.R.sq, na.rm = TRUE))
  }
  cat(sprintf("[fallback] %s: power = %d (%s)\n", tag, power, note))

  ## ---- adjacency-based dissimilarity + average-linkage clustering ----
  adj <- (1 + cor(datExpr, use = "p")) / 2          # signed hybrid adjacency
  diss <- 1 - adj
  hc <- hclust(as.dist(diss), method = "average")

  if (has_dtc) {
    colors <- dynamicTreeCut::cutreeDynamic(dendro = hc, distM = as.matrix(diss),
                                            deepSplit = 2, pamRespectsDendro = FALSE,
                                            minClusterSize = MIN_MODULE_SIZE,
                                            method = "hybrid")
    labels <- c("grey", paste0("M", seq_len(max(colors))))
    colors <- labels[colors + 1]
  } else {
    cl <- cutree(hc, h = quantile(hc$height, 0.90))
    ## renumber and drop small clusters
    tab <- table(cl)
    keep <- names(tab)[tab >= MIN_MODULE_SIZE]
    colors <- ifelse(cl %in% keep, paste0("M", cl), "grey")
  }
  mods <- sort(setdiff(unique(colors), "grey"))
  cat(sprintf("[fallback] %s: %d modules; sizes: %s\n", tag, length(mods),
              paste(mods, as.integer(table(colors)[mods]), sep = "=", collapse = ", ")))

  ## ---- module eigengene = PC1 per module, sign-aligned to the module mean ----
  ME <- sapply(mods, function(m) {
    sub <- datExpr[, colors == m, drop = FALSE]
    pc <- prcomp(sub, center = TRUE, scale. = FALSE)$x[, 1]
    if (cor(pc, rowMeans(scale(sub))) < 0) pc <- -pc
    as.numeric(pc)
  })
  rownames(ME) <- rownames(datExpr); colnames(ME) <- mods

  traits <- data.frame(
    drought  = as.integer(info$treatment == "drought"),
    genotype = as.integer(info$genotype  == "Maresi"),
    tissue   = as.integer(info$tissue    == "root"),
    row.names = info$sample)[rownames(datExpr), , drop = FALSE]

  mt <- do.call(rbind, lapply(names(traits), function(tr) {
    ct <- cor(ME, traits[[tr]], use = "p")
    data.frame(module = rownames(ct), trait = tr, cor = as.numeric(ct),
               p = as.numeric(corPvalueStudent(ct, nrow(traits))), row.names = NULL)
  }))
  mt$fdr <- p.adjust(mt$p, method = "BH")

  kME <- cor(datExpr, ME, use = "p")                # module membership
  gene_tab <- data.frame(gene = colnames(datExpr), module = colors, row.names = NULL)
  gene_tab$kME_module <- sapply(seq_len(nrow(gene_tab)), function(i) {
    m <- gene_tab$module[i]; if (m == "grey") NA_real_ else kME[i, m]
  })
  gene_tab$variance <- apply(datExpr, 2, var)
  gene_tab$GS_drought  <- cor(datExpr, traits$drought, use = "p")
  gene_tab$GS_genotype <- cor(datExpr, traits$genotype, use = "p")
  gene_tab$GS_tissue   <- cor(datExpr, traits$tissue, use = "p")
  gene_tab$is_consensus <- gene_tab$gene %in% cons_genes
  gene_tab$is_top100    <- gene_tab$gene %in% top100_genes

  hub <- do.call(rbind, lapply(mods, function(m) {
    s <- gene_tab[gene_tab$module == m, ]; s <- s[order(-s$kME_module), ]; head(s, 20)
  }))

  N <- nrow(gene_tab)
  cross <- do.call(rbind, lapply(mods, function(m) {
    inmod <- gene_tab$gene[gene_tab$module == m]; n <- length(inmod)
    do <- function(set, label) {
      K <- sum(gene_tab$gene %in% set); k <- sum(inmod %in% set)
      data.frame(module = m, module_size = n, gene_set = label, set_in_network = K,
                 overlap = k, expected = n * K / N,
                 fold_enrichment = ifelse(n * K > 0, k / (n * K / N), NA_real_),
                 p_hypergeom = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
                 row.names = NULL)
    }
    rbind(do(cons_genes, "drought_consensus_1197"), do(top100_genes, "top100_drought_candidates"))
  }))
  cross$fdr <- p.adjust(cross$p_hypergeom, method = "BH")

  ## ---- plots ----
  pretty <- c(drought = "Drought vs Control", genotype = "Genotype (Maresi vs CamB)",
              tissue = "Tissue (Root vs Leaf)")
  M <- matrix(NA_real_, nrow = length(mods), ncol = ncol(traits),
              dimnames = list(mods, pretty[colnames(traits)])); P <- M
  for (i in seq_len(nrow(mt))) { M[mt$module[i], pretty[mt$trait[i]]] <- mt$cor[i]
                                 P[mt$module[i], pretty[mt$trait[i]]] <- mt$p[i] }
  png(file.path(OUT, paste0("module_trait_heatmap_", tag, "_hclust.png")),
      width = 1500, height = max(500, 110 * nrow(M) + 300), res = 130)
  par(mar = c(7, 8, 4, 2))
  image(seq_len(ncol(M)), seq_len(nrow(M)), t(M[nrow(M):1, , drop = FALSE]),
        col = colorRampPalette(c("blue", "white", "red"))(50), zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "",
        main = paste0("Module-trait relationships (", tag, ", hclust fallback)\n",
                      "n = ", nrow(datExpr), " samples"))
  axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(nrow(M)), labels = rev(rownames(M)), las = 1, cex.axis = 0.8)
  for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M)))
    text(j, nrow(M) - i + 1, sprintf("%.2f\n(%.1e)", M[i, j], P[i, j]), cex = 0.6)
  box(); dev.off()

  png(file.path(OUT, paste0("module_dendrogram_", tag, "_hclust.png")), width = 1500, height = 800, res = 130)
  plot(hc, labels = FALSE, hang = 0.03, main = paste0(tag, " hclust dendrogram (power = ", power, ")"))
  dev.off()

  list(power = power, power_note = note, moduleSizes = setNames(as.integer(table(colors)[mods]), mods),
       modTrait = mt, geneTab = gene_tab, hub = hub, cross = cross, colors = colors,
       hc = hc, MEs = ME, sft = sft_tab, nSamples = nrow(datExpr), nGenes = ncol(datExpr))
}

d103 <- prep_dataset("GSE103278"); res103 <- run_fallback(d103$mat, d103$info, "GSE103278")
d128 <- prep_dataset("GSE128048"); res128 <- run_fallback(d128$mat, d128$info, "GSE128048")

repl <- do.call(rbind, lapply(names(res103$moduleSizes), function(m) {
  g <- res103$geneTab$gene[res103$geneTab$module == m]
  g <- intersect(g, rownames(d128$mat)); if (length(g) < 10) return(NULL)
  z <- t(scale(t(d128$mat[g, , drop = FALSE]))); z[is.na(z)] <- 0
  sc <- colMeans(z)
  tr <- as.integer(d128$info$treatment[match(names(sc), d128$info$sample)] == "drought")
  ct <- cor(sc, tr, use = "p")
  r1 <- res103$modTrait[res103$modTrait$module == m & res103$modTrait$trait == "drought", ]
  data.frame(module = m, n_genes_tested = length(g), r_GSE103278 = r1$cor, p_GSE103278 = r1$p,
             r_GSE128048 = as.numeric(ct),
             p_GSE128048 = as.numeric(corPvalueStudent(ct, length(tr))),
             direction_concordant = sign(r1$cor) == sign(ct), row.names = NULL)
}))

saveRDS(list(res103 = res103, res128 = res128, repl = repl, method = "hclust"),
        file.path(OUT, "hclust_results.rds"))
w <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "NA")
w(res103$geneTab, "module_genes_GSE103278_hclust.tsv")
w(res128$geneTab, "module_genes_GSE128048_hclust.tsv")
w(res103$hub, "hub_genes_GSE103278_hclust.tsv")
w(res128$hub, "hub_genes_GSE128048_hclust.tsv")
w(res103$modTrait, "module_trait_GSE103278_hclust.tsv")
w(res128$modTrait, "module_trait_GSE128048_hclust.tsv")
w(res103$cross, "module_consensus_overlap_GSE103278_hclust.tsv")
w(res128$cross, "module_consensus_overlap_GSE128048_hclust.tsv")
w(res103$sft, "soft_threshold_fit_GSE103278_hclust.tsv")
w(repl, "module_replication_GSE128048_hclust.tsv")
cat("[fallback] done\n")
