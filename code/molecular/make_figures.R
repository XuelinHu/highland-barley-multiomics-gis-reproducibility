# Regenerate paper figures from the corrected reanalysis (limma v2), base R graphics.
# PCA: prcomp on quantile-normalised log2 intensities. Volcano: mean logFC of
# the four genotype/tissue contrasts vs F-stat FDR. Consensus: cross-dataset
# concordance of direction-consistent probes.
suppressPackageStartupMessages(library(limma))
root <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "")
inp <- file.path(root, 'processed/reanalysis')
out <- file.path(inp, 'figures')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pal <- c(CamB_leaf = '#1f77b4', CamB_root = '#2ca02c', Maresi_leaf = '#ff7f0e', Maresi_root = '#d62728')

datasets <- if ('--consensus-only' %in% commandArgs(trailingOnly = TRUE)) character() else c('GSE103278', 'GSE128048')
for (acc in datasets) {
  obj <- readRDS(file.path(inp, paste0(acc, '.rds')))
  meta <- read.delim(file.path(inp, paste0(acc, '_design.tsv')), check.names = FALSE)
  meta$group <- paste(meta$genotype, meta$tissue, sep = '_')
  pca <- prcomp(t(obj$expression$E), center = TRUE, scale. = TRUE)
  d <- as.data.frame(pca$x); d$group <- meta$group; d$trt <- meta$treatment
  v <- round(100 * summary(pca)$importance[2, 1:2], 1)
  png(file.path(out, paste0(acc, '_PCA.png')), width = 7, height = 5.2, units = 'in', res = 150)
  par(mar = c(4, 4, 3, 6))
  plot(d$PC1, d$PC2, col = pal[d$group], pch = ifelse(d$trt == 'drought', 17, 16),
       cex = 1.6, asp = 1,
       xlab = paste0('PC1 (', v[1], '%)'), ylab = paste0('PC2 (', v[2], '%)'),
       main = paste(acc, 'sample PCA (quantile-normalised log2)'))
  legend('topright', inset = c(-0.28, 0), xpd = TRUE, legend = names(pal),
         col = pal, pch = 16, bty = 'n', cex = 0.85, title = 'group')
  legend('bottomright', inset = c(-0.28, 0), xpd = TRUE, legend = c('control', 'drought'),
         pch = c(16, 17), bty = 'n', cex = 0.85, title = 'treatment')
  dev.off()

  fstat <- topTableF(obj$fit, number = Inf, adjust.method = 'BH', sort.by = 'F')
  mean_fc <- rowMeans(obj$fit$coefficients)
  volc <- data.frame(logFC = mean_fc, FDR = fstat[rownames(obj$fit$coefficients), 'adj.P.Val'])
  sig <- with(volc, FDR < 0.05 & abs(logFC) >= 1)
  png(file.path(out, paste0(acc, '_volcano.png')), width = 6.4, height = 5, units = 'in', res = 150)
  par(mar = c(4, 4, 3, 1))
  plot(volc$logFC, -log10(volc$FDR), pch = 16, cex = 0.35, col = ifelse(sig, '#b22222', 'grey75'),
       xlab = 'mean log2FC across 4 contrasts', ylab = expression(-log[10]~FDR),
       main = paste(acc, 'mean drought effect'),
       sub = sprintf('%d significant probes', sum(sig)))
  abline(h = -log10(0.05), lty = 2); abline(v = c(-1, 1), lty = 2)
  dev.off()
  message(acc, ': figures written, sig=', sum(sig))
}

cons <- read.delim(file.path(inp, 'drought_consensus_probes.tsv'))
top <- read.delim(file.path(inp, 'top100_drought_candidates.tsv'))
stopifnot('ProbeName' %in% names(cons), 'probe' %in% names(top))
is_top <- cons$ProbeName %in% top$probe
stopifnot(sum(is_top) == nrow(top), nrow(top) == 100)
rng <- range(c(cons$gse103_logFC, cons$gse128_logFC))
png(file.path(out, 'consensus_effects.png'), width = 6.4, height = 5.6, units = 'in', res = 300)
par(mar = c(4, 4, 3, 1))
plot(cons$gse103_logFC, cons$gse128_logFC, pch = 16, cex = 0.35, col = 'grey40',
     xlim = rng, ylim = rng, asp = 1,
     xlab = 'GSE103278 mean log2FC', ylab = 'GSE128048 mean log2FC',
     main = 'Direction-consistent drought probes')
abline(0, 1, lty = 2, col = 'grey60')
points(cons$gse103_logFC[is_top], cons$gse128_logFC[is_top],
       pch = 16, cex = 0.8, col = '#b22222')
mtext(sprintf('%d probes, %d top-100 highlighted', nrow(cons), sum(is_top)), cex = 0.8)
dev.off()
message('consensus figure written')
