#!/usr/bin/env Rscript
# Compare memberships on the fixed 8,000-gene universes of the saved networks.
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(4)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript module_membership_sensitivity.R EVIDENCE_ROOT OUTPUT_DIR")
base <- file.path(args[1], "processed/reanalysis")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saved <- readRDS(file.path(base, "wgcna/wgcna_results.rds"))
consensus <- read.delim(file.path(base, "drought_consensus_genes.tsv"))$gene_id
top100 <- read.delim(file.path(base, "top100_drought_candidates.tsv"))$gene_id
overlaps <- list()
membership <- list()
counter <- 0L

for (key in c("res103", "res128")) {
  original <- saved[[key]]
  dataset <- if (key == "res103") "GSE103278" else "GSE128048"
  expression <- original$datExpr
  genes <- colnames(expression)
  ref <- original$geneTab
  refcolors <- ref$module[match(genes, ref$gene)]
  stopifnot(!anyNA(refcolors), length(genes) == 8000)
  drought <- original$traits$drought
  refstats <- original$modTrait[original$modTrait$trait == "drought", ]
  refstats <- refstats[refstats$fdr < .05 & refstats$module != "grey", ]
  for (power in c(12, 18, 26)) {
    if (power == 26) {
      colors <- refcolors
    } else {
      set.seed(42)
      net <- blockwiseModules(expression, power = power, networkType = "signed",
        TOMType = "signed", minModuleSize = 30, mergeCutHeight = .25,
        deepSplit = 2, pamRespectsDendro = FALSE, reassignThreshold = 1e-6,
        numericLabels = TRUE, maxBlockSize = ncol(expression) + 100,
        verbose = 0, saveTOMs = FALSE)
      colors <- labels2colors(net$colors)
    }
    membership[[paste(dataset, power)]] <- data.frame(dataset, power, gene = genes, module = colors)
    mes <- moduleEigengenes(expression, colors)$eigengenes
    cors <- as.numeric(cor(mes, drought, use = "p"))
    names(cors) <- sub("^ME", "", colnames(mes))
    for (i in seq_len(nrow(refstats))) {
      r <- refstats[i, ]
      a <- genes[refcolors == r$module]
      candidates <- names(cors)[sign(cors) == sign(r$cor) & names(cors) != "grey"]
      if (!length(candidates)) next
      jac <- sapply(candidates, function(m) {
        b <- genes[colors == m]
        length(intersect(a, b)) / length(union(a, b))
      })
      chosen <- candidates[which.max(jac)]
      b <- genes[colors == chosen]
      counter <- counter + 1L
      overlaps[[counter]] <- data.frame(dataset, reference_power = 26, power,
        reference_module = r$module, matched_module = chosen,
        reference_drought_r = r$cor, matched_drought_r = cors[[chosen]],
        reference_size = length(a), matched_size = length(b),
        intersection = length(intersect(a, b)), union = length(union(a, b)),
        jaccard = max(jac), reference_retention = length(intersect(a, b)) / length(a),
        reference_consensus = sum(a %in% consensus),
        retained_consensus = sum(intersect(a, b) %in% consensus),
        reference_top100 = sum(a %in% top100),
        retained_top100 = sum(intersect(a, b) %in% top100))
    }
    cat(dataset, "power", power, "completed\n")
    flush.console()
  }
}
write.table(do.call(rbind, overlaps), file.path(out, "wgcna_module_membership_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, membership), file.path(out, "wgcna_memberships_by_power.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "WGCNA_revision_sessionInfo.txt"))
