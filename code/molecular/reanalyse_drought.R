# Reprocess single-colour Agilent data; no reuse of deprecated results.
suppressPackageStartupMessages(library(limma))
root <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "")
out <- file.path(root, 'processed/reanalysis')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
for (acc in c('GSE103278', 'GSE128048')) {
  meta <- read.delim(file.path(root, 'processed/metadata', paste0(acc, '_phenotype.tsv')), check.names=FALSE)
  src <- meta$source_name_ch1
  meta$genotype <- trimws(sub(',.*', '', src))
  meta$tissue <- ifelse(grepl('roots', src), 'root', 'leaf')
  meta$treatment <- ifelse(grepl('drought', src), 'drought', 'control')
  files <- list.files(file.path(root, 'processed/raw_extracted', acc), full.names=TRUE)
  files <- vapply(meta$geo_accession, function(id) {
    found <- files[startsWith(basename(files), paste0(id, '_'))]
    stopifnot(length(found)==1); found
  }, '')
  raw <- read.maimages(files, source='agilent', green.only=TRUE)
  keep <- raw$genes$ControlType == 0
  stopifnot(any(keep), all(is.finite(raw$E)))
  corrected <- backgroundCorrect(raw, method='normexp', offset=50)
  expr <- normalizeBetweenArrays(corrected, method='quantile')
  expr <- expr[keep, ]
  colnames(expr) <- meta$geo_accession
  group <- factor(paste(meta$genotype, meta$tissue, meta$treatment, sep='_'))
  stopifnot(all(table(group)==3))
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  comparisons <- unlist(lapply(c('CamB','Maresi'), function(g)
    paste0(g, '_', c('leaf','root'), '_drought-', g, '_', c('leaf','root'), '_control')))
  contrasts <- makeContrasts(contrasts=comparisons, levels=design)
  fit <- eBayes(contrasts.fit(lmFit(expr,design),contrasts), trend=TRUE)
  for (i in seq_along(comparisons)) {
    tab <- topTable(fit, coef=i, number=Inf, adjust.method='BH', sort.by='P')
    write.table(tab,file.path(out,paste0(acc,'_',i,'_limma.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  }
  write.table(meta,file.path(out,paste0(acc,'_design.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  saveRDS(list(expression=expr, design=design, contrasts=contrasts, fit=fit), file.path(out,paste0(acc,'.rds')))
  message(acc, ': finished four genotype/tissue contrasts')
}
writeLines(capture.output(sessionInfo()), file.path(out,'R_sessionInfo.txt'))
