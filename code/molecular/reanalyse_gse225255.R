# Correct GSE225255 reanalysis: sample matching from deposited phenotype metadata.
# Strain (BCS8/BCS24) x treatment (ck/drought), 3 reps; contrasts per strain and mean effect.
suppressPackageStartupMessages(library(limma))
root <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "")
out <- file.path(root, 'processed/reanalysis')
dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(file.path(root, 'processed/metadata', 'GSE225255_phenotype.tsv'),
                   check.names = FALSE, row.names = 1)
meta$strain <- sub('BCS', 'W', meta$characteristics_ch1.0.strain)  # BCS8 -> W8
meta$treatment <- meta$characteristics_ch1.2.treatment
stopifnot(all(meta$treatment %in% c('ck', 'drought')))

path <- file.path(root, 'geo_raw/GSE225255/GSE225255_1_genes_fpkm_expression.txt.gz')
raw <- read.delim(gzfile(path), check.names = FALSE, row.names = 1)
expr <- raw[, paste0('FPKM.', sub('^BCS', 'W', meta$title))]  # BCS8CK1 -> FPKM.W8CK1
colnames(expr) <- rownames(meta)
stopifnot(all(!is.na(expr)))
stopifnot(ncol(expr) == nrow(meta) && nrow(meta) == 12)

# filter very lowly expressed genes before testing
keep <- rowSums(expr >= 1) >= 3
message('genes kept: ', sum(keep), ' / ', nrow(expr))

group <- factor(paste(meta$strain, meta$treatment, sep = '_'))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrasts <- makeContrasts(
  W8_drought_vs_ck  = W8_drought - W8_ck,
  W24_drought_vs_ck = W24_drought - W24_ck,
  mean_drought_vs_ck = (W8_drought + W24_drought - W8_ck - W24_ck) / 2,
  levels = design)
fit <- eBayes(contrasts.fit(lmFit(log2(expr[keep, ] + 1), design), contrasts), trend = TRUE)
for (coef in colnames(contrasts)) {
  tab <- topTable(fit, coef = coef, number = Inf, adjust.method = 'BH', sort.by = 'P')
  tab <- cbind(gene_id = rownames(tab), tab)
  write.table(tab, file.path(out, paste0('GSE225255_', coef, '_limma.tsv')),
              sep = '\t', quote = FALSE, row.names = FALSE)
}
write.table(meta, file.path(out, 'GSE225255_design.tsv'), sep = '\t', quote = FALSE)
saveRDS(list(expression = expr, design = design, contrasts = contrasts, fit = fit),
        file.path(out, 'GSE225255.rds'))
message('GSE225255: finished (', sum(keep), ' genes)')
