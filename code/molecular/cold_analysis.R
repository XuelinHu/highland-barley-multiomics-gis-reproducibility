# Low-temperature validation: within-platform limma analyses of the four cold
# series, from the official series matrices (downloaded 2026-09-11).
#   GSE41513 (GPL16130, 2-colour): log2 ratios vs 0H -> one-sample tests per group.
#   GSE27821/GSE27822 (GPL1340): linear Affy -> log2, treatment contrasts.
#   GSE47882 (GPL17286): linear custom Agilent -> log2, barley only, day contrasts.
suppressPackageStartupMessages(library(limma))
root <- file.path(Sys.getenv("BARLEY_DATA_ROOT", "evidence"), "")
raw <- file.path(root, 'geo_raw')
out <- file.path(root, 'processed/reanalysis/cold')
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path))
  trow <- sub('^!Sample_title\t', '', lines[startsWith(lines, '!Sample_title')])
  titles <- gsub('"', '', unlist(strsplit(trow, '\t')))
  titles <- titles[nzchar(titles)]
  start <- which(lines == '!series_matrix_table_begin') + 1
  mat <- read.delim(gzfile(path), skip = start - 1, header = TRUE, quote = '"',
                    check.names = FALSE, comment.char = '!')
  ids <- mat[[1]]; mat <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(mat) <- 'numeric'; rownames(mat) <- ids
  stopifnot(ncol(mat) == length(titles))
  list(mat = mat, titles = titles)
}

write_tables <- function(fit, prefix, coefs) {
  for (cn in coefs) {
    tab <- topTable(fit, coef = cn, number = Inf, adjust.method = 'BH', sort.by = 'P')
    write.table(tab, file.path(out, paste0(prefix, '_', cn, '_limma.tsv')),
                sep = '\t', quote = FALSE, row.names = TRUE)
  }
}

## ---- GSE41513: one-sample tests of log2 ratios vs 0H -----------------------
sm <- read_series_matrix(file.path(raw, 'GSE41513/GSE41513_series_matrix.txt.gz'))
tissue <- sub('Cold (shoot|root) time course.*', '\\1', sm$titles)
time <- sub('.*time course (\\d+)H/0H.*', '\\1H', sm$titles)
grp <- factor(paste(tissue, time, sep = '_'))
tab <- table(grp)
test_grps <- names(tab)[tab >= 2]          # 24H groups have a single replicate
keep <- grp %in% test_grps
g2 <- droplevels(grp[keep])
design <- model.matrix(~ 0 + g2); colnames(design) <- levels(g2)
E <- sm$mat[, keep]
fit <- eBayes(lmFit(E, design))
write_tables(fit, 'GSE41513', colnames(design))
write.table(data.frame(title = sm$titles, group = as.character(grp)),
            file.path(out, 'GSE41513_groups.tsv'), sep = '\t', quote = FALSE)
message('GSE41513 done: groups tested = ', paste(colnames(design), collapse = ', '))

## ---- GSE27821: vernalization/cold contrasts (n = 3 per group) ---------------
sm <- read_series_matrix(file.path(raw, 'GSE27821/GSE27821_series_matrix.txt.gz'))
trt <- sub('treatment: (.*) - rep\\d+', '\\1', sm$titles)
grp <- factor(make.names(gsub('[- ]', '', trt)))
E <- log2(sm$mat + 1)
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
contrasts <- makeContrasts(Shorttermcold_vs_Nonvernalized = Shorttermcold - Nonvernalized,
                           Vernalized_vs_Nonvernalized = Vernalized - Nonvernalized,
                           X1daypostvernalized_vs_Vernalized = X1daypostvernalized - Vernalized,
                           levels = design)
fit <- eBayes(contrasts.fit(lmFit(E, design), contrasts))
write_tables(fit, 'GSE27821', colnames(contrasts))
message('GSE27821 done')

## ---- GSE27822: post-vernalized vs non-vernalized (n = 3 per group) ----------
sm <- read_series_matrix(file.path(raw, 'GSE27822/GSE27822_series_matrix.txt.gz'))
trt <- sub('treatment: (.*) - rep\\d+.*', '\\1', sm$titles)
grp <- factor(gsub('[- ]', '', trt))
E <- log2(sm$mat + 1)
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
contrasts <- makeContrasts(Postvernalized_vs_Nonvernalized = Postvernalized - Nonvernalized,
                           levels = design)
fit <- eBayes(contrasts.fit(lmFit(E, design), contrasts))
write_tables(fit, 'GSE27822', colnames(contrasts))
message('GSE27822 done')

## ---- GSE47882: barley subset, day contrasts blocked by genotype -------------
sm <- read_series_matrix(file.path(raw, 'GSE47882/GSE47882_series_matrix.txt.gz'))
sp <- sub('^(\\w+)_genotype:.*', '\\1', sm$titles)
gen <- sub('^\\w+_genotype:(\\S+)_day\\d+.*', '\\1', sm$titles)
day <- factor(sub('.*_day(\\d+)_.*', 'd\\1', sm$titles))
keep <- sp == 'barley'
E <- log2(sm$mat[, keep] + 1)
design <- model.matrix(~ 0 + day[keep] + gen[keep])
colnames(design) <- c(levels(day), paste0('geno_', levels(factor(gen[keep]))[-1]))
contrasts <- makeContrasts(d1_vs_d0 = d1 - d0, d7_vs_d0 = d7 - d0, d21_vs_d0 = d21 - d0,
                           levels = design)
fit <- eBayes(contrasts.fit(lmFit(E, design), contrasts))
write_tables(fit, 'GSE47882', colnames(contrasts))
write.table(data.frame(title = sm$titles[keep], genotype = gen[keep], day = as.character(day[keep])),
            file.path(out, 'GSE47882_barley_groups.tsv'), sep = '\t', quote = FALSE)
message('GSE47882 done')
writeLines(capture.output(sessionInfo()), file.path(out, 'R_sessionInfo.txt'))
