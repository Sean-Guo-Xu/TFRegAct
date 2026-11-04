# From the package root: Rscript dev/generate_smoothing_example.R <Seurat.RData> <activity.csv> EP300
# R quantitative panels: all scored macrophages, optional KNN smoothing, exploratory cell-level tests.
suppressPackageStartupMessages({library(TFRegAct); library(Seurat); library(ggplot2); library(patchwork)})
theme_set(theme_classic(base_size = 11, base_family = 'sans'))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[3] == 'EP300')
a <- read.csv(args[2])
e <- new.env(); load(args[1], envir = e)
obj <- subset(e$pbmc, cells = a$cell)
rm(e); gc()
j <- match(colnames(obj), a$cell)
stopifnot(!anyNA(j), !anyDuplicated(a$cell), ncol(obj) == 14627L,
          all(obj[[]]$cell == 'Macrophage'), identical(as.character(obj$sample), a$condition_label[j]))
obj$EP300_activity_A <- a$activity_mean[j]
activity_fit <- list(seurat_object = obj)

# Run the three short examples shown in the README.
readme <- readLines('README.md')
start <- match('<!-- smoothing-example:start -->', readme)
end <- match('<!-- smoothing-example:end -->', readme)
chunks <- strsplit(paste(readme[(start + 1L):(end - 1L)], collapse = '\n'), '```r\n', fixed = TRUE)[[1]][-1]
blocks <- lapply(chunks, function(x) strsplit(x, '\n```', fixed = TRUE)[[1]][1])
stopifnot(length(blocks) == 3L)
for (block in blocks) suppressMessages(eval(parse(text = block)))

# Verify that simplifying the figures has not changed the scores or tests.
source_data <- read.csv('docs/figures/smoothing-source-data.csv')
source_index <- match(colnames(obj), source_data$cell)
stopifnot(!anyNA(source_index), all(is.finite(obj$EP300_activity_smooth)),
  isTRUE(all.equal(obj$EP300_activity_smooth, source_data$smoothed[source_index], check.attributes = FALSE)),
  isTRUE(all.equal(obj$EP300_activity_A, source_data$original[source_index], check.attributes = FALSE)))
expected <- read.csv('docs/figures/ep300-cell-level-ttests.csv')
stopifnot(isTRUE(all.equal(test$p.value, expected$p[expected$method == 'KNN-smoothed'], tolerance = 1e-12)))
box_stats <- ggplot_build(p_box)$data[[1]]
for (group in 1:2) {
  values <- md$EP300_activity_smooth[md$sample == c('Normal', 'AAA')[group]]
  row <- box_stats[box_stats$group == group, ]
  stopifnot(isTRUE(all.equal(as.numeric(row[1, c('lower', 'middle', 'upper')]),
                            as.numeric(quantile(values, c(.25, .5, .75))), tolerance = 1e-12)))
}

export <- function(p, name, width, height) {
  prefix <- file.path('docs/figures', name)
  ggsave(paste0(prefix, '.png'), p, width = width, height = height, dpi = 600, bg = 'white')
  grDevices::cairo_pdf(paste0(prefix, '.pdf'), width = width, height = height, family = 'sans')
  print(p); dev.off()
  svglite::svglite(paste0(prefix, '.svg'), width = width, height = height)
  print(p); dev.off()
}
export(p_feature, 'smoothing-featureplot', 6, 5.4)
export(p_box, 'ep300-cell-level-boxplot', 5, 5)
writeLines(c('Exact README examples executed on all 14627 macrophages; no sampling or score clipping.',
  'Smoothing scores and cell-level p-values match the previous saved results.',
  'Both box quartiles verified on the untransformed scores.',
  'FeaturePlot uses a linear colour scale; original UMAP coordinates retained.',
  'Cell-level Welch tests are exploratory; cells within samples and smoothed scores are dependent.',
  'R only; 600-dpi PNG and editable PDF/SVG; centered titles and no captions.'),
  'docs/figures/smoothing-QA.txt')
print(test)
