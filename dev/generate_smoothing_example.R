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
expected_p <- expected$p[match(c('Original', 'KNN-smoothed'), expected$method)]
observed_p <- unname(vapply(tests, function(test) test$p.value, numeric(1)))
stopifnot(isTRUE(all.equal(observed_p, expected_p, tolerance = 1e-12)))
feature_build <- ggplot_build(p_feature)
box_build <- ggplot_build(p_box)
stopifnot(all(table(feature_data$smoothing) == ncol(obj)),
  length(unique(feature_build$layout$layout$SCALE_X)) == 1L,
  length(unique(feature_build$layout$layout$SCALE_Y)) == 1L,
  length(unique(box_build$layout$layout$SCALE_Y)) == 1L)
box_stats <- box_build$data[[1]]
for (panel in seq_along(panel_order)) {
  for (group in 1:2) {
    values <- box_data$Activity[box_data$smoothing == panel_order[panel] &
      box_data$sample == levels(box_data$sample)[group]]
    row <- box_stats[as.integer(box_stats$PANEL) == panel & box_stats$group == group, ]
    stopifnot(nrow(row) == 1L,
      isTRUE(all.equal(as.numeric(row[1, c('lower', 'middle', 'upper')]),
                       as.numeric(quantile(values, c(.25, .5, .75))), tolerance = 1e-12)))
  }
}

export <- function(p, name, width, height) {
  prefix <- file.path('docs/figures', name)
  ggsave(paste0(prefix, '.png'), p, width = width, height = height, dpi = 600, bg = 'white')
  grDevices::cairo_pdf(paste0(prefix, '.pdf'), width = width, height = height, family = 'sans')
  print(p); dev.off()
  svglite::svglite(paste0(prefix, '.svg'), width = width, height = height)
  print(p); dev.off()
}
export(p_feature, 'smoothing-featureplot', 8, 4.5)
export(p_box, 'ep300-boxplot-before-after', 8, 4.6)
writeLines(c('The README examples use all 14627 macrophages without sampling or score clipping.',
  'Original and smoothed scores and both cell-level p-values match the saved source data.',
  'Boxplot quartiles verified on the untransformed scores. Panels share a linear y-axis.',
  'UMAP panels retain the same coordinates and share a linear colour scale.',
  'Cell-level Welch tests are exploratory because cells within samples and smoothed scores are dependent.',
  'R exports include 600-dpi PNG and editable PDF/SVG.'),
  'docs/figures/smoothing-QA.txt')
print(tests)
