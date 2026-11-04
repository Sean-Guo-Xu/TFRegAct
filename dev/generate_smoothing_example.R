# Run from the package root with: Rscript dev/generate_smoothing_example.R
#   <Seurat.RData> <exported_activity.csv> <target_tf>
# Figure contract: illustrate score post-processing using real inferred activity.
# Quantitative grid, R only, all supplied cells; exploratory sample-level Welch t-tests.
suppressPackageStartupMessages({library(TFRegAct); library(Seurat); library(ggplot2); library(patchwork)})
theme_set(theme_classic(base_size = 11, base_family = 'sans'))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, file.exists(args[1]), file.exists(args[2]))
e <- new.env(); load(args[1], envir = e)
stopifnot(exists('pbmc', envir = e, inherits = FALSE))
a <- read.csv(args[2])
stopifnot(!anyDuplicated(a$cell), all(is.finite(a$activity_mean)), all(a$activity_mean >= 0))
obj <- subset(e$pbmc, cells = a$cell)
rm(e); gc()
cells <- colnames(obj)
j <- match(cells, a$cell)
stopifnot(ncol(obj) == nrow(a), !anyNA(j),
          all(obj[[]]$cell == 'Macrophage'),
          identical(as.character(obj$sample), as.character(a$condition_label[j])))
obj[[paste0(args[3], '_activity_A')]] <- stats::setNames(a$activity_mean[j], cells)
activity_fit <- list(target_tf = args[3], seurat_object = obj)
rm(obj, a); gc()

# Execute the exact three example blocks embedded in the README.
readme <- readLines('README.md')
first <- which(readme == '<!-- smoothing-example:start -->')
last <- which(readme == '<!-- smoothing-example:end -->')
stopifnot(length(first) == 1L, length(last) == 1L, first < last)
example <- readme[seq.int(first + 1L, last - 1L)]
inside <- FALSE; blocks <- list(); block <- character()
for (line in example) {
  if (identical(line, '```r')) {
    stopifnot(!inside); inside <- TRUE; block <- character()
  } else if (identical(line, '```') && inside) {
    inside <- FALSE; blocks[[length(blocks) + 1L]] <- block
  } else if (inside) block <- c(block, line)
}
stopifnot(!inside, length(blocks) == 3L)
for (block in blocks) eval(parse(text = block), envir = .GlobalEnv)
stopifnot(length(smooth_score) == length(raw_score),
          identical(names(smooth_score), names(raw_score)),
          all(is.finite(smooth_score)), all(smooth_score >= 0),
          nrow(box_data) == 2L * length(cells), !anyNA(box_data$condition))
stopifnot(nrow(sample_means) == 2L * length(unique(activity_seurat$orig.ident)),
          all(is.finite(ttest_results$p)), all(ttest_results$p >= 0 & ttest_results$p <= 1),
          all(ttest_results$n_Normal == 3L), all(ttest_results$n_AAA == 7L))
stopifnot(is.null(p_feature$patches$annotation$subtitle),
          is.null(p_feature$patches$annotation$caption),
          is.null(p_box$labels$subtitle), is.null(p_box$labels$caption),
          identical(p_box$theme$plot.title$hjust, 0.5))

# Check every plotted cell, raw score and common colour/coordinate scale.
for (i in seq_along(feature_panels)) {
  expected_column <- c(raw_column, smooth_column)[i]
  stopifnot(setequal(rownames(feature_panels[[i]]$data), cells))
  expected <- activity_seurat[[expected_column]][rownames(feature_panels[[i]]$data), 1]
  stopifnot(isTRUE(all.equal(feature_panels[[i]]$data[[expected_column]], expected, check.attributes = FALSE)))
}
builds <- lapply(feature_panels, ggplot_build)
ranges <- lapply(builds, function(b) b$layout$panel_params[[1]][c('x.range', 'y.range')])
stopifnot(isTRUE(all.equal(ranges[[1]], ranges[[2]])))
box_statistics <- ggplot_build(p_box)$data[[1]]
for (panel in 1:2) for (group in 1:2) {
  values <- box_data$activity[box_data$method == levels(box_data$method)[panel] &
                               box_data$condition == levels(box_data$condition)[group]]
  row <- box_statistics[as.integer(box_statistics$PANEL) == panel & box_statistics$group == group, ]
  stopifnot(nrow(row) == 1L,
            isTRUE(all.equal(as.numeric(row[1, c('lower', 'middle', 'upper')]),
                             as.numeric(quantile(values, c(.25, .5, .75))), tolerance = 1e-12)))
}

out <- file.path('docs', 'figures'); dir.create(out, recursive = TRUE, showWarnings = FALSE)
embedding <- Embeddings(activity_seurat, 'umap')[cells, 1:2, drop = FALSE]
source_data <- data.frame(cell = cells, sample = as.character(activity_seurat$sample),
  orig.ident = activity_seurat$orig.ident, batch = activity_seurat$batch,
  UMI = as.numeric(umi[cells]), UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2],
  original = as.numeric(raw_score[cells]), smoothed = as.numeric(smooth_score[cells]))
write.csv(source_data, file.path(out, 'smoothing-source-data.csv'), row.names = FALSE)
summary_table <- do.call(rbind, lapply(split(source_data, source_data$sample), function(d)
  data.frame(condition = d$sample[1], cells = nrow(d),
    original_mean = mean(d$original), smoothed_mean = mean(d$smoothed),
    original_median = median(d$original), smoothed_median = median(d$smoothed),
    original_sd = sd(d$original), smoothed_sd = sd(d$smoothed))))
write.csv(summary_table, file.path(out, 'smoothing-summary.csv'), row.names = FALSE)
write.csv(sample_means, file.path(out, 'smoothing-sample-means.csv'), row.names = FALSE)
write.csv(ttest_results, file.path(out, 'smoothing-ttests.csv'), row.names = FALSE)

export <- function(plot, filename, width, height) {
  ggsave(file.path(out, paste0(filename, '.png')), plot, width = width, height = height,
         units = 'in', dpi = 600, bg = 'white')
  grDevices::cairo_pdf(file.path(out, paste0(filename, '.pdf')), width = width, height = height, family = 'sans')
  print(plot); grDevices::dev.off()
  svglite::svglite(file.path(out, paste0(filename, '.svg')), width = width, height = height)
  print(plot); grDevices::dev.off()
}
export(p_feature, 'smoothing-featureplot', 10.4, 5.4)
export(p_box, 'smoothing-boxplot', 8, 5)
writeLines(c('Purpose: demonstrate optional UMI-adaptive KNN post-processing with exploratory sample-level tests.',
  paste('TF:', tf, '| Cells:', length(cells), '| Coordinate dimensions:', ncol(coordinate), '| k: 15'),
  'All scored macrophages retained; no subsampling or value clipping.',
  'KNN search uses integrated.cca; existing UMAP coordinates are used only for display.',
  'Neighbours can cross samples/batches; no grouping constraint is applied.',
  'Raw scores retained; no model refitting and no uncertainty recalculation.',
  'FeaturePlot panels verified for cell identities, exact score values and common UMAP ranges.',
  'FeaturePlot has a shared square-root colour scale, without changing scores.',
  'Boxplot statistics use raw scores; square-root coordinate transform is applied afterwards.',
  'All four box median/quartile summaries verified against the raw numeric values.',
  'All cells contribute to boxes and outlier points; diamond markers show sample means.',
  'Two-sided unpaired Welch t-tests use orig.ident mean activity: 3 Normal and 7 AAA sample IDs.',
  'Nominal p-values are displayed without multiple-testing correction; both tests are exploratory.',
  'Smoothed scores borrow across samples; smoothed sample means are not strictly independent.',
  'Full t statistics, Welch degrees of freedom, mean differences and 95% confidence intervals exported.',
  'Both figure titles centered; subtitles and bottom captions removed.',
  'R-only rendering; 600-dpi PNG and editable vector PDF/SVG.',
  'README code blocks were executed directly to generate these panels.'),
  file.path(out, 'smoothing-QA.txt'))
writeLines(sub('[[:blank:]]+$', '', capture.output(sessionInfo())),
  file.path(out, 'smoothing-sessionInfo.txt'))
print(summary_table)
print(ttest_results)
cat('README smoothing demonstration completed.\n')
