# From the package root: Rscript dev/generate_ep300_cell_comparison.R
# Use the cached scores and the short boxplot example in the README; no model refit.
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_classic(base_size = 11, base_family = 'sans'))
d <- read.csv('docs/figures/smoothing-source-data.csv')
md <- data.frame(sample = d$sample, EP300_activity_A = d$original,
                 EP300_activity_smooth = d$smoothed)
readme <- readLines('README.md')
start <- match('<!-- smoothing-example:start -->', readme)
end <- match('<!-- smoothing-example:end -->', readme)
chunks <- strsplit(paste(readme[(start + 1L):(end - 1L)], collapse = '\n'), '```r\n', fixed = TRUE)[[1]][-1]
block <- strsplit(chunks[3], '\n```', fixed = TRUE)[[1]][1]
block <- sub('md <- obj[[]]\n', '', block, fixed = TRUE)
invisible(eval(parse(text = block)))
expected <- read.csv('docs/figures/ep300-cell-level-ttests.csv')
methods <- c('Before smoothing' = 'Original', 'After smoothing' = 'KNN-smoothed')
score_columns <- c('Before smoothing' = 'original', 'After smoothing' = 'smoothed')
expected_p <- setNames(expected$p[match(methods, expected$method)], names(methods))
actual_p <- vapply(tests[names(methods)], function(result) result$p.value, numeric(1))
stopifnot(nrow(d) == 14627L, !anyDuplicated(d$cell),
          nrow(box_data) == 2L * nrow(d),
          identical(levels(box_data$smoothing), names(methods)),
          isTRUE(all.equal(actual_p, expected_p, tolerance = 1e-12)))

built <- ggplot_build(p_box)
layout <- built$layout$layout
box_stats <- built$data[[1L]]
stopifnot(length(built$layout$panel_scales_y) == 1L,
          built$layout$panel_scales_y[[1L]]$get_transformation()$name == 'identity')
groups <- c('Normal', 'AAA')
quartiles <- c('lower', 'middle', 'upper')
for (panel in names(methods)) {
  group_counts <- table(box_data$sample[box_data$smoothing == panel])
  stopifnot(identical(as.integer(group_counts[groups]), c(6085L, 8542L)))
  panel_id <- layout$PANEL[layout$smoothing == panel]
  panel_stats <- box_stats[box_stats$PANEL == panel_id, ]
  actual_quartiles <- as.matrix(panel_stats[order(panel_stats$x), quartiles])
  expected_quartiles <- t(vapply(groups, function(group) {
    quantile(d[[score_columns[[panel]]]][d$sample == group],
             probs = c(0.25, 0.5, 0.75), names = FALSE)
  }, numeric(3)))
  stopifnot(isTRUE(all.equal(unname(actual_quartiles), unname(expected_quartiles),
                            tolerance = 1e-12)))
  message(sprintf('%s: Welch p = %.15g', panel, actual_p[[panel]]))
  print(data.frame(sample = groups, actual_quartiles, row.names = NULL))
}

prefix <- 'docs/figures/ep300-boxplot-before-after'
width <- 8
height <- 4.6
ggsave(paste0(prefix, '.png'), p_box, width = width, height = height, dpi = 600, bg = 'white')
grDevices::cairo_pdf(paste0(prefix, '.pdf'), width = width, height = height, family = 'sans')
print(p_box); dev.off()
svglite::svglite(paste0(prefix, '.svg'), width = width, height = height)
print(p_box); dev.off()
