# Run from the package root: Rscript dev/generate_ep300_cell_comparison.R
# Exploratory change of analysis unit: cells, using previously inferred EP300 scores.
# Quantitative grid; R-only; original and smoothed panels; all 14627 cells retained.
Sys.setenv(LC_ALL = 'C')
suppressPackageStartupMessages(library(ggplot2))
input <- 'docs/figures/smoothing-source-data.csv'
out <- 'docs/figures'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(input)
stopifnot(nrow(d) == 14627L, !anyDuplicated(d$cell),
  all(is.finite(d$original)), all(is.finite(d$smoothed)),
  all(d$original >= 0), all(d$smoothed >= 0),
  all(d$sample %in% c('Normal', 'AAA')))
d$condition <- factor(d$sample, levels = c('Normal', 'AAA'))
box_data <- rbind(
  data.frame(cell = d$cell, orig.ident = d$orig.ident, condition = d$condition,
    method = 'Original', activity = d$original),
  data.frame(cell = d$cell, orig.ident = d$orig.ident, condition = d$condition,
    method = 'KNN-smoothed', activity = d$smoothed))
box_data$method <- factor(box_data$method, levels = c('Original', 'KNN-smoothed'))
results <- do.call(rbind, lapply(levels(box_data$method), function(m) {
  x <- box_data[box_data$method == m, ]
  normal <- x$activity[x$condition == 'Normal']
  aaa <- x$activity[x$condition == 'AAA']
  test <- t.test(aaa, normal, alternative = 'two.sided', paired = FALSE, var.equal = FALSE)
  data.frame(method = m, analysis_unit = 'cell', n_Normal = length(normal), n_AAA = length(aaa),
    mean_Normal = mean(normal), mean_AAA = mean(aaa),
    difference_AAA_minus_Normal = mean(aaa) - mean(normal),
    t = unname(test$statistic), df = unname(test$parameter), p = test$p.value,
    ci_low = test$conf.int[1], ci_high = test$conf.int[2])
}))
stopifnot(all(results$n_Normal == 6085L), all(results$n_AAA == 8542L))
sample_results <- read.csv('docs/figures/smoothing-ttests.csv')
results$sample_level_p <- sample_results$p[match(results$method, sample_results$method)]
results$method <- factor(results$method, levels = levels(box_data$method))
results$label <- paste0('Cell-level t-test: p = ', formatC(results$p, digits = 2, format = 'e'))
results$bracket_y <- max(box_data$activity) * 1.06

p <- ggplot(box_data, aes(condition, activity, fill = condition)) +
  geom_boxplot(width = 0.5, outlier.size = 0.25, outlier.alpha = 0.15, linewidth = 0.35) +
  geom_segment(data = results, inherit.aes = FALSE,
    aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y), linewidth = 0.35) +
  geom_segment(data = results, inherit.aes = FALSE,
    aes(x = 1, xend = 1, y = bracket_y, yend = bracket_y * 0.96), linewidth = 0.35) +
  geom_segment(data = results, inherit.aes = FALSE,
    aes(x = 2, xend = 2, y = bracket_y, yend = bracket_y * 0.96), linewidth = 0.35) +
  geom_text(data = results, inherit.aes = FALSE,
    aes(x = 1.5, y = bracket_y * 1.09, label = label), size = 3.4) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(values = c(Normal = '#4D8299', AAA = '#C47F62')) +
  scale_y_continuous(limits = c(0, max(results$bracket_y) * 1.2),
    expand = expansion(mult = c(0, 0.02))) +
  coord_transform(y = 'sqrt') +
  labs(title = 'EP300 activity: cell-level comparison', x = NULL, y = 'Inferred activity') +
  theme_classic(base_size = 11, base_family = 'sans') +
  theme(legend.position = 'none', strip.background = element_blank(),
    plot.title = element_text(face = 'bold', hjust = 0.5))

# The test uses untransformed activity; the square root applies only to display.
built <- ggplot_build(p)
for (panel in 1:2) for (group in 1:2) {
  x <- box_data$activity[as.integer(box_data$method) == panel &
    as.integer(box_data$condition) == group]
  row <- built$data[[1]][as.integer(built$data[[1]]$PANEL) == panel &
    built$data[[1]]$group == group, ]
  stopifnot(nrow(row) == 1L,
    isTRUE(all.equal(as.numeric(row[1, c('lower', 'middle', 'upper')]),
      as.numeric(quantile(x, c(.25, .5, .75))), tolerance = 1e-12)))
}
stopifnot(nrow(box_data) == 2L * nrow(d), !anyNA(results$p),
  all(results$p > 0 & results$p <= 1), is.null(p$labels$subtitle), is.null(p$labels$caption))
write.csv(results, file.path(out, 'ep300-cell-level-ttests.csv'), row.names = FALSE)
# Cell-level source values are already supplied in smoothing-source-data.csv.
ggsave(file.path(out, 'ep300-cell-level-boxplot.png'), p, width = 8, height = 5,
  units = 'in', dpi = 600, bg = 'white')
grDevices::cairo_pdf(file.path(out, 'ep300-cell-level-boxplot.pdf'), width = 8, height = 5, family = 'sans')
print(p); grDevices::dev.off()
svglite::svglite(file.path(out, 'ep300-cell-level-boxplot.svg'), width = 8, height = 5)
print(p); grDevices::dev.off()
writeLines(c('EP300 exploratory cell-level comparison; existing model fits and scores unchanged.',
  'All cells included: Normal=6085, AAA=8542; no subsampling or score clipping.',
  'Two-sided unpaired Welch t-tests use untransformed cell activity values.',
  'Two nominal p-values without multiple-testing correction; full statistics and 95% CIs exported.',
  'The intervals and p-values assume independent cells and ignore sample clustering.',
  'Cells from the same biological sample are not independent biological replicates.',
  'Pooled KNN smoothing adds dependence. Treat both cell-level tests as exploratory.',
  'Box median and IQR; whiskers 1.5 IQR; all outlying cells displayed.',
  'Square-root display transform applied after box statistics; all four quartile summaries verified.',
  'Centered title; no subtitle or bottom caption. R-only, 600-dpi PNG and editable PDF/SVG.',
  paste('Source:', input), 'Original sample-level test results remain available and unchanged.'),
  file.path(out, 'ep300-cell-level-QA.txt'))
writeLines(sub('[[:blank:]]+$', '', capture.output(sessionInfo())),
  file.path(out, 'ep300-cell-level-sessionInfo.txt'))
print(results[, setdiff(names(results), c('label', 'bracket_y'))], digits = 8)
