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
stopifnot(nrow(d) == 14627L, !anyDuplicated(d$cell),
          isTRUE(all.equal(test$p.value, expected$p[expected$method == 'KNN-smoothed'], tolerance = 1e-12)))
prefix <- 'docs/figures/ep300-cell-level-boxplot'
ggsave(paste0(prefix, '.png'), p_box, width = 5, height = 5, dpi = 600, bg = 'white')
grDevices::cairo_pdf(paste0(prefix, '.pdf'), width = 5, height = 5, family = 'sans')
print(p_box); dev.off()
svglite::svglite(paste0(prefix, '.svg'), width = 5, height = 5)
print(p_box); dev.off()
print(test)
