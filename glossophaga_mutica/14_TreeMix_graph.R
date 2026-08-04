# TreeMix analysis: model comparison (m = 0-5), tree/residual plots, and
# extraction of node coordinates (drift values) underlying the m = 0 figure.
# Requires plotting_funcs.R, distributed with TreeMix
# (https://bitbucket.org/nygcresearch/treemix/wiki/Home) -- not a CRAN/
# Bioconductor package, must be downloaded separately.

library(RColorBrewer)
library(R.utils)
source("/path/to/your/treemix/plotting_funcs.R")

BASE <- "/path/to/your/directory"
setwd(BASE)

# Variance explained per m value
calc_variance <- function(prefix) {
  cov <- read.table(gzfile(paste0(prefix, ".cov.gz")), as.is = TRUE)
  m_cov <- read.table(gzfile(paste0(prefix, ".modelcov.gz")), as.is = TRUE)
  var_obj <- var(cov[lower.tri(cov)])
  v_res <- var(cov[lower.tri(cov)] - m_cov[lower.tri(m_cov)])
  return(1 - v_res/var_obj)
}

variance_results <- data.frame(
  m = 0:5,
  Variance_Explained = sapply(0:5, function(x) calc_variance(paste0("mutica_tm_m", x)))
)
print(variance_results)

# Plot trees
pdf("Supplementary_Figure_TreeMix_Trees.pdf", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))
for (m in 0:5) {
  prefix <- paste0("mutica_tm_m", m)
  plot_tree(prefix)
  title(paste0("m = ", m, "\nVar: ", round(variance_results$Variance_Explained[m+1], 4)))
}
dev.off()

# Plot residuals
pdf("Supplementary_Figure_TreeMix_Residuals.pdf", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(5, 5, 4, 2))
for (m in 0:5) {
  prefix <- paste0("mutica_tm_m", m)
  plot_resid(prefix, "pop_order.txt")
  title(paste0("Residuals (m = ", m, ")"))
}
dev.off()

# Extract node coordinates underlying the m = 0 tree plot (Supplementary
# Figure 3). plot_tree() returns the exact x/y coordinates it draws; the
# x value is the cumulative drift parameter from the root to each node,
# i.e. the same value represented on the plot's x-axis. Not reported as a
# table in the manuscript, but kept here so exact per-population drift
# values are reproducible/available on request (e.g. for reviewers).
source("/path/to/your/treemix/plotting_funcs.R")
result <- plot_tree("mutica_tm_m0")
print(result$d[, c("x", "y")])
