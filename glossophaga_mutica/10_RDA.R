library(vegan)
library(terra)

# Fixed seed for reproducibility -- anova.cca's permutation test showed
# run-to-run variation in Pr(>F) despite reporting an "entire set" of
# permutations; pin the seed so results are identical across reruns.
set.seed(123)

BASE <- "Your directory"
setwd(BASE)

# Map raw sample/locality names to publication codes
code_map <- c(
  "Concha"    = "ConPa",
  "Ixtapa"    = "IxtPa",
  "Poana"     = "PoaGo",
  "Quilamula" = "QuiGo",
  "Hobonil"   = "HobYu",
  "Nizanda"   = "NizOa",
  "Nayarit"   = "IslNa"
)

# Genetic distances
gen_dist <- as.dist(read.table("mutica_ibs.mdist"))

# Coordinates -- sample identifiers are read directly from this file
# rather than hardcoded, so downstream labels always match the actual row
# order used to build points_spat/env_data below.
coords <- read.table("coordinates.csv", header = TRUE)
sample_codes <- code_map[gsub("\\.bam$", "", coords$sample)]

# Local WorldClim layers
bioclim <- rast(list.files(
  path = file.path(BASE, "climate/wc2.1_2.5m"),
  pattern = "\\.tif$",
  full.names = TRUE
))
points_spat <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
env_data <- terra::extract(x = bioclim, y = points_spat, buffer = 5000, fun = mean, ID = FALSE)
env_data_scaled <- scale(env_data)

# 2. Environmental PCA (retain first 3 axes)

env_pca <- prcomp(env_data_scaled, center = FALSE, scale. = FALSE)
env_pca_scores <- env_pca$x[, 1:3]
variance_explained <- summary(env_pca)$importance[2, 1:3] * 100
cumulative_variance <- sum(variance_explained)

cat("Variance explained by PC1:", round(variance_explained[1], 1), "%\n")
cat("Variance explained by PC2:", round(variance_explained[2], 1), "%\n")
cat("Variance explained by PC3:", round(variance_explained[3], 1), "%\n")
cat("Cumulative variance explained by PC1-3:", round(cumulative_variance, 1), "%\n")

cat("\nPC loadings (PC1):\n")
print(sort(env_pca$rotation[, 1]))


# 3. Preliminary check: is there meaningful environmental heterogeneity
#    among the sampled localities in the first place?
# Rationale: RDA/permutation tests ask whether genetic variation tracks
# environmental variation. That question is only meaningful if the
# sampled sites actually differ climatically -- otherwise a weak or
# non-significant environmental effect in the RDA is uninformative
# (it cannot distinguish "no local adaptation" from "no real climatic
# gradient was available to detect an effect against"). With only 7
# sampling localities (one individual per locality, no replication),
# a formal within-vs-between-group test (e.g. PERMANOVA) is not
# possible; instead we characterize environmental heterogeneity
# descriptively and quantitatively before proceeding to the RDA.


# Range and spread of each raw bioclim variable across sites
env_data_df <- as.data.frame(env_data)
rownames(env_data_df) <- sample_codes

env_summary <- data.frame(
  variable = colnames(env_data_df),
  min      = sapply(env_data_df, min),
  max      = sapply(env_data_df, max),
  range    = sapply(env_data_df, function(x) max(x) - min(x)),
  mean     = sapply(env_data_df, mean),
  sd       = sapply(env_data_df, sd),
  cv_pct   = sapply(env_data_df, function(x) 100 * sd(x) / abs(mean(x)))
)

cat("\n--- Environmental heterogeneity among sampled localities (raw bioclim variables) ---\n")
print(env_summary, row.names = FALSE)
write.csv(env_summary, "env_heterogeneity_summary.csv", row.names = FALSE)

# Overall multivariate environmental dispersion among sites (pairwise
# distances in scaled bioclim space; used both as a single-number summary
# here and, as a full matrix, for reference/reporting below)
env_dist <- dist(env_data_scaled)
env_matrix <- as.matrix(env_dist)
rownames(env_matrix) <- colnames(env_matrix) <- sample_codes

env_dist_values <- as.vector(env_dist)
cat("\n--- Overall pairwise environmental distance (scaled bioclim PC space) ---\n")
cat(sprintf("Mean pairwise environmental distance: %.2f\n", mean(env_dist_values)))
cat(sprintf("Range: %.2f -- %.2f\n", min(env_dist_values), max(env_dist_values)))
cat("\nPairwise environmental distances:\n")
print(round(env_matrix, 2))

# Visual check: environmental PCA biplot of sampling sites -- a quick,
# honest way to *show* (not just assert) that sites occupy distinct
# regions of climate space, for the supplementary material.
pdf("env_pca_sites.pdf", width = 6, height = 5)
plot(env_pca_scores[, 1], env_pca_scores[, 2],
     xlab = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
     pch = 19, col = "steelblue", cex = 1.5,
     main = "Environmental PCA of sampled localities")
text(env_pca_scores[, 1], env_pca_scores[, 2], labels = sample_codes, pos = 3, cex = 0.8)
dev.off()

cat("\nSaved: env_heterogeneity_summary.csv and env_pca_sites.pdf\n")
cat("Inspect these before proceeding to the RDA -- do sites show real climatic separation?\n")

# 4. RDA: partition genetic variance among environmental and spatial
#    components


# Geographic distance - first PCNM axis only
coords_matrix <- coords[, c("longitude", "latitude")]
geo_dist_matrix <- terra::distance(coords_matrix, unit = "km")
geo_dist <- as.dist(geo_dist_matrix)
geo_pcnm <- pcnm(geo_dist)
geo_pcnm_1 <- scores(geo_pcnm)[, 1, drop = FALSE]

# Genetic distances as PCoA for RDA
gen_pcoa <- cmdscale(gen_dist, k = nrow(coords) - 1, eig = TRUE)

# RDA full model: genetics ~ geography + environment
rda_full <- rda(gen_pcoa$points ~ geo_pcnm_1 + env_pca_scores)
cat("\nRDA full model: Genetics ~ Geography + Environment\n")
print(rda_full)
set.seed(123)
anova(rda_full)

# RDA partial model: genetics ~ environment | geography
rda_partial <- rda(gen_pcoa$points ~ env_pca_scores + Condition(geo_pcnm_1))
cat("\nRDA partial model: Genetics ~ Environment | Geography\n")
print(rda_partial)
set.seed(123)
anova(rda_partial)

# Write model summaries to file for reproducibility (not just console output)
sink("rda_model_summary.txt")
cat("RDA full model: Genetics ~ Geography + Environment\n")
print(rda_full)
set.seed(123)
print(anova(rda_full))
cat("\nRDA partial model: Genetics ~ Environment | Geography\n")
print(rda_partial)
set.seed(123)
print(anova(rda_partial))
sink()