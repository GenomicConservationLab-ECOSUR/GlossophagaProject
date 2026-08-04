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

# Load genetic distances
gen_dist <- as.dist(read.table("mutica_ibs.mdist")) # From the 08_RDA_prep.sh script.

# Load coordinates -- sample identifiers are read directly from this file
# rather than hardcoded, so downstream labels always match the actual row
# order used to build points_spat/env_data below.
coords <- read.table("coordinates.csv", header = TRUE) # You can obtain the coordinates for this study in Sup. table 1
sample_codes <- code_map[gsub("\\.bam$", "", coords$sample)]

# Load local WorldClim layers
bioclim <- rast(list.files(
  path = file.path(BASE, "climate/wc2.1_2.5m"),
  pattern = "\\.tif$",
  full.names = TRUE
))
points_spat <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
env_data <- terra::extract(x = bioclim, y = points_spat, buffer = 5000, fun = mean, ID = FALSE)
env_data_scaled <- scale(env_data)

# Environmental PCA - retain first 3 axes
env_pca <- prcomp(env_data_scaled, center = FALSE, scale. = FALSE)
env_pca_scores <- env_pca$x[, 1:3]
variance_explained <- summary(env_pca)$importance[2, 1:3] * 100
cumulative_variance <- sum(variance_explained)

cat("Variance explained by PC1:", round(variance_explained[1], 1), "%\n")
cat("Variance explained by PC2:", round(variance_explained[2], 1), "%\n")
cat("Variance explained by PC3:", round(variance_explained[3], 1), "%\n")
cat("Cumulative variance explained by PC1-3:", round(cumulative_variance, 1), "%\n")

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

# PC loadings and pairwise environmental distances
cat("\nPC loadings (PC1):\n")
print(sort(env_pca$rotation[, 1]))

env_dist <- dist(env_data_scaled)
env_matrix <- as.matrix(env_dist)
rownames(env_matrix) <- colnames(env_matrix) <- sample_codes
cat("\nPairwise environmental distances:\n")
print(round(env_matrix, 2))
