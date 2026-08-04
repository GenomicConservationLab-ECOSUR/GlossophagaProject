setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026")
library(terra)
library(factoextra)

localities <- read.table("coords.txt", header = TRUE, sep = "\t")

# The coordinates should be these:

# Name	Longitude	Latitude
# concha	-92.159444	15.004194
# nizanda	-94.995303	16.654892
# quilamula	-99.019915	18.510747
# hobonil	-89.02743	20.00728
# poana	-92.757623	17.535882
# ixtapa	-101.562611	17.657333

# Exclude Islas Marias/Nayarit, this locality is analyzed separately (Island vs. Mainland comparison)
# and should not be included in the environmental-cluster PCA/UPGMA at all,
# so that the resulting mainland clustering (including Nizanda's status)
# is a genuine result of clustering the 6 mainland samples alone, not a
# leftover grouping from a 7-sample run with one point removed after the fact.

# Load environmental layers
bio_path <- "/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/env_layers"
bio_files <- list.files(bio_path, pattern = "wc2.1_2.5m_bio_.*\\.tif$", full.names = TRUE)
env_stack <- rast(bio_files)

names(env_stack) <- gsub("wc2.1_2.5m_bio_", "bio", names(env_stack))
target_vars <- c("bio2", "bio3", "bio8", "bio9", "bio13", "bio14", "bio15", "bio18", "bio19")
env_subset <- env_stack[[target_vars]]

ext_values <- extract(env_subset, localities[, c("Longitude", "Latitude")])
pca_ready <- na.omit(cbind(localities, ext_values))
pca_matrix <- pca_ready[, target_vars]
rownames(pca_matrix) <- pca_ready$Name

res.pca <- prcomp(pca_matrix, scale. = TRUE)

fviz_pca_var(res.pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +
  labs(title = "Environmental Variable Contributions")

fviz_pca_biplot(res.pca,
                geom.ind = c("point", "text"),
                label = "ind",
                col.ind = "black",
                pointshape = 19,
                pointsize = 2,
                geom.var = c("arrow", "text"),
                col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,
                alpha.var = 0.5) +
  theme_minimal() +
  labs(title = "G. mutica Environmental Biplot")

# Mahalanobis distance on first 2 PCs
pc_scores <- res.pca$x[, 1:2]
cov_pc <- cov(pc_scores)
mahal_dist <- as.dist(apply(pc_scores, 1, function(i) {
  apply(pc_scores, 1, function(j) sqrt(mahalanobis(i, j, cov_pc)))
}))

# Hierarchical clustering (UPGMA) -- now on 6 mainland samples only
hc_env <- hclust(mahal_dist, method = "average")

plot(hc_env, main = "",
     xlab = "G. mutica", sub = "", ylab = "Distance")

# k=5 clusters among 6 samples -- one pair merges, the rest stand alone
rect.hclust(hc_env, k = 5, border = "red")

fviz_pca_ind(res.pca, geom.ind = c("point", "text"), label = "ind",
             repel = TRUE, habillage = as.factor(pca_ready$Name),
             pointshape = 19, pointsize = 2) +
  theme_minimal() + theme(legend.position = "none")

# pixy_map generation stays the same -- no need to change, since it
# already reflects the intended 5-cluster mainland grouping
pixy_map <- data.frame(
  sample_id = c("Concha.bam", "Poana.bam", "Quilamula.bam",
                "Hobonil.bam", "Ixtapa.bam", "Nizanda.bam"),
  population_id = c("Cluster_Concha", "Cluster_Poana", "Cluster_Quilamula",
                    "Cluster_Hobo_Ixta", "Cluster_Hobo_Ixta", "Cluster_Nizanda")
)
write.table(pixy_map,
            file = "master_map_env.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
