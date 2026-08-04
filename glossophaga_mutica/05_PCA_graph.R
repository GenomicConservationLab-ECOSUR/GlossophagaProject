setwd("Your/directory")
# Run in this in your R installation to produce a PCA figure

library(ggplot2)
library(ggrepel)

# Load PCA eigenvectors
pca <- read.table("mutica_pca_results.eigenvec", header = FALSE) # From 04_plink.sh results.
colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# Map raw BAM-derived sample names to publication codes
code_map <- c(
  "Concha"    = "ConPa",
  "Ixtapa"    = "IxtPa",
  "Poana"     = "PoaGo",
  "Quilamula" = "QuiGo",
  "Hobonil"   = "HobYu",
  "Nizanda"   = "NizOa",
  "Nayarit"   = "IslNa"   # BAM/sample ID; locality is Islas Marias, Nayarit
)

pca$Sample <- sub("\\.bam$", "", pca$IID)
pca$Sample <- code_map[pca$Sample]

# Eigenvalues / variance explained
eigenvals <- scan("mutica_pca_results.eigenval")
total_variance <- sum(eigenvals)
var_explained <- (eigenvals / total_variance) * 100

cat("Variance explained by PCA axes:\n")
for (i in 1:5) {
  cat(sprintf("PC%d: %.2f%%\n", i, var_explained[i]))
}

# PCA plot with codes as labels
ggplot(pca, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 2, color = "#377eb8") +
  geom_text_repel(
    size = 4,
    color = "black",
    segment.color = "grey40",
    segment.size = 0.3,
    box.padding = 0.8
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  )
