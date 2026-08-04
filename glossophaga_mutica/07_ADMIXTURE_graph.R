setwd("Your directory")

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Ancestry bar plot (K = 5) ---
K_value <- 5
Q <- read.table(paste0("mutica_final_pruned_popgen.", K_value, ".Q")) # From the Plink analysis
fam <- read.table("mutica_final_pruned_popgen.fam")

df <- cbind(fam, Q)
colnames(df) <- c("FamilyID", "IndID", "PatID", "MatID", "Sex", "Phenotype",
                  paste0("Cluster", 1:K_value))

# Map raw BAM-derived sample names to publication codes
code_map <- c(
  "Concha"    = "ConPa",
  "Ixtapa"    = "IxtPa",
  "Poana"     = "PoaGo",
  "Quilamula" = "QuiGo",
  "Hobonil"   = "HobYu",
  "Nizanda"   = "NizOa",
  "Nayarit"   = "IslNa"
)

df$IndID <- gsub("\\.bam$", "", df$IndID)
df$IndID <- code_map[df$IndID]

# Geographic order for plotting
geo_order <- c("IslNa", "PoaGo", "QuiGo", "ConPa", "IxtPa", "NizOa", "HobYu")
df$IndID <- factor(df$IndID, levels = geo_order)

df_long <- df %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Cluster",
               values_to = "Ancestry")

admix_plot <- ggplot(df_long, aes(x = IndID, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "italic"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = NULL, x = NULL, y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3")

admix_plot


# Write per-individual ancestry proportions to file for reference in the text
write.csv(df, "admixture_K5_ancestry_proportions.csv", row.names = FALSE)

# --- Cross-validation error plot ---
cv <- read.table("CV.txt", sep = ":", header = FALSE)
cv$K <- as.numeric(gsub("[^0-9]", "", cv$V1))
cv$CVerror <- as.numeric(cv$V2)  # "-nan" (K=6, failed to converge) becomes NA here

cv_valid <- cv[is.finite(cv$CVerror), ]

cv_plot <- ggplot(cv_valid, aes(x = K, y = CVerror)) +
  geom_line() +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "Cross-Validation Error") +
  scale_x_continuous(breaks = cv$K) +
  ggtitle("") +
  geom_text(aes(label = round(CVerror, 3)), vjust = -0.7, size = 4)

cv_plot

cat("Optimal K (lowest CV error):", cv_valid$K[which.min(cv_valid$CVerror)],
    "with CV error =", min(cv_valid$CVerror), "\n")
cat("K values that failed to converge:", cv$K[!is.finite(cv$CVerror)], "\n")
