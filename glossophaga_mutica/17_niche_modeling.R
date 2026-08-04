# AUTHOR: Jesús Antonio Rocamontes Morales
# E-MAIL: jesus.rocamontes@uaz.edu.mx
# INSTITUTION: Universidad Autónoma de Zacatecas
# SCRIPT DESCRIPTION: Species distribution modeling for G. mutica across
# present, LIG, and LGM climate scenarios using kuenm2/maxnet (pure R,
# no MaxEnt GUI/Java required).

setwd("/run/media/rocamontes/EcosurMutica/niche/niche_modeling")

# Packages ----------------------------------------------------------------
packages <- c(
  "rgbif", "dplyr", "CoordinateCleaner", "spThin",
  "sf", "raster", "terra", "usdm", "corrplot", "kuenm2",
  "geodata", "ggplot2", "viridis", "predicts",
  "rnaturalearth", "rnaturalearthdata"
)
lapply(packages, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
})

# Folder setup --------------------------------------------------------------
folders_to_create <- c(
  "layers/raw/LIG_CCSM3",      # PaleoClim LIG (~130 ka), native 2.5 arc-min
  "layers/raw/LGM_CCSM4",      # WorldClim 1.4 LGM (~21 ka), native 2.5 arc-min
  "diagnostics",
  "results"
)
lapply(folders_to_create, dir.create, recursive = TRUE, showWarnings = FALSE)
message("All project folders created.")

# Occurrence data -----------------------------------------------------------
force_occ_rerun  <- FALSE
occ_cache_path   <- "diagnostics/G_mutica_occurrences.csv"

if (!force_occ_rerun && file.exists(occ_cache_path)) {
  occ_cached <- read.csv(occ_cache_path)
  df_final   <- data.frame(
    decimalLongitude = occ_cached$lon,
    decimalLatitude  = occ_cached$lat,
    species          = occ_cached$species
  )
  message(sprintf("Loaded cached occurrences (%d records) from %s",
                  nrow(df_final), occ_cache_path))
} else {
  datosGbif <- occ_search(
    scientificName = "Glossophaga mutica",
    limit          = 13000,
    hasCoordinate  = TRUE
  )
  df <- as.data.frame(datosGbif$data) %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
           basisOfRecord == "PRESERVED_SPECIMEN") %>%
    filter(!duplicated(paste(decimalLatitude, decimalLongitude)))
  df$id <- paste0("id_", seq_len(nrow(df)))
  
  # Coordinate cleaning (flag/exclude problematic coordinates)
  clean_res <- clean_coordinates(
    df, lon = "decimalLongitude", lat = "decimalLatitude", species = "species"
  )
  df_clean <- df[clean_res$.summary, ]
  message(sprintf("Records after coordinate cleaning: %d", nrow(df_clean)))
  
  # Spatial thinning (4 km minimum distance, 10 replicates, keep best)
  set.seed(42)
  thinned <- spThin::thin(
    loc.data = df_clean,
    lat.col  = "decimalLatitude",
    long.col = "decimalLongitude",
    spec.col = "species",
    thin.par = 4, reps = 10,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    verbose = FALSE
  )
  best_rep    <- which.max(sapply(thinned, nrow))
  thin_coords <- thinned[[best_rep]]
  df_final <- df_clean[
    match(
      paste0(round(thin_coords$Latitude, 6), round(thin_coords$Longitude, 6)),
      paste0(round(df_clean$decimalLatitude, 6), round(df_clean$decimalLongitude, 6))
    ), ]
  df_final <- df_final[!is.na(df_final$decimalLatitude), ]
  message(sprintf("Records after thinning: %d", nrow(df_final)))
  
  write.csv(
    data.frame(species = "G_mutica",
               lon = df_final$decimalLongitude,
               lat = df_final$decimalLatitude),
    occ_cache_path,
    row.names = FALSE
  )
  message("Occurrences cached to ", occ_cache_path)
}

# Download environmental layers ---------------------------------------------
ext_m <- ext(-116, -70, 2, 32)

# Present (WorldClim 2.1, native 2.5 arc-minutes)
wc_pres_raw <- geodata::worldclim_global(
  var  = "bio",
  res  = 2.5,
  path = "layers/raw"
)
env_pres <- crop(wc_pres_raw, ext_m)
env_pres <- mask(env_pres, env_pres[[1]] > -Inf)
env_pres <- trim(env_pres)
names(env_pres) <- paste0("bio", 1:19)
rm(wc_pres_raw); gc()

# LIG (~130 ka, CCSM3 via PaleoClim; Brown et al. 2018; www.paleoclim.org)
# Downloaded at native 2.5 arc-minute resolution (matching present-day and
# LGM layers). Files are named bio_1.tif ... bio_19.tif; sort numerically,
# not alphabetically (alphabetic order gives bio_1, bio_10, bio_11...).
lig_files <- list.files("layers/raw/LIG_CCSM3", pattern = "\\.tif$", full.names = TRUE)
if (length(lig_files) >= 19) {
  lig_nums  <- as.integer(gsub(".*bio_(\\d+)\\.tif$", "\\1", lig_files))
  lig_files <- lig_files[order(lig_nums)][1:19]
  env_lig   <- rast(lig_files)
  # PaleoClim LIG layers carry a "naked" WGS84 CRS string without the EPSG:4326
  # code, which causes compareGeom() to fail against present/LGM layers even
  # though the projection is identical. Force EPSG:4326 explicitly.
  crs(env_lig) <- "EPSG:4326"
  env_lig <- crop(env_lig, ext_m)
  # Resample to exactly align grid cells with env_pres. Both are nominally
  # 2.5 arc-minutes, but grid origins can differ slightly between sources;
  # this guarantees cell-for-cell alignment before masking (same rationale
  # as the LGM resample step below).
  env_lig <- resample(env_lig, env_pres, method = "bilinear")
  env_lig <- mask(env_lig, env_pres[[1]])
  env_lig <- trim(env_lig)
  names(env_lig) <- paste0("bio", 1:19)
  # PaleoClim LIG layers store temperature variables (bio1, bio2, bio4-bio11,
  # EXCLUDING bio3) as C x 10, while WorldClim 2.1 (present) stores them
  # directly in C. Divide by 10 so units match the present-day layers.
  # bio3 (Isothermality, 0-100 index) is already on the correct scale.
  temp_vars <- setdiff(paste0("bio", 1:11), "bio3")
  env_lig[[temp_vars]] <- env_lig[[temp_vars]] / 10
  message("LIG layers loaded (native 2.5 arc-min, CCSM3), resampled to align with present-day grid. Temperature variables rescaled (bio3 left unscaled).")
  gc()
} else {
  stop(paste0(
    "LIG layers not found or incomplete in layers/raw/LIG_CCSM3/\n",
    "Download from: http://www.paleoclim.org\n",
    "  -> LIG -> 2.5 arc-minutes -> extract 19 .tif files to layers/raw/LIG_CCSM3/"
  ))
}

# LGM (~21 ka), WorldClim 1.4, CCSM4, native 2.5 arc-minutes
# Download URL: https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_2-5m.zip
# Extract to: layers/raw/LGM_CCSM4/
lgm_files <- list.files("layers/raw/LGM_CCSM4", pattern = "\\.tif$", full.names = TRUE)
if (length(lgm_files) < 19) {
  stop(paste0(
    "LGM CCSM4 layers not found or incomplete in layers/raw/LGM_CCSM4/\n",
    "Download from: https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_2-5m.zip"
  ))
}
# Files are named cclgmbi1.tif ... cclgmbi19.tif; sort numerically, not
# alphabetically (alphabetic order gives bi1, bi10, bi11... which is wrong).
lgm_nums  <- as.integer(gsub(".*?(\\d+)\\.tif$", "\\1", basename(lgm_files)))
lgm_files <- lgm_files[order(lgm_nums)][1:19]
env_lgm   <- rast(lgm_files)
# WorldClim 1.4 LGM layers carry a "naked" WGS84 CRS string without the
# EPSG:4326 code -- force it explicitly to match present and LIG layers.
crs(env_lgm) <- "EPSG:4326"
env_lgm <- crop(env_lgm, ext_m)
# Resample to exactly align grid cells with env_pres. Both are nominally
# 2.5 arc-minutes, but grid origins can differ slightly between sources;
# this guarantees cell-for-cell alignment before masking.
env_lgm <- resample(env_lgm, env_pres, method = "bilinear")
env_lgm <- mask(env_lgm, env_pres[[1]])
env_lgm <- trim(env_lgm)
names(env_lgm) <- paste0("bio", 1:19)
# WorldClim 1.4 stores temperature variables (bio1, bio2, bio4-bio11,
# EXCLUDING bio3) as C x 10. Divide by 10 to match WorldClim 2.1 units.
# bio3 (Isothermality) is already on the correct 0-100 scale.
temp_vars <- setdiff(paste0("bio", 1:11), "bio3")
env_lgm[[temp_vars]] <- env_lgm[[temp_vars]] / 10
message("LGM CCSM4 layers loaded (native 2.5 arc-min), resampled to align with present-day grid. Temperature variables rescaled (bio3 left unscaled).")
gc()

# Variable extraction, outlier diagnostics, and VIF selection ---------------
# Extract present-climate values at occurrence points (all 19 variables)
puntos_sf <- st_as_sf(df_final,
                      coords = c("decimalLongitude", "decimalLatitude"),
                      crs = 4326)
env_vals  <- terra::extract(env_pres, vect(puntos_sf), df = TRUE)
env_vals  <- na.omit(env_vals[, -1])
colnames(env_vals) <- paste0("bio", 1:ncol(env_vals))

# PCA on the full, unfiltered set of 19 bioclimatic variables
pca_init  <- prcomp(env_vals, scale. = TRUE)
pc1_var   <- summary(pca_init)$importance[2, 1]

# Outlier diagnostics: Mahalanobis distance from the first 3 axes of the
# FULL-variable PCA (pca_init), matching the Methods description. This
# must run on the unfiltered variable set, computed BEFORE VIF filtering
# removes any variables, and is independent of which variables VIF later
# retains.
pca_scores_full <- pca_init$x[, 1:3]
m_dist          <- mahalanobis(pca_scores_full, colMeans(pca_scores_full), cov(pca_scores_full))
limit_95        <- qchisq(0.95, df = 3)
message(sprintf("Outlier diagnostic (full 19-variable PCA): %d of %d occurrence points exceed the 95%% Mahalanobis threshold",
                sum(m_dist > limit_95), length(m_dist)))

# VIF threshold set from PC1 variance explained (full-variable PCA)
vif_thresh <- ifelse(pc1_var > 0.7, 5, 10)
message(sprintf("PC1 explains %.1f%% of variance -- VIF threshold set to %d",
                pc1_var * 100, vif_thresh))

# VIF stepwise elimination
vif_res   <- vifstep(env_vals, th = vif_thresh)
vars_keep <- vif_res@results$Variables
message(sprintf("Variables retained after VIF filtering (%d): %s",
                length(vars_keep), paste(vars_keep, collapse = ", ")))

# Apply selection to all time periods
final_pres <- env_pres[[vars_keep]]
final_lig  <- env_lig[[vars_keep]]
final_lgm  <- env_lgm[[vars_keep]]

# Diagnostic PDF --------------------------------------------------------
# PCA on retained (post-VIF) variables, for comparison/visualization only
# -- not used for the outlier diagnostic above.
env_sel_vals <- env_vals[, vars_keep, drop = FALSE]
pca_sel      <- prcomp(env_sel_vals, scale. = TRUE)

pdf("diagnostics/G_mutica_ENM_Diagnostics.pdf", width = 10, height = 8)

# Variable clustering (all 19 vars)
plot(hclust(dist(pca_init$rotation[, 1:3])),
     main = "Variable Clustering (all 19 bioclim)")

# Correlation matrix (retained vars)
corrplot(cor(env_sel_vals), method = "number", type = "upper",
         title = "Correlation Matrix (retained variables)",
         mar = c(0, 0, 2, 0))

# PCA biplots
par(mfrow = c(1, 2))
biplot(pca_init, main = "PCA -- All variables")
biplot(pca_sel,  main = "PCA -- Retained variables")
par(mfrow = c(1, 1))

# Mahalanobis outlier plot (full 19-variable PCA, computed above)
plot(m_dist,
     main = "Mahalanobis Distances (95% CI threshold, full variable set)",
     pch  = 20, ylab = "Distance",
     col  = ifelse(m_dist > limit_95, "red", "black"))
abline(h = limit_95, col = "red", lty = 2)
legend("topright",
       legend = c("Normal", sprintf("Outlier (n=%d)", sum(m_dist > limit_95))),
       col    = c("black", "red"), pch = 20)

dev.off()
message("Diagnostics PDF saved to diagnostics/G_mutica_ENM_Diagnostics.pdf")

# Performs regularization + feature class tuning via maxnet (pure R, no
# Java/MaxEnt GUI required), evaluates candidates via partial ROC, omission
# rate, and AICc, then projects the selected model to past climate scenarios.
# Citation: kuenm2 (Trindade et al. 2026).
#
# CACHING: calibration() takes ~60-90 min. Set force_kuenm2_rerun <- TRUE
# only if occurrence data or variable selection has changed; otherwise the
# cached result is loaded and this step is skipped.
force_kuenm2_rerun <- FALSE
kuenm2_cache_path  <- "diagnostics/kuenm2_calibration.rds"

# Step 1 -- Prepare data: occurrences, background, layers, and tuning grid
occ_df <- data.frame(
  species   = "G_mutica",
  longitude = df_final$decimalLongitude,
  latitude  = df_final$decimalLatitude
)
calib_data <- prepare_data(
  algorithm        = "maxnet",
  occ              = occ_df,
  species          = "species",
  x                = "longitude",
  y                = "latitude",
  raster_variables = final_pres,
  n_background     = 10000,
  features         = c("l", "lq", "lqh", "lqhp", "lqhpt"),
  r_multiplier     = seq(0.5, 4, by = 0.5),
  partition_method = "kfolds",
  n_partitions     = 5,
  seed             = 42
)

# Step 2 -- Calibration: fits all candidate models and evaluates them
# Selection criteria: significant pROC (p < 0.05), OR10, and lowest AICc
if (!force_kuenm2_rerun && file.exists(kuenm2_cache_path)) {
  cal_results <- readRDS(kuenm2_cache_path)
  message("Loaded cached kuenm2 calibration from ", kuenm2_cache_path,
          " -- skipping tuning step.")
} else {
  cal_results <- calibration(
    data             = calib_data,
    error_considered = 10,
    parallel         = FALSE,
    verbose          = TRUE
  )
  saveRDS(cal_results, kuenm2_cache_path)
  message("kuenm2 calibration cached to ", kuenm2_cache_path)
}

summary(cal_results)
write.csv(cal_results$calibration_results,
          "diagnostics/kuenm2_calibration_results.csv",
          row.names = FALSE)

# Performance plot -- OR10 vs AICc for all candidate models
cal_df <- cal_results$calibration_results$All_results
if (!is.null(cal_df) && nrow(cal_df) > 0) {
  cal_summary <- cal_df %>%
    group_by(ID, R_multiplier, Features, AICc) %>%
    summarise(OR10_mean = mean(Omission_rate_at_10, na.rm = TRUE), .groups = "drop")
  p_cal <- ggplot(cal_summary, aes(x = AICc, y = OR10_mean)) +
    geom_point(alpha = 0.3, color = "#4DA6C8", size = 1.5) +
    geom_point(data = cal_results$selected_models,
               aes(x = AICc, y = Omission_rate_at_10.mean),
               color = "#0A4C8A", size = 3, shape = 17) +
    labs(title = "kuenm2 Calibration -- All candidate models",
         subtitle = "Blue triangle = selected model(s)",
         x = "AICc", y = "Mean omission rate at 10%") +
    theme_classic()
  ggsave("diagnostics/kuenm2_calibration_plot.pdf",
         plot = p_cal, width = 10, height = 6)
  message("kuenm2 calibration plot saved to diagnostics/kuenm2_calibration_plot.pdf")
} else {
  message("Calibration plot skipped -- results table not in expected format.")
}

# Step 3 -- Fit the selected best model(s)
fitted_models <- fit_selected(cal_results)
message(sprintf("Best model(s) fitted: %d selected candidate(s)",
                nrow(cal_results$selected_models)))

print(cal_results$selected_models[, c("R_multiplier", "Features",
                                      "Omission_rate_at_10.mean",
                                      "Mean_AUC_ratio_at_10.mean",
                                      "pval_pROC_at_10.mean",
                                      "AICc", "dAIC")])
write.csv(cal_results$selected_models,
          "diagnostics/kuenm2_selected_models.csv",
          row.names = FALSE)

# Step 4 -- Project to present and past climate scenarios
suit_pres <- predict_selected(fitted_models,
                              new_variables = final_pres,
                              consensus     = "mean",
                              type          = "cloglog")$General_consensus
suit_lig  <- predict_selected(fitted_models,
                              new_variables = final_lig,
                              consensus     = "mean",
                              type          = "cloglog")$General_consensus
suit_lgm  <- predict_selected(fitted_models,
                              new_variables = final_lgm,
                              consensus     = "mean",
                              type          = "cloglog")$General_consensus
message("All projections complete (present, LIG, LGM).")
writeRaster(suit_lgm, "results/LGM_CCSM4_suitability.tif", overwrite = TRUE)

# Model evaluation summary ---------------------------------------------------
best_params <- cal_results$selected_models
best_rm     <- best_params$R_multiplier[1]
best_fc     <- best_params$Features[1]
best_proc   <- best_params$Mean_AUC_ratio_at_10.mean[1]
best_or10   <- best_params$Omission_rate_at_10.mean[1]
best_aicc   <- best_params$AICc[1]

model_summary <- data.frame(
  Parameter = c(
    "Occurrences (after thinning)",
    "Number of variables retained",
    "Variables retained",
    "VIF threshold used",
    "Best regularization multiplier (RM)",
    "Best feature classes (FC)",
    "Partial ROC ratio (mean)",
    "OR10 (10th percentile omission rate)",
    "AICc (best model)",
    "LGM GCM",
    "LIG GCM"
  ),
  Value = c(
    nrow(df_final), length(vars_keep),
    paste(vars_keep, collapse = ", "),
    vif_thresh, best_rm, best_fc,
    round(best_proc, 3),
    round(best_or10, 3),
    round(best_aicc, 2),
    "CCSM4 (WorldClim 1.4, 2.5 arc-min)",
    "CCSM3 (PaleoClim, 2.5 arc-min)"
  )
)
write.csv(model_summary, "diagnostics/model_evaluation_summary.csv", row.names = FALSE)
print(model_summary)

# Binary maps, habitat dynamics, and area calculations -----------------------
# Threshold computed directly from model predictions at occurrence points.
suit_at_occ     <- terra::extract(suit_pres, vect(puntos_sf), df = TRUE)[, 2]
suit_at_occ     <- na.omit(suit_at_occ)
threshold_10pct <- quantile(suit_at_occ, probs = 0.10)
message(sprintf("10th percentile Cloglog threshold: %.4f", threshold_10pct))

make_binary <- function(r, thresh) {
  classify(r, matrix(c(-Inf, thresh, 0,
                       thresh, Inf, 1),
                     ncol = 3, byrow = TRUE))
}
bin_pres <- make_binary(suit_pres, threshold_10pct)
bin_lig  <- make_binary(suit_lig,  threshold_10pct)
bin_lgm  <- make_binary(suit_lgm,  threshold_10pct)

# Habitat dynamics: past vs. present
# 0 = unsuitable both; 1 = past only (contraction); 2 = present only (expansion)
# 3 = suitable in both (stable)
calc_dynamics <- function(bin_past, bin_curr) {
  dyn <- classify(
    bin_past * 2 + bin_curr,
    matrix(c(0, 0,  1, 1,  2, 3,  3, 2), ncol = 2, byrow = TRUE)
  )
  levels(dyn) <- data.frame(
    value = 0:3,
    label = c("Unsuitable", "Contraction", "Stable", "Expansion")
  )
  dyn
}
dyn_lgm <- calc_dynamics(bin_lgm, bin_pres)
dyn_lig <- calc_dynamics(bin_lig, bin_pres)

calc_area_km2 <- function(dyn_raster, label) {
  cell_areas <- cellSize(dyn_raster, unit = "km")
  cats       <- c(0, 1, 2, 3)
  cat_labels <- c("Unsuitable", "Contraction", "Stable", "Expansion")
  areas <- sapply(cats, function(cat) {
    mask_r <- ifel(dyn_raster == cat, cell_areas, NA)
    global(mask_r, "sum", na.rm = TRUE)$sum
  })
  data.frame(Period = label, Category = cat_labels, Area_km2 = round(areas, 2))
}
area_lgm <- calc_area_km2(dyn_lgm, "LGM (CCSM4)")
area_lig <- calc_area_km2(dyn_lig, "LIG (CCSM3)")
niche_table <- rbind(area_lgm, area_lig)
write.csv(niche_table, "results/niche_dynamics_area_table.csv", row.names = FALSE)
print(niche_table)

# Four-class categorical suitability rasters
classify_suitability <- function(r, thresh) {
  classify(r,
           matrix(c(
             0,      thresh, 0,
             thresh, 0.5,    1,
             0.5,    0.75,   2,
             0.75,   1.0,    3
           ), ncol = 3, byrow = TRUE))
}
cat_pres <- classify_suitability(suit_pres, threshold_10pct)
cat_lig  <- classify_suitability(suit_lig,  threshold_10pct)
cat_lgm  <- classify_suitability(suit_lgm,  threshold_10pct)
writeRaster(cat_pres, "results/categorical_present.tif",   overwrite = TRUE)
writeRaster(cat_lig,  "results/categorical_LIG_CCSM3.tif", overwrite = TRUE)
writeRaster(cat_lgm,  "results/categorical_LGM_CCSM4.tif", overwrite = TRUE)
writeRaster(dyn_lgm, "results/dynamics_LGM_vs_present.tif", overwrite = TRUE)
writeRaster(dyn_lig, "results/dynamics_LIG_vs_present.tif", overwrite = TRUE)

# MESS -- Multivariate Environmental Similarity Surface ---------------------
# Identifies areas where the COMBINATION of environmental variables in past
# projections is non-analog relative to the training data. Negative values
# indicate non-analog climates where extrapolation is less reliable
# (Elith et al. 2010).
ref_vals <- env_vals[, vars_keep, drop = FALSE]
mess_lig <- predicts::mess(final_lig, ref_vals)
mess_lgm <- predicts::mess(final_lgm, ref_vals)
writeRaster(mess_lig, "results/MESS_LIG_CCSM3.tif", overwrite = TRUE)
writeRaster(mess_lgm, "results/MESS_LGM_CCSM4.tif", overwrite = TRUE)

pdf("diagnostics/MESS_analysis.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(mess_lig, main = "MESS -- LIG (CCSM3)", col = rev(viridis(100)))
plot(mess_lgm, main = "MESS -- LGM (CCSM4)", col = rev(viridis(100)))
par(mfrow = c(1, 1))
dev.off()

nonanalog_lig <- mess_lig < 0
nonanalog_lgm <- mess_lgm < 0
nonanalog_area_lig <- global(ifel(bin_lig == 1 & nonanalog_lig, cellSize(bin_lig, unit = "km"), NA), "sum", na.rm = TRUE)$sum
nonanalog_area_lgm <- global(ifel(bin_lgm == 1 & nonanalog_lgm, cellSize(bin_lgm, unit = "km"), NA), "sum", na.rm = TRUE)$sum
total_suitable_lgm <- global(ifel(bin_lgm == 1, cellSize(bin_lgm, unit = "km"), NA), "sum", na.rm = TRUE)$sum
total_suitable_lig <- global(ifel(bin_lig == 1, cellSize(bin_lig, unit = "km"), NA), "sum", na.rm = TRUE)$sum

mess_summary <- data.frame(
  Parameter = c(
    "Total suitable area LGM (km2)",
    "Total suitable area LIG (km2)",
    "Non-analog area within LGM suitable range (km2)",
    "Non-analog area within LIG suitable range (km2)",
    "Proportion non-analog (LGM, %)",
    "Proportion non-analog (LIG, %)"
  ),
  Value = c(
    round(total_suitable_lgm, 1),
    round(total_suitable_lig, 1),
    round(nonanalog_area_lgm, 1),
    round(nonanalog_area_lig, 1),
    round(100 * nonanalog_area_lgm / total_suitable_lgm, 1),
    round(100 * nonanalog_area_lig / total_suitable_lig, 1)
  )
)
write.csv(mess_summary, "diagnostics/MESS_summary.csv", row.names = FALSE)
print(mess_summary)

# Categorical suitability maps: clean, individual panels (present, LIG, LGM)
# No titles or legends -- intended as publication-ready base maps.

# White = unsuitable, then 3 shades of blue
suit_colors <- c(
  "#F0F0F0",   # 0 -- unsuitable
  "#A8D5E2",   # 1 -- low suitability
  "#4DA6C8",   # 2 -- high suitability
  "#0A4C8A"    # 3 -- very high suitability
)

ext_plot       <- ext(-116, -70, 2, 32)
coast_crop     <- crop(ne_coastline(scale = "medium", returnclass = "sv"), ext_plot)
countries_crop <- crop(ne_countries(scale = "medium", returnclass = "sv"), ext_plot)

plot_clean <- function(cat_raster, filename) {
  cat_raster <- crop(cat_raster, ext_plot)
  # When a raster contains only ONE of the four categories (e.g. all cells
  # unsuitable), terra's plot() with explicit breaks can map that single
  # value to the wrong color. Force the full 0-3 range to always be present
  # by setting four corner pixels (outside the visible data extent) to
  # values 0-3, guaranteeing consistent color mapping regardless of which
  # categories are actually present.
  corner_cells <- cellFromRowCol(cat_raster, row = 1, col = 1:4)
  cat_raster[corner_cells] <- 0:3
  
  pdf(file.path("results", filename), width = 7, height = 5.5, useDingbats = FALSE)
  par(mar = c(3, 3, 0.5, 0.5), oma = c(0, 0, 0, 0))
  plot(cat_raster,
       col      = suit_colors,
       breaks   = c(-0.5, 0.5, 1.5, 2.5, 3.5),
       legend   = FALSE,
       axes     = TRUE,
       box      = TRUE,
       xlim     = c(-116, -70),
       ylim     = c(2, 32),
       cex.axis = 0.7,
       col.axis = "#444444")
  lines(coast_crop,     col = "#555555", lwd = 0.7)
  lines(countries_crop, col = "#888888", lwd = 0.3, lty = 2)
  dev.off()
}

plot_clean(cat_pres, "clean_present.pdf")
plot_clean(cat_lig,  "clean_LIG_CCSM3.pdf")
plot_clean(cat_lgm,  "clean_LGM_CCSM4.pdf")
message("Clean suitability maps saved to results/ (present, LIG, LGM).")
message("=== All done. Outputs in results/ and diagnostics/ ===")
