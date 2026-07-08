

####--------- 1. Environment Setup & Package Loading ---------####

# Check current R session information (version, platform, loaded packages)
sessionInfo()

# Clear the workspace to ensure a clean analysis environment
rm(list = ls())

# Load required packages
library(ENMeval)    # MaxEnt hyperparameter tuning and model evaluation
library(sf)         # Spatial data handling (shapefiles, geometry operations)
library(data.table) # Efficient data manipulation for large datasets
library(tidyverse)  # Data wrangling and visualization (dplyr, ggplot2, tidyr included)
library(usdm)       # Multicollinearity evaluation (VIF calculation for variable selection)
library(ecospat)    # ENM/Niche analysis tools 

# Fix random seed to ensure reproducibility across the entire workflow
# (e.g., background sampling, k-fold partitioning)
set.seed(777)


####--------- 2. Data Loading & Preprocessing ---------####


# Some stream sections near the border with North Korea lack MK-PRISM climate data,
# resulting in NA values for climatic variables.
# In addition, NA or Inf values may occasionally occur during the calculation of physical environmental variables.
# To remove these missing or abnormal values in a consistent manner,
# filtering was performed using p_bio01, a representative climatic variable, as the reference.
# Most NA or Inf values in physical environmental variables occur during the calculation of coefficients of variation (CV),
# rather than during the calculation of means or medians.
# These values will be removed or corrected separately in later analytical steps
# if additional effects are detected.


# Read Data 
# Load fish occurrence data
presencedata <- st_read("./data_raw/R.steindachneri_occ.gpkg", quiet = FALSE) %>%
  filter(BBSNCD %in% c("10", "13")) %>%
  st_transform(4326) %>%
  dplyr::select(-any_of(c("cv_stst", "cv_ct_slp", "cv_ctpo", "cv_ctst", "cv_stslp"))) %>%
  mutate(str_ID = as.character(str_ID), slp_CV = std_stslp / mean_stslp) %>%
  rename(ws_area = wsarea, elev_mean = meanstdem, slp_mean = mean_stslp, st_power_mean = mean_stpo, urban = per_ct_urb, agri = per_ct_agr, forest = per_ct_for) %>%
  distinct(str_ID, .keep_all = TRUE) %>%
  filter(!is.na(p_bio01))


# Load stream network data
streamnetwork <- st_read("./data_raw/stream_network_Rhynchocypris.gpkg", quiet = FALSE) %>%
  st_transform(4326) %>% # convert coordinate system into wgs1984 (epsg4326) 
  dplyr::select(-c(cv_stst, cv_ct_slp, cv_ctpo, cv_ctst, cv_stslp)) %>%
  mutate(slp_CV = std_stslp / mean_stslp) %>%
  rename(ws_area = wsarea, elev_mean= meanstdem, slp_mean= mean_stslp, st_power_mean= mean_stpo,
         urban= per_ct_urb, agri= per_ct_agr, forest= per_ct_for) %>% 
  distinct(str_ID, .keep_all = TRUE) %>%
  filter(!is.na(p_bio01))  


# Select stream segments for prediction and visualization
# (retain only small-basins (SBSNCD) where occurrences are recorded)
streamnetwork2 <- streamnetwork %>%
  filter(SBSNCD %in% presencedata$SBSNCD)

streamnetwork3 <- streamnetwork %>%
  filter(SBSNCD %in% presencedata$SBSNCD) %>%
  st_drop_geometry()


# Background dataset for distribution prediction (MaxEnt model projection)
bgdata <- streamnetwork %>%
  filter(!str_ID %in% presencedata$str_ID) %>%   # Remove presence segments
  filter(SBSNCD %in% presencedata$SBSNCD)        # Keep only small-basins within study area


# Load nationwide watershed boundary shapefile
# Used as a base map layer for visualization
korea_watershed <- st_read("./data_raw/WKMBBSN.gpkg") %>%
  st_transform(4326)   # Convert CRS to WGS84 (EPSG:4326)


# Loadboundary for the Han River watershed
# Extract only target basins (BBSNCD: 10, 13) for model visualization and analysis
han_watershed <- korea_watershed %>%
  st_transform(4326) %>%         # CRS transformation
  filter(BBSNCD %in% c("10","13"))  # Select study basins (Han River system)


####--------- 3. BioClim Variable Description ---------####

# Bio Data Explanation #

#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (×100)
#BIO4 = Temperature Seasonality (standard deviation ×100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#bio17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#bio19 = Precipitation of Coldest Quarter


####--------- 4-1 Variable Selection & Correlation Analysis ---------####

# Check multicollinearity among predictor variables using VIF
# High VIF indicates multicollinearity.
# Using highly correlated predictors may cause unstable parameter estimates
# (VIF > 10 generally indicates strong multicollinearity)
# We adopted a strict threshold (VIF < 4) when selecting final variables to avoid multicollinearity issues.


# Refined VIF check after variable reduction
# (Only variables satisfying VIF < 4 were retained for modeling)
usdm::vif(dplyr::select(streamnetwork3, slp_mean, slp_CV, ws_area,
                        forest, weir_dens,
                        p_bio05, p_bio12, p_bio13))


####--------- 4-2 Pearson Correlation & Significance Test ---------####

# install.packages("GGally") # Install the 'GGally' package
library(GGally)  # used for correlation matrix visualization

# Select variables for correlation analysis
cor1 <- streamnetwork3 %>%
  dplyr::select(elev_mean, slp_mean, slp_CV, ws_area,
                forest, agri, weir_dens,
                p_bio01, p_bio02, p_bio03, p_bio04, p_bio05,
                p_bio06, p_bio07, p_bio08, p_bio09, p_bio10,
                p_bio11, p_bio12, p_bio13, p_bio14, p_bio15,
                p_bio16, p_bio17, p_bio18, p_bio19) %>%
  rename_with(
    ~ gsub("^p_bio", "Bio", ., ignore.case = TRUE),
    starts_with("p_bio")
  )


# Variable order
vars <- c("elev_mean", "ws_area", "slp_mean", "slp_CV", 
          "forest", "agri", "weir_dens",
          "Bio01", "Bio02", "Bio03", "Bio04", "Bio05",
          "Bio06", "Bio07", "Bio08", "Bio09", "Bio10",
          "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
          "Bio16", "Bio17", "Bio18", "Bio19")

x <- cor1 %>%
  dplyr::select(all_of(vars)) 

# Calculate correlation coefficients and p-values
rmat <- cor(x, use = "complete.obs")
pmat <- sapply(seq_along(vars), function(i) {
  sapply(seq_along(vars), function(j) {
    cor.test(x[[i]], x[[j]])$p.value
  })
})

df <- expand.grid(Var1 = vars, Var2 = vars) %>%
  mutate(
    i = match(Var1, vars),
    j = match(Var2, vars),
    r = rmat[cbind(i, j)],
    p = pmat[cbind(i, j)],
    sig = case_when(
      p < .001 ~ "***",
      p < .01  ~ "**",
      p < .05  ~ "*",
      TRUE     ~ ""
    ),
    Var1 = factor(Var1, levels = vars),
    Var2 = factor(Var2, levels = vars)
  ) %>%
  filter(i <= j)   

ggplot(df, aes(Var1, Var2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4, vjust = 0.6) +
  geom_text(aes(label = sig), size = 4, vjust = -0.6) +
  
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  
  coord_fixed() +
  
  scale_x_discrete(
    position = "bottom",
    expand = c(0, 0),
    drop = FALSE
  ) +
  
  scale_y_discrete(
    limits = rev(vars),
    expand = c(0, 0),
    drop = FALSE
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    panel.grid = element_blank(),
    plot.caption = element_text(
      size = 10,
      hjust = 0,
      margin = margin(t = 8),
      colour = "black"
    )
  ) +
  
  labs(
    caption = "Significance levels: * p < 0.05, ** p < 0.01, *** p < 0.001"
  )




####--------- 5. Occurrence & Background Data Preparation ---------####

# Prepare occurrence dataset for MaxEnt modeling 
occ <- data.frame(dplyr::select(presencedata, 
                                longitude= mean_lon, latitude= mean_lat, 
                                slp_mean, slp_CV, ws_area, 
                                forest, weir_dens,
                                p_bio05, p_bio12, p_bio13)) %>% # Select only the variables needed for modeling
  dplyr::select(-geom) %>% # remove geometry
  as_tibble() %>% # make as a tibble format
  na.omit() # remove missing data


# check occurrence point location
occ_data_fig <- st_as_sf(occ, coords= c("longitude", "latitude"), crs= 4326) # Convert occurrence data to spatial points 


# Species Occurrence Visualization 
occurrence_plot <- ggplot() +
  
  ## (1) Basin boundary 
  geom_sf(
    data = han_watershed,
    colour = "gray30",
    fill   = NA,
    linewidth = 0.8
  ) +
  
  ## (2) Stream network 
  geom_sf(
    data = streamnetwork,
    aes(linewidth = st_order),
    color = "gray75",
    alpha = 0.7
  ) +
  scale_linewidth_continuous(
    range = c(0.3, 1.2),
    guide = "none"
  ) +
  
  ## (3) Occurrence points
  geom_sf(
    data = occ_data_fig,
    shape = 21,
    size  = 1.5,
    fill  = "red",   
    color = "black",
    stroke = 0.4
  ) +
  
  ## (4) Coordinate system 
  coord_sf(expand = FALSE) +
  
  ## (5) Theme 
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.title      = element_blank(),
    plot.caption    = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  )


occurrence_plot



# Background concept #
# In MaxEnt (presence-only SDM), 'background' points represent the available environmental space
# within the species' accessible area. They do not indicate absence, but serve as environmental
# contrast against presence points to estimate the species' niche.
# Proper background selection is crucial because biased or overly broad background areas may lead
# to overprediction or inflated suitability. Hence, background points were sampled only within the
# basin where occurrence records exist.


# Background size in MaxEnt is commonly set to ~10,000 (default) or 10× presence records,
# especially in raster-based studies where environmental space is continuous.
# However, our study uses stream-segment units (limited total segments = 8,443 in Han River system),
# and modeling was restricted to sub-basins accessible to the species.
# Therefore, background = available segments within presence basins instead of a fixed large number.

# Prepare background dataset for MaxEnt modeling
bgdata2 <- data.frame(dplyr::select(bgdata, 
                                     longitude= mean_lon, latitude= mean_lat, 
                                     slp_mean, slp_CV, ws_area, 
                                     forest, weir_dens,
                                     p_bio05, p_bio12, p_bio13)) %>% # Select relevant variables for background
  as.data.frame() %>%
  dplyr::select(-geom) %>% # remove geometry
  as_tibble() %>% # make as a tibble format
  na.omit() # remove missing data


# Background Sampling for MaxEnt #
# Randomly sample background points used as environmental availability
if(nrow(bgdata2)<973) # Use all segments if fewer than target sample size
{
  samplesize = nrow(bgdata2)
} else {
  samplesize = 973 # Adjustable sampling size depending on study design
}
bg = bgdata2[sample(1:nrow(bgdata2), samplesize),] # Random sampling of background segments



# check occurrence and background point location
bg_data_fig <- st_as_sf(bg, coords= c("longitude", "latitude"), crs= 4326) # Convert background points to spatial points

# Species Occurrence vs Background Visualization
occ_bg_plot <- ggplot() +
  # Basin boundaries
  geom_sf(data = han_watershed, colour = "gray30", fill = NA, linewidth = 0.5) +
  
  # Stream network with thickness proportional to stream order
  geom_sf(data = streamnetwork, aes(linewidth = st_order), color = "steelblue", alpha = 0.6) +
  scale_linewidth_continuous(range = c(0.3, 1.5), guide = "none") +
  
  # Background points
  geom_sf(data = bg_data_fig, shape = 21, size = 1.5, fill = "blue", color = "black", stroke = 0.3, alpha = 0.9) +
  
  # Occurrence points
  geom_sf(data = occ_data_fig, shape = 21, size = 1.5, fill = "red", color = "black", stroke = 0.3, alpha = 0.9) +
  
  # Coordinate axes
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(min(st_coordinates(korea_watershed)[,1])),
                                  ceiling(max(st_coordinates(korea_watershed)[,1])), by = 1)) +
  scale_y_continuous(breaks = seq(floor(min(st_coordinates(korea_watershed)[,2])),
                                  ceiling(max(st_coordinates(korea_watershed)[,2])), by = 1)) +
  
  # Labels and theme
  labs(title = "Species Occurrence vs Background",
       caption = "Red points: Occurrence\nBlue points: Background\nStream thickness represents stream order") +
  
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, face = "italic"),
        plot.caption = element_text(size = 9, color = "gray30"),
        axis.title = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none")

occ_bg_plot


####--------- 6. MaxEnt Model Tuning & Evaluation ---------####

# Feature Class (FC) = Set of feature transformations used to shape response curves
#  - L: Linear, Q: Quadratic, H: Hinge, P: Product, T: Threshold
#  - More complex FC combinations increase model flexibility but may lead to overfitting
#
# Regularization Multiplier (RM) = Penalty controlling model complexity
#  - Low RM → detailed model, risk of overfitting
#  - High RM → smoother response, better generalization

# MaxEnt modeling using ENMevaluate
# - Tune hyperparameters: Feature classes (FC) & Regularization multiplier (RM)
# - FC tested: L, LQ, H, LQH, LQHP, LQHPT
# - RM tested: 1 to 5

# Set tuning parameters for the model (Hyperparameter setting)
tune.args <- list(fc= c("L", "LQ","H", "LQH","LQHP","LQHPT"),   rm = 1:5) # model tuning parameter range


# Run MaxEnt with spatial "block" cross-validation
# Random k-fold may inflate model performance due to spatial proximity,
# whereas block partitioning ensures spatial independence, providing a more realistic evaluation.
mod_res <- ENMevaluate(
  occ,
  bg = bg,
  algorithm = "maxent.jar",
  tune.args = tune.args,
  partitions = "block",       
  parallel = FALSE
)


# Check performance summary for all FC × RM model candidates
results_all <- eval.results(mod_res)
head(results_all)

results_all <- results_all %>%
  mutate(
    delta_auc_val = max(auc.val.avg) - auc.val.avg,
    delta_cbi_val = max(cbi.val.avg) - cbi.val.avg
  )


# Visualize key evaluation metrics
# - or.mtp : omission rate at minimum training presence
# - auc.val: validation AUC (discrimination ability)
# - cbi.val: Continuous Boyce Index (prediction reliability)
evalplot.stats(e = mod_res,
               stats = c("or.mtp","auc.val","cbi.val"),
               color = "fc", x.var = "rm", error.bars = FALSE)

# Model selection and evaluation 

# Extract evaluation results for all candidate FC × RM models
sel_mod <- eval.results(mod_res)

# Select the optimal model with the lowest AICc (delta.AICc = 0)
#  - AICc balances goodness-of-fit and model complexity
#  - delta.AICc = 0 indicates the best-supported model among candidates
sel_opt.seq <- sel_mod %>%
  filter(delta.AICc == 0)

# Display the selected optimal model
sel_opt.seq  


# R package for running MaxEnt 
library(dismo)


# Define environmental variable used in MaxEnt
env_vars <- c("slp_mean", "slp_CV", "ws_area", 
              "forest", "weir_dens", "p_bio05", "p_bio12", "p_bio13")

# Combine environmental values of presence (occ) and background (bg)
# MaxEnt compares environmental conditions between presence vs background
env_layers <- rbind(occ[, env_vars], bg[, env_vars])

# Create presence/background indicator vector
# presence = 1, background = 0
# number of rows must match env_layers (critical)
p_indicator <- c(rep(1, nrow(occ)), rep(0, nrow(bg)))


# Run MaxEnt model using dismo::maxent
mod_maxent <- maxent(
  x = env_layers,        # Environmental variable matrix (presence + background)
  p = p_indicator,       # Presence (1) / Background (0) indicator
  args = c(
    "replicates=10",               # 10 replicates = 10-fold cross-validation
    "replicatetype=crossvalidate", # Apply cross-validation
    "randomseed=false",            # If TRUE, random seed changes each run (FALSE keeps consistent split)  
    # --- Regularization & Feature Class configuration ---
    "betamultiplier=3",            # Regularization Multiplier = 5 (controls model complexity)
    "autofeature=false",
    "linear=true",
    "quadratic=true",
    "hinge=true",
    "product=false",
    "threshold=false",
    # --- Output settings ---
    "jackknife=true",              # Produce jackknife variable importance test
    "responsecurves=true"         # Generate response curves
  ),
  path = "./data_outcome/R.steindachneri_maxent_result"
)

# Save trained MaxEnt model as .rds
saveRDS(mod_maxent, file = "./data_outcome/R.steindachneri_maxent_result.rds")

# Load saved model
mod_maxent <- readRDS("./data_outcome/R.steindachneri_maxent_result.rds")


mod_maxent


####--------- 7. Current Distribution Prediction (dismo::maxent replicates) ---------####

# (1) Build prediction input table (one row = one stream segment)
proj_with_id <- streamnetwork2 %>%
  st_drop_geometry() %>%
  as.data.frame() %>% 
  dplyr::select(
    str_ID,
    longitude = mean_lon,
    latitude  = mean_lat,
    slp_mean, slp_CV, ws_area,
    forest, weir_dens,
    p_bio05, p_bio12, p_bio13
  ) %>%
  na.omit() %>%
  distinct(str_ID, .keep_all = TRUE)

# (2) Prepare input data for prediction
pred_input <- proj_with_id[, c(
  "slp_mean", "slp_CV", "ws_area",
  "forest", "weir_dens",
  "p_bio05", "p_bio12", "p_bio13"
)] %>% as.data.frame()

# (3) Predict suitability for each MaxEnt replicate (n = 10)
pred_list <- lapply(mod_maxent@models, function(m) {
  predict(
    object = m,
    x = pred_input,
    args = c("outputformat=cloglog")  # cloglog scale
  )
})

# (4) Compute mean suitability across replicates
pred_mean <- rowMeans(do.call(cbind, pred_list), na.rm = TRUE)

# (5) Combine predictions with segment IDs
pred_df <- proj_with_id %>%
  as.data.frame() %>%                      
  dplyr::mutate(probability = as.numeric(pred_mean)) %>% 
  dplyr::select(str_ID, probability)       

# (6) Attach suitability value to segment ID
streamnetwork_pred_present <- streamnetwork %>%
  left_join(pred_df, by = "str_ID")

# (7) Plot predicted present habitat suitability
present_pred_fig <- ggplot() +
  geom_sf(data = han_watershed, colour = "gray30", fill = NA, linewidth = 0.8) +
  geom_sf(
    data = streamnetwork_pred_present,
    aes(color = probability, linewidth = st_order),
    show.legend = c(color = TRUE, linewidth = FALSE)
  ) +
  scale_color_gradientn(
    colors = c("#a6cee3", "#fdae61", "#d7191c"),   # light blue → orange → red
    name = "Habitat suitability",
    limits = c(0, 1),
    na.value = "gray90"
  ) +
  scale_linewidth_continuous(range = c(0.3, 1.2), guide = "none") +
  coord_sf(expand = FALSE) +
  labs(title = "Predicted Habitat Suitability (Historical)") +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.3, "cm"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# (8) Render the map
present_pred_fig


####--------- MaxEnt prediction for all future SSP scenarios ---------####


# Multi-replicate mean prediction for all SSP scenarios 
library(dplyr)

# Shared non-climatic environmental variables
base_vars <- c("slp_mean", "slp_CV", "ws_area", "forest", "weir_dens")

# Scenario-specific climate variables (BIO05, BIO12, BIO13)
# s1=SSP1-2.6, s2=SSP2-4.5, s5=SSP5-8.5
# 21_40=2021–2040 ... 81_100=2081–2100
scenarios <- list(
  s1_21_40 = c("s1_24_b05", "s1_24_b12", "s1_24_b13"),
  s1_41_60 = c("s1_46_b05", "s1_46_b12", "s1_46_b13"),
  s1_61_80 = c("s1_68_b05", "s1_68_b12", "s1_68_b13"),
  s1_81_100 = c("s1_81_b05", "s1_81_b12", "s1_81_b13"),
  
  s2_21_40 = c("s2_24_b05", "s2_24_b12", "s2_24_b13"),
  s2_41_60 = c("s2_46_b05", "s2_46_b12", "s2_46_b13"),
  s2_61_80 = c("s2_68_b05", "s2_68_b12", "s2_68_b13"),
  s2_81_100 = c("s2_81_b05", "s2_81_b12", "s2_81_b13"),
  
  s5_21_40 = c("s5_24_b05", "s5_24_b12", "s5_24_b13"),
  s5_41_60 = c("s5_46_b05", "s5_46_b12", "s5_46_b13"),
  s5_61_80 = c("s5_68_b05", "s5_68_b12", "s5_68_b13"),
  s5_81_100 = c("s5_81_b05", "s5_81_b12", "s5_81_b13")
)

# container for output objects
predict_results <- list()

# Run prediction for each SSP × period combination
for (name in names(scenarios)) {
  vars <- scenarios[[name]]
  message("Running prediction for: ", name)
  
  # (1) Construct segment-wise future environment table
  proj_with_id <- streamnetwork2 %>%
    st_drop_geometry() %>%
    as.data.frame() %>%
    dplyr::select(
      str_ID,
      longitude = mean_lon,
      latitude  = mean_lat,
      all_of(base_vars),
      p_bio05 = vars[1],
      p_bio12 = vars[2],
      p_bio13 = vars[3]
    ) %>%
    stats::na.omit() %>%
    dplyr::distinct(str_ID, .keep_all = TRUE)
  
  # (2) Matrix for prediction input
  pred_input <- proj_with_id[, c(
    "slp_mean", "slp_CV", "ws_area",
    "forest", "weir_dens",
    "p_bio05", "p_bio12", "p_bio13"
  )] %>% as.data.frame()
  
  # (3) Predict from each replicate model (10 models)
  pred_list <- lapply(mod_maxent@models, function(m) {
    predict(m, x = pred_input, args = c("outputformat=cloglog"))
  })
  
  # (4) Mean suitability across replicates
  pred_mean <- rowMeans(do.call(cbind, pred_list), na.rm = TRUE)
  
  # (5) Combine segment ID + predicted probability
  pred_df <- proj_with_id %>%
    dplyr::mutate(probability = as.numeric(pred_mean)) %>%
    dplyr::select(str_ID, probability)
  
  # (6) Attach prediction to spatial stream network
  streamnetwork_pred <- streamnetwork %>%
    dplyr::left_join(pred_df, by = "str_ID")
  
  # (7) Save results
  predict_results[[name]] <- streamnetwork_pred
}

# Export each scenario object for later visualization
streamnetwork_pred_s1_21_40 <- predict_results[["s1_21_40"]]
streamnetwork_pred_s1_41_60 <- predict_results[["s1_41_60"]]
streamnetwork_pred_s1_61_80 <- predict_results[["s1_61_80"]]
streamnetwork_pred_s1_81_100 <- predict_results[["s1_81_100"]]

streamnetwork_pred_s2_21_40 <- predict_results[["s2_21_40"]]
streamnetwork_pred_s2_41_60 <- predict_results[["s2_41_60"]]
streamnetwork_pred_s2_61_80 <- predict_results[["s2_61_80"]]
streamnetwork_pred_s2_81_100 <- predict_results[["s2_81_100"]]

streamnetwork_pred_s5_21_40 <- predict_results[["s5_21_40"]]
streamnetwork_pred_s5_41_60 <- predict_results[["s5_41_60"]]
streamnetwork_pred_s5_61_80 <- predict_results[["s5_61_80"]]
streamnetwork_pred_s5_81_100 <- predict_results[["s5_81_100"]]


# Prediction Map Function 
plot_pred_map <- function(stream_data, prob_col, title_text) {
  ggplot() +
    # (1) Basin boundary for spatial context
    geom_sf(data = han_watershed, colour = "gray30", fill = NA, linewidth = 0.8) +
    
    # (2) Predicted suitability (segment-based)
    geom_sf(data = stream_data, aes(color = !!sym(prob_col), linewidth = st_order), show.legend = FALSE) +
    
    # (3) Color gradient: low (cool) → high (warm)
    scale_color_gradientn(
      colors = c("#a6cee3", "#fdae61", "#d7191c"),
      name = "Habitat suitability",
      limits = c(0, 1),
      na.value = "gray90"
    ) +
    
    # (4) Stream order-based linewidth
    scale_linewidth_continuous(range = c(0.3, 1.2)) +
    
    # (5) Coordinate system and extent
    coord_sf(expand = FALSE) +
    scale_x_continuous(
      breaks = seq(floor(min(st_coordinates(korea_watershed)[,1])),
                   ceiling(max(st_coordinates(korea_watershed)[,1])), by = 1)
    ) +
    scale_y_continuous(
      breaks = seq(floor(min(st_coordinates(korea_watershed)[,2])),
                   ceiling(max(st_coordinates(korea_watershed)[,2])), by = 1)
    ) +
    
    # (6) Titles and labels
    labs(title = title_text) +
    
    # (7) Theme identical to hist_pred_fig
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# SSP126 
ssp1_21_40_pred_fig <- plot_pred_map(streamnetwork_pred_s1_21_40, "probability", "SSP126 2021–2040 Prediction")
ssp1_41_60_pred_fig <- plot_pred_map(streamnetwork_pred_s1_41_60, "probability", "SSP126 2041–2060 Prediction")
ssp1_61_80_pred_fig <- plot_pred_map(streamnetwork_pred_s1_61_80, "probability", "SSP126 2061–2080 Prediction")
ssp1_81_100_pred_fig <- plot_pred_map(streamnetwork_pred_s1_81_100, "probability", "SSP126 2081–2100 Prediction")

# SSP245
ssp2_21_40_pred_fig <- plot_pred_map(streamnetwork_pred_s2_21_40, "probability", "SSP245 2021–2040 Prediction")
ssp2_41_60_pred_fig <- plot_pred_map(streamnetwork_pred_s2_41_60, "probability", "SSP245 2041–2060 Prediction")
ssp2_61_80_pred_fig <- plot_pred_map(streamnetwork_pred_s2_61_80, "probability", "SSP245 2061–2080 Prediction")
ssp2_81_100_pred_fig <- plot_pred_map(streamnetwork_pred_s2_81_100, "probability", "SSP245 2081–2100 Prediction")

# SSP585
ssp5_21_40_pred_fig <- plot_pred_map(streamnetwork_pred_s5_21_40, "probability", "SSP585 2021–2040 Prediction")
ssp5_41_60_pred_fig <- plot_pred_map(streamnetwork_pred_s5_41_60, "probability", "SSP585 2041–2060 Prediction")
ssp5_61_80_pred_fig <- plot_pred_map(streamnetwork_pred_s5_61_80, "probability", "SSP585 2061–2080 Prediction")
ssp5_81_100_pred_fig <- plot_pred_map(streamnetwork_pred_s5_81_100, "probability", "SSP585 2081–2100 Prediction")


# Multiple Figure 
library(ggpubr)


# Arrange all figures (SSP scenarios) in a 4x4 grid
ggarrange(ssp1_21_40_pred_fig, ssp1_41_60_pred_fig, ssp1_61_80_pred_fig, ssp1_81_100_pred_fig,
          ssp2_21_40_pred_fig, ssp2_41_60_pred_fig, ssp2_61_80_pred_fig, ssp2_81_100_pred_fig,
          ssp5_21_40_pred_fig, ssp5_41_60_pred_fig, ssp5_61_80_pred_fig, ssp5_81_100_pred_fig,
          ncol = 4, nrow = 3, align= "hv")




####--------- 8. MaxSSS Threshold & Binary Classification ---------####

# Variables used for threshold evaluation
pred_vars <- c("slp_mean","slp_CV","ws_area",
               "forest","weir_dens",
               "p_bio05","p_bio12","p_bio13")

# Extract environmental matrices for presence/background
occ_pred_input <- occ[, pred_vars] %>% as.data.frame()
bg_pred_input  <- bg[,  pred_vars] %>% as.data.frame()


# (1) Compute MaxSSS threshold for each replicate model
thr_list <- sapply(seq_along(mod_maxent@models), function(i) {
  
  # replicate model
  m_i <- mod_maxent@models[[i]]
  
  # Predict occurrence/background suitability values
  occ_i <- predict(m_i, occ_pred_input,
                   args = "outputformat=cloglog")
  bg_i  <- predict(m_i, bg_pred_input,
                   args = "outputformat=cloglog")
  
  # Combine predictions and presence/absence vector
  Pred_i <- c(occ_i, bg_i)
  Sp_i   <- c(rep(1, length(occ_i)),
              rep(0, length(bg_i)))
  
  # Return MaxSSS threshold 
  ecospat::ecospat.max.tss(Pred_i, Sp_i)$max.threshold
})

# (2) Use the mean of replicate thresholds as final MaxSSS
maxsss_thr <- mean(thr_list, na.rm = TRUE)
maxsss_thr_sd   <- sd(thr_list, na.rm = TRUE)

maxsss_thr
maxsss_thr_sd


# Binary classification using MaxSSS (1=Suitable / 0=Unsuitable)
streamnetwork_pred_binary <- streamnetwork_pred_present %>%
  mutate(suitability = ifelse(probability >= maxsss_thr, 1, 0))

# Visualization of historical binary prediction
present_binary_fig <- ggplot() +
  # (1) Basin boundary
  geom_sf(data = han_watershed, colour = "gray30", fill = NA, linewidth = 0.8) +
  
  # (2) Binary suitability (0 = unsuitable, 1 = suitable)
  geom_sf(
    data = streamnetwork_pred_binary,
    aes(color = as.factor(suitability), linewidth = st_order),
    show.legend = c(color = TRUE, linewidth = FALSE)  
  ) +
  
  # (3) Two-category color scheme
  scale_color_manual(
    values = c("0" = "#88deeb", "1" = "#d7191c"),
    name = "Habitat suitability",
    labels = c("Unsuitable", "Suitable"),
    drop = FALSE,
    na.value = "gray90"
  ) +
  
  # (4) Stream line scaling (no legend for st_order)
  scale_linewidth_continuous(range = c(0.3, 1.2), guide = "none") +
  
  # (5) Coordinate system
  coord_sf(expand = FALSE) +
  
  # (6) Theme and formatting
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    axis.title = element_blank(),          
    panel.grid = element_blank(),           
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA)
  )

present_binary_fig




# Future suitability binary map plotting

plot_future_binary <- function(stream_data, threshold, title_text) {
  stream_data <- stream_data %>%
    mutate(suitability = ifelse(probability >= threshold, 1, 0))

  ggplot() +
    # (1) Basin boundary
    geom_sf(data = han_watershed, colour = "gray30", fill = NA, linewidth = 0.8) +

    # (2) Binary suitability
    geom_sf(
      data = stream_data,
      aes(color = as.factor(suitability), linewidth = st_order),
      show.legend = FALSE
    ) +

    # (3) Color scheme and NA
    scale_color_manual(
      values = c("0" = "#88deeb", "1" = "#d7191c"),
      name = "Habitat suitability",
      labels = c("Unsuitable", "Suitable"),
      drop = FALSE,
      na.value = "gray90"
    ) +

    # (4) Stream line scaling
    scale_linewidth_continuous(range = c(0.3, 1.2), guide = "none") +

    # (5) Coordinate system
    coord_sf(expand = FALSE) +

    # (6) Theme identical to hist_pred_binary_fig
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_blank(),
      plot.caption = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA)
    )
}


####--------- 10. SSP Scenario & Period Mapping ---------####

# Maps SSP to short code used in variable naming
ssp_code_map <- c(SSP126 = "s1", SSP245 = "s2", SSP585 = "s5")

# Display/labeling format for figures/tables
ssp_name_map <- c(SSP126 = "SSP1", SSP245 = "SSP2", SSP585 = "SSP5")

# Future scenario time slices for projection
future_scenarios <- list(
  SSP126 = c("21_40", "41_60", "61_80", "81_100"),
  SSP245 = c("21_40", "41_60", "61_80", "81_100"),
  SSP585 = c("21_40", "41_60", "61_80", "81_100")
)

# Year labels for visualization output
period_year_map <- c(
  "21_40" = "2021–2040",
  "41_60" = "2041–2060",
  "61_80" = "2061–2080",
  "81_100" = "2081–2100"
)


####--------- 11. Binary Prediction Loop ---------####

for (sc in names(future_scenarios)) {
  prefix <- ssp_code_map[sc]   # s1, s2, s5
  ssp_name <- ssp_name_map[sc] # SSP1, SSP2, SSP5
  
  for (period in future_scenarios[[sc]]) {
    stream_var <- paste0("streamnetwork_pred_", prefix, "_", period)
    year_text  <- period_year_map[period]
    title_text <- paste0(ssp_name, " ", year_text, " Prediction (Binary)")
    binary_df <- get(stream_var) %>%
      mutate(suitability = ifelse(probability >= maxsss_thr, 1, 0))
    assign(paste0("streamnetwork_pred_", prefix, "_", period, "_binary"), binary_df)
    assign(paste0(ssp_name, "_", period, "_binary_fig"),
           plot_future_binary(
             stream_data = binary_df,
             threshold  = maxsss_thr,
             title_text = title_text
           ))
  }
}

# Display figures 
SSP1_21_40_binary_fig
SSP1_41_60_binary_fig
SSP1_61_80_binary_fig
SSP1_81_100_binary_fig

SSP2_21_40_binary_fig
SSP2_41_60_binary_fig
SSP2_61_80_binary_fig
SSP2_81_100_binary_fig


SSP5_21_40_binary_fig
SSP5_41_60_binary_fig
SSP5_61_80_binary_fig
SSP5_81_100_binary_fig


library(ggpubr)

# 4x4 grid: Historical + 12 future scenarios
ggarrange(
  SSP1_21_40_binary_fig, SSP1_41_60_binary_fig, SSP1_61_80_binary_fig, SSP1_81_100_binary_fig,
  SSP2_21_40_binary_fig, SSP2_41_60_binary_fig, SSP2_61_80_binary_fig, SSP2_81_100_binary_fig,
  SSP5_21_40_binary_fig, SSP5_41_60_binary_fig, SSP5_61_80_binary_fig, SSP5_81_100_binary_fig,
  ncol = 4, nrow = 3,
  align = "hv",
  common.legend = TRUE,
  legend = "right"
)


####--------- 12. Export Binary Suitability Results in Long Format ---------####


# Function to convert binary prediction object to long-format dataframe
# 1  = suitable within occurrence SBSNCDs
# 0  = unsuitable within occurrence SBSNCDs
# NA = outside occurrence SBSNCDs; not projected
make_binary_long <- function(data, species_name, ssp_scenario_name, time_name) {
  data %>%
    st_drop_geometry() %>%
    transmute(
      species = species_name,
      str_ID = str_ID,
      ssp_scenario = ssp_scenario_name,
      time = time_name,
      occurrence = suitability
    )
}

# Present result
binary_present <- make_binary_long(
  data = streamnetwork_pred_binary,
  species_name = "R.steindachneri",
  ssp_scenario_name = "present",
  time_name = "present"
)

# Future results
binary_future <- bind_rows(
  make_binary_long(streamnetwork_pred_s1_21_40_binary, "R.steindachneri", "SSP126", "2021-2040"),
  make_binary_long(streamnetwork_pred_s1_41_60_binary, "R.steindachneri", "SSP126", "2041-2060"),
  make_binary_long(streamnetwork_pred_s1_61_80_binary, "R.steindachneri", "SSP126", "2061-2080"),
  make_binary_long(streamnetwork_pred_s1_81_100_binary, "R.steindachneri", "SSP126", "2081-2100"),
  
  make_binary_long(streamnetwork_pred_s2_21_40_binary, "R.steindachneri", "SSP245", "2021-2040"),
  make_binary_long(streamnetwork_pred_s2_41_60_binary, "R.steindachneri", "SSP245", "2041-2060"),
  make_binary_long(streamnetwork_pred_s2_61_80_binary, "R.steindachneri", "SSP245", "2061-2080"),
  make_binary_long(streamnetwork_pred_s2_81_100_binary, "R.steindachneri", "SSP245", "2081-2100"),
  
  make_binary_long(streamnetwork_pred_s5_21_40_binary, "R.steindachneri", "SSP585", "2021-2040"),
  make_binary_long(streamnetwork_pred_s5_41_60_binary, "R.steindachneri", "SSP585", "2041-2060"),
  make_binary_long(streamnetwork_pred_s5_61_80_binary, "R.steindachneri", "SSP585", "2061-2080"),
  make_binary_long(streamnetwork_pred_s5_81_100_binary, "R.steindachneri", "SSP585", "2081-2100")
)

# Combine present + future
R_steindachneri_binary_suitability_long_result <- bind_rows(
  binary_present,
  binary_future
)

# Check result
R_steindachneri_binary_suitability_long_result %>%
  count(ssp_scenario, time, occurrence, drop = FALSE)

# Export CSV
write.csv(
  R_steindachneri_binary_suitability_long_result,
  "./data_outcome/R.steindachneri_binary_suitability_long_result.csv",
  row.names = FALSE,
  na = "NA"
)


####--------- 13. TSS, Sensitivity, and Specificity Calculation using PresenceAbsence ---------####

library(PresenceAbsence)

n_models <- length(mod_maxent@models)

tss_list <- numeric(n_models)
sensitivity_list <- numeric(n_models)
specificity_list <- numeric(n_models)

for (i in seq_len(n_models)) {
  
  m_i <- mod_maxent@models[[i]]
  
  occ_i <- predict(
    m_i,
    occ_pred_input,
    args = "outputformat=cloglog"
  )
  
  bg_i <- predict(
    m_i,
    bg_pred_input,
    args = "outputformat=cloglog"
  )
  
  df_eval <- data.frame(
    id = 1:(length(occ_i) + length(bg_i)),
    observed = c(
      rep(1, length(occ_i)),
      rep(0, length(bg_i))
    ),
    predicted = c(occ_i, bg_i)
  )
  
  thr_i <- thr_list[i]
  
  cm <- PresenceAbsence::cmx(
    df_eval,
    threshold = thr_i
  )
  
  sens_val <- PresenceAbsence::sensitivity(cm)$sensitivity
  spec_val <- PresenceAbsence::specificity(cm)$specificity
  
  sensitivity_list[i] <- sens_val
  specificity_list[i] <- spec_val
  tss_list[i] <- sens_val + spec_val - 1
}

mean_tss <- mean(tss_list, na.rm = TRUE)
sd_tss <- sd(tss_list, na.rm = TRUE)

mean_sensitivity <- mean(sensitivity_list, na.rm = TRUE)
sd_sensitivity <- sd(sensitivity_list, na.rm = TRUE)

mean_specificity <- mean(specificity_list, na.rm = TRUE)
sd_specificity <- sd(specificity_list, na.rm = TRUE)

mean_tss
sd_tss

mean_sensitivity
sd_sensitivity

mean_specificity
sd_specificity


# Check replicate-level evaluation metrics
threshold_eval_table <- tibble(
  replicate = seq_len(n_models),
  maxsss_threshold = thr_list,
  sensitivity = sensitivity_list,
  specificity = specificity_list,
  TSS = tss_list
)

threshold_eval_table


# Summary table for reporting
threshold_eval_summary <- tibble(
  metric = c(
    "MaxSSS threshold",
    "Sensitivity at MaxSSS threshold",
    "Specificity at MaxSSS threshold",
    "TSS at MaxSSS threshold"
  ),
  mean = c(
    mean(maxsss_thr, na.rm = TRUE),
    mean_sensitivity,
    mean_specificity,
    mean_tss
  ),
  sd = c(
    maxsss_thr_sd,
    sd_sensitivity,
    sd_specificity,
    sd_tss
  )
)

threshold_eval_summary


####--------- 14. AUC Calculation (10 replicates) ---------####

results_csv <- read.csv(
  "./data_outcome/R.steindachneri_maxent_result/maxentResults.csv"
)

auc_list <- results_csv$Test.AUC

mean_auc <- mean(auc_list, na.rm = TRUE)
sd_auc   <- sd(auc_list, na.rm = TRUE)

mean_auc
sd_auc


####--------- 15. Model selection appendix (ΔAICc) ---------####
library(dplyr)
library(ggplot2)

# Extract tuning results from ENMeval output
tuning_df <- mod_res@results %>%
  dplyr::select(fc, rm, delta.AICc) %>%
  mutate(
    rm = as.numeric(rm),
    fc = factor(fc, levels = c("L", "H", "LQ", "LQH", "LQHP", "LQHPT"))
  )

# Define feature class color palette
fc_colors <- c(
  "L"     = "#F8746B",
  "H"     = "#B79F00",
  "LQ"    = "#00B72E",
  "LQH"   = "#00BCC2",
  "LQHP"  = "#5F9AFF",
  "LQHPT" = "#F562E3"
)

# Plot ΔAICc across regularization multipliers
p_model_selection <- ggplot(tuning_df, aes(rm, delta.AICc, color = fc, group = fc)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_color_manual(values = fc_colors) +
  scale_x_continuous(breaks = sort(unique(tuning_df$rm))) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  labs(
    x = "Regularization Multiplier",
    y = expression(Delta*AIC[c]),
    color = "Feature Class"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = c(-5, max(tuning_df$delta.AICc)))

p_model_selection


####--------- 16. Occurrence and background segment visualization ---------####

# Full stream network for visualization.
# p_bio01 NA segments are NOT removed here.
streamnetwork_vis <- st_read("./data_raw/stream_network_Rhynchocypris.gpkg") %>%
  filter(BBSNCD %in% c("10", "13")) %>%
  st_transform(4326) %>%
  mutate(
    str_ID = as.character(str_ID),
    st_order = as.numeric(st_order)
  ) %>%
  distinct(str_ID, .keep_all = TRUE)


# SBSNCD boundary.
sbsn_watershed <- st_read("./data_raw/WKMSBSN.gpkg") %>%
  st_transform(4326) %>%
  filter(as.character(BBSNCD) %in% c("10", "13"))


# Extract occurrence segments.
occ_seg <- streamnetwork_vis %>%
  filter(str_ID %in% presencedata$str_ID)


# Extract background segments.
bg_seg <- streamnetwork_vis %>%
  filter(str_ID %in% bgdata$str_ID)


# Occurrence segment plot.
occ_segment_plot_clean <- ggplot() +
  
  # Full stream network.
  geom_sf(
    data = streamnetwork_vis,
    aes(linewidth = st_order),
    colour = "gray68",
    alpha = 0.45
  ) +
  
  # SBSNCD boundary.
  geom_sf(
    data = sbsn_watershed,
    fill = NA,
    colour = "#0050be",
    linewidth = 0.75,
    alpha = 0.55
  ) +
  
  # BBSNCD boundary.
  geom_sf(
    data = han_watershed,
    fill = NA,
    colour = "black",
    linewidth = 1
  ) +
  
  # Occurrence segments.
  geom_sf(
    data = occ_seg,
    colour = "red",
    linewidth = 0.64,
    alpha = 0.95,
    lineend = "round"
  ) +
  
  scale_linewidth_continuous(
    range = c(0.6, 1.15),
    guide = "none"
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.caption = element_text(size = 9, colour = "gray30"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "none"
  )


# Background segment plot.
bg_segment_plot_clean <- ggplot() +
  
  # Full stream network.
  geom_sf(
    data = streamnetwork_vis,
    aes(linewidth = st_order),
    colour = "gray68",
    alpha = 0.45
  ) +
  
  
  # SBSNCD boundary.
  geom_sf(
    data = sbsn_watershed,
    fill = NA,
    colour = "#0050be",
    linewidth = 0.8,
    alpha = 0.55
  ) +
  
  # BBSNCD boundary.
  geom_sf(
    data = han_watershed,
    fill = NA,
    colour = "black",
    linewidth = 1.2
  ) +
  
  
  # Background segments.
  geom_sf(
    data = bg_seg,
    colour = "forestgreen",
    linewidth = 0.64,
    alpha = 0.95,
    lineend = "round"
  ) +
  
  scale_linewidth_continuous(
    range = c(0.6, 1.15),
    guide = "none"
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.caption = element_text(size = 9, colour = "gray30"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "none"
  )


# Display plots.
occ_segment_plot_clean
bg_segment_plot_clean
