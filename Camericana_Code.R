Sys.setenv(LANG = "en")

# ==============================================================================
# PORTABLE WORKING DIRECTORY SETUP
# ==============================================================================
# Automatically sets the working directory to the location of this script when 
# run inside RStudio or Positron. Place this script alongside all input files 
# in the root directory ("Camericana_DataRepository").
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(glmmTMB)
library(ggeffects)
library(sjPlot)

# ==============================================================================
# PART 1: GEOGRAPHIC FRAMEWORK & SPATIAL MASTER BUILD
# ==============================================================================
# Constructing a unified UK spatial grid using electoral boundaries and 
# standardizing human population density across nations.
# ==============================================================================
cat("\n--- CHECKING FOR SPATIAL MASTER ---\n")
spatial_rds <- "Camericana_UK_Spatial_Master.rds"

if (file.exists(spatial_rds)) {
  uk_spatial_master <- readRDS(spatial_rds)
  cat("Loaded existing Camericana_UK_Spatial_Master.rds\n")
} else {
  cat("Building Spatial Master from shapefiles and calculating areas...\n")
  
  # A. England & Wales (MSOA 2021)
  ew_sf <- sf::st_read("/vsizip/Camericana_England&Wales_Shapefile.zip/MSOA_2021_EW_BGC_V3.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>%
    dplyr::mutate(Area_SqKm = as.numeric(sf::st_area(geometry)) / 1000000) %>%
    dplyr::select(Zone_Code = MSOA21CD, Zone_Name = MSOA21NM, Area_SqKm)
  
  ew_pop <- read.csv("Camericana_England&Wales_Population.csv", sep = ";", check.names = FALSE) %>%
    dplyr::mutate(Total_Pop = as.numeric(stringr::str_replace_all(EW_MSOA2021_PopEstimate2022, "[ ,]", ""))) %>%
    dplyr::select(Zone_Code = EW_MSOA2021_Code, Total_Pop)
  
  ew_master <- ew_sf %>% dplyr::left_join(ew_pop, by = "Zone_Code")
  
  # B. Scotland (IZ 2022)
  scot_sf <- sf::st_read("/vsizip/Camericana_Scotland_Shapefile.zip/SG_IntermediateZoneBdry_2022_MHW.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>%
    dplyr::mutate(Area_SqKm = as.numeric(sf::st_area(geometry)) / 1000000) %>%
    dplyr::select(Zone_Code = IZCode, Zone_Name = IZName, Area_SqKm)
  
  scot_pop <- read.csv("Camericana_Scotland_Population.csv", sep = ";", check.names = FALSE) %>%
    dplyr::mutate(Total_Pop = as.numeric(stringr::str_replace_all(S_IZ2022_PopEstimate2022, "[ ,]", ""))) %>%
    dplyr::select(Zone_Code = S_IZ2022_Code, Total_Pop)
  
  scot_master <- scot_sf %>% dplyr::left_join(scot_pop, by = "Zone_Code")
  
  # C. Northern Ireland (SDZ 2021)
  ni_sf <- sf::st_read("/vsizip/Camericana_NorthernIreland_Shapefile.zip/SDZ2021.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>%
    dplyr::mutate(Area_SqKm = as.numeric(sf::st_area(geometry)) / 1000000) %>%
    dplyr::select(Zone_Code = SDZ2021_cd, Zone_Name = SDZ2021_nm, Area_SqKm)
  
  ni_pop <- read.csv("Camericana_NorthernIreland_Population.csv", sep = ";", check.names = FALSE) %>%
    dplyr::mutate(Total_Pop = as.numeric(stringr::str_replace_all(NI_SDZ2021_PopEstimate2021, "[ ,]", ""))) %>%
    dplyr::select(Zone_Code = NI_SDZ2021_Code, Total_Pop)
  
  ni_master <- ni_sf %>% dplyr::left_join(ni_pop, by = "Zone_Code")
  
  # D. Bind spatial geometries and standardize population density for modeling
  uk_spatial_master <- dplyr::bind_rows(ew_master, scot_master, ni_master) %>%
    dplyr::filter(!is.na(Total_Pop)) %>%
    dplyr::mutate(
      Pop_Density_SqKm = Total_Pop / Area_SqKm,
      Pop_Density_Z = scale(Pop_Density_SqKm)[,1]
    )
    
  saveRDS(uk_spatial_master, spatial_rds)
  cat("Saved new Camericana_UK_Spatial_Master.rds\n")
}

# ==============================================================================
# PART 2: BIOLOGICAL & CLIMATE DATA PROCESSING
# ==============================================================================
# Extracts baseline historical temperatures (1991-2020) and integrates
# presence-only occurrence records for target and calibration species.
# ==============================================================================
cat("\n--- PROCESSING BIOLOGICAL & CLIMATE DATA ---\n")

# Load and average the 12-month temperature raster into a single annual baseline
temp_sf <- sf::st_read("/vsizip/Camericana_Temperature_Shapefile.zip/Monthly_Temperature_Observations_1991-2020.shp", quiet = TRUE) %>%
  sf::st_transform(27700)

grid_sf <- temp_sf %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    baseline_1991_2020 = mean(c_across(starts_with("tas")), na.rm = TRUE) 
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(baseline_1991_2020)

# Export tabulated historical temperature measurements as a CSV for open reference
temp_csv_file <- "Camericana_Temperature_HistoricalMeasurements.csv"
if (!file.exists(temp_csv_file)) {
  # Dynamically identify the ID column (checking common variants like OBJECTID, FID, Id)
  id_col <- intersect(c("OBJECTID", "objectid", "Objectid", "Id", "FID", "fid", "ID"), names(temp_sf))[1]
  
  temp_csv_data <- sf::st_drop_geometry(temp_sf) %>%
    dplyr::mutate(
      baseline_1991_2020 = temp_sf %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(val = mean(c_across(starts_with("tas")), na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        dplyr::pull(val)
    )
  
  # Select and rename the ID column to "SquareID", followed by 12 tas columns and the yearly average
  if (!is.na(id_col)) {
    temp_csv_data <- temp_csv_data %>% 
      dplyr::rename(SquareID = !!id_col) %>%
      dplyr::select(SquareID, starts_with("tas"), baseline_1991_2020)
  } else {
    # Fallback if no standard ID is matched: create a sequential SquareID
    temp_csv_data <- temp_csv_data %>% 
      dplyr::mutate(SquareID = dplyr::row_number()) %>%
      dplyr::select(SquareID, starts_with("tas"), baseline_1991_2020)
  }
  
  write.csv(temp_csv_data, temp_csv_file, row.names = FALSE)
  cat("Exported Camericana_Temperature_HistoricalMeasurements.csv with SquareID\n")
}

# Load target (C. americana) and calibration (P. coccinea, C. septempunctata) records
ca_pts <- read.csv("Camericana_Records_Camericana.csv", sep = ";") %>%
  sf::st_as_sf(coords = c("Ca_Lon", "Ca_Lat"), crs = 4326) %>% sf::st_transform(27700)

pc_pts <- read.csv("Camericana_Records_Pcoccinea.csv", sep = ";") %>%
  sf::st_as_sf(coords = c("Pc_Lon", "Pc_Lat"), crs = 4326) %>% sf::st_transform(27700)

c7_pts <- read.csv("Camericana_Records_Cseptempunctata.csv", sep = ";") %>%
  sf::st_as_sf(coords = c("C7_Lon", "C7_Lat"), crs = 4326) %>% sf::st_transform(27700)

# ==============================================================================
# CLIMATE EXTRACTION: AREA-WEIGHTED INTERSECTION
# ==============================================================================
cat("\n--- CALCULATING BLENDED TEMPERATURES ---\n")
climate_rds <- "Camericana_Zone_Climate_Blended.rds"

if (file.exists(climate_rds)) {
  zone_climate <- readRDS(climate_rds)
  cat("Loaded existing Camericana_Zone_Climate_Blended.rds (Skipped intersection!)\n")
} else {
  cat("Calculating area-weighted temperatures...\n")
  
  valid_uk_spatial <- sf::st_make_valid(dplyr::select(uk_spatial_master, Zone_Code))
  valid_grid <- sf::st_make_valid(dplyr::select(grid_sf, baseline_1991_2020))
  
  zone_climate_intersect <- sf::st_intersection(valid_uk_spatial, valid_grid)
  
  # Calculate exact temperature based on proportional overlap of climate grid
  zone_climate <- zone_climate_intersect %>%
    dplyr::mutate(intersect_area = as.numeric(sf::st_area(geometry))) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(Zone_Code) %>%
    dplyr::summarise(
      Temp_Hist = weighted.mean(baseline_1991_2020, w = intersect_area, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(Temp_Hist)) %>%
    dplyr::mutate(Temp_Hist_Z = scale(Temp_Hist)[,1])
  
  saveRDS(zone_climate, climate_rds)
  rm(valid_uk_spatial, valid_grid, zone_climate_intersect)
  gc()
}

# Aggregate biological records to census area bounds (1994-2021)
ca_zone <- sf::st_join(ca_pts, uk_spatial_master, join = sf::st_intersects) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = Ca_Year) %>%
  dplyr::filter(Year >= 1994 & Year <= 2021 & !is.na(Zone_Code)) %>%
  dplyr::group_by(Zone_Code, Year) %>%
  dplyr::summarise(Count_Ca = dplyr::n(), .groups = "drop")

pc_zone <- sf::st_join(pc_pts, uk_spatial_master, join = sf::st_intersects) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = Pc_Year) %>%
  dplyr::filter(Year >= 1994 & Year <= 2021 & !is.na(Zone_Code)) %>%
  dplyr::group_by(Zone_Code, Year) %>%
  dplyr::summarise(Count_Pc = dplyr::n(), .groups = "drop")

c7_zone <- sf::st_join(c7_pts, uk_spatial_master, join = sf::st_intersects) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = C7_Year) %>%
  dplyr::filter(Year >= 1994 & Year <= 2021 & !is.na(Zone_Code)) %>%
  dplyr::group_by(Zone_Code, Year) %>%
  dplyr::summarise(Count_C7 = dplyr::n(), .groups = "drop")

# ==============================================================================
# PART 3: MASTER PANEL & GIS SUMMARIES 
# ==============================================================================
# Constructs the zero-inflated longitudinal panel for statistical modelling 
# and the flattened summary layer for GIS visualization.
# ==============================================================================
cat("\n--- BUILDING MASTER PANELS & SUMMARIES ---\n")

active_zones <- unique(c(ca_zone$Zone_Code, pc_zone$Zone_Code, c7_zone$Zone_Code))

# Build longitudinal panel integrating all covariates
master_panel <- tidyr::expand_grid(Zone_Code = active_zones, Year = 1994:2021) %>%
  dplyr::left_join(ca_zone, by = c("Zone_Code", "Year")) %>%
  dplyr::left_join(pc_zone, by = c("Zone_Code", "Year")) %>%
  dplyr::left_join(c7_zone, by = c("Zone_Code", "Year")) %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>%
  dplyr::left_join(sf::st_drop_geometry(uk_spatial_master)[, c("Zone_Code", "Pop_Density_SqKm", "Pop_Density_Z")], by = "Zone_Code") %>%
  dplyr::mutate(
    Count_Ca = ifelse(is.na(Count_Ca), 0, Count_Ca),
    Count_Pc = ifelse(is.na(Count_Pc), 0, Count_Pc),
    Count_C7 = ifelse(is.na(Count_C7), 0, Count_C7),
    Presence_Ca = ifelse(Count_Ca > 0, 1, 0),
    Presence_Pc = ifelse(Count_Pc > 0, 1, 0),
    Presence_C7 = ifelse(Count_C7 > 0, 1, 0),
    Year_Centered = Year - 1994
  ) %>%
  dplyr::filter(!is.na(Temp_Hist_Z) & !is.na(Pop_Density_Z))

# Sub-panel strictly isolated for zero-truncated abundance modelling
abundance_panel <- master_panel %>% dplyr::filter(Count_Ca > 0)

# Cumulative totals for GIS mapping and Density Coefficient calculation
ca_totals <- ca_zone %>% dplyr::group_by(Zone_Code) %>% dplyr::summarise(Total_Ca = sum(Count_Ca), .groups = "drop")
pc_totals <- pc_zone %>% dplyr::group_by(Zone_Code) %>% dplyr::summarise(Total_Pc = sum(Count_Pc), .groups = "drop")
c7_totals <- c7_zone %>% dplyr::group_by(Zone_Code) %>% dplyr::summarise(Total_C7 = sum(Count_C7), .groups = "drop")

zone_gis_summary <- sf::st_drop_geometry(uk_spatial_master) %>%
  dplyr::left_join(ca_totals, by = "Zone_Code") %>%
  dplyr::left_join(pc_totals, by = "Zone_Code") %>%
  dplyr::left_join(c7_totals, by = "Zone_Code") %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>%
  dplyr::mutate(
    Total_Ca = ifelse(is.na(Total_Ca), 0, Total_Ca),
    Total_Pc = ifelse(is.na(Total_Pc), 0, Total_Pc),
    Total_C7 = ifelse(is.na(Total_C7), 0, Total_C7),
    Density_Coefficient = Total_Ca / Total_Pop * 1000
  ) %>%
  dplyr::select(Zone_Code, Total_Ca, Total_Pc, Total_C7, Density_Coefficient, Temp_Hist, Temp_Hist_Z)

rm(ew_master, ew_pop, ew_sf, ni_master, ni_pop, ni_sf, scot_master, scot_pop, scot_sf, 
   ca_pts, pc_pts, c7_pts, temp_sf, grid_sf)
gc() 

# ==============================================================================
# SUMMARY STATISTICS FOR MANUSCRIPT TEXT
# ==============================================================================
cat("\n--- RAW STANDARD DEVIATIONS FOR TEXT --- \n")

sd_temp <- sd(master_panel$Temp_Hist, na.rm = TRUE)
sd_pop  <- sd(master_panel$Pop_Density_SqKm, na.rm = TRUE)

cat(sprintf("1 Standard Deviation of Temperature: %.3f °C\n", sd_temp))
cat(sprintf("1 Standard Deviation of Population Density: %.1f people per SqKm\n", sd_pop))

# ==============================================================================
# PART 4: STATISTICAL ANALYSIS (H2: PRESENCE MODEL)
# ==============================================================================
cat("\n--- CALCULATING INITIAL ESTABLISHMENT THRESHOLDS (B-GLM) ---\n")
presence_rds <- "Camericana_Presence_Model.rds"

if (file.exists(presence_rds)) {
  Camericana_Presence_Model <- readRDS(presence_rds)
  cat("Loaded existing Camericana_Presence_Model.rds\n")
} else {
  Camericana_Presence_Model <- glm(
    Presence_Ca ~ 
      Year_Centered +       # Chronological expansion over time
      Pop_Density_Z +       # Primary proxy for urban infrastructure (concrete penalty)
      Temp_Hist_Z +         # Baseline thermal limits for establishment
      Presence_Pc +         # Calibration covariate: Natural woodland observer effort
      Presence_C7,          # Calibration covariate: Domestic garden observer effort
    
    family = binomial(link = "logit"), 
    data = master_panel
  )
  saveRDS(Camericana_Presence_Model, presence_rds)
  cat("Saved new Camericana_Presence_Model.rds\n")
}

# ==============================================================================
# PART 5: STATISTICAL ANALYSIS (H2: ABUNDANCE MODEL)
# ==============================================================================
cat("\n--- CALCULATING DRIVERS OF ABUNDANCE (TNB-GLMM) ---\n")
abundance_rds <- "Camericana_Abundance_Model.rds"

if (file.exists(abundance_rds)) {
  Camericana_Abundance_Model <- readRDS(abundance_rds)
  cat("Loaded existing Camericana_Abundance_Model.rds\n")
} else {
  Camericana_Abundance_Model <- glmmTMB::glmmTMB(
    Count_Ca ~ 
      Year_Centered *       # 3-way interaction assessing the rate of population 
      Pop_Density_Z *       # expansion across varying human densities and 
      Temp_Hist_Z +         # regional thermal carrying capacities
      Count_Pc +            # Calibration covariate: Volume of woodland observers
      Count_C7 +            # Calibration covariate: Volume of garden observers
      (1 | Zone_Code),      # Random intercept handling unmeasured geographic variance
    
    data = abundance_panel, 
    family = glmmTMB::truncated_nbinom2 # Handles zero-truncation and overdispersion
  )
  saveRDS(Camericana_Abundance_Model, abundance_rds)
  cat("Saved new Camericana_Abundance_Model.rds\n")
}

# ==============================================================================
# H2: DENSITY-DEPENDENT ABUNDANCE PLOT (INTERACTIONS)
# ==============================================================================
predictions_rds <- "Camericana_Abundance_Predictions.rds"

if (file.exists(predictions_rds)) {
  Camericana_Abundance_Predictions <- readRDS(predictions_rds)
  cat("Loaded existing Camericana_Abundance_Predictions.rds\n")
} else {
  # Simulate theoretical outcomes holding covariates strictly at zero
  Camericana_Abundance_Predictions <- ggeffects::ggpredict(
    Camericana_Abundance_Model, # The fitted truncated negative binomial GLMM
    terms = c(                  # The variables to cross-tabulate and plot
      "Year_Centered",             # X-axis: Chronological timeline of the invasion
      "Pop_Density_Z [-1, 0, 1]",  # Grouping lines: Rural (-1 SD), Average (0 SD), and Urban (+1 SD) densities
      "Temp_Hist_Z [-1, 0, 1]"     # Facet panels: Cooler (-1 SD), Average (0 SD), and Warmer (+1 SD) climates
    ),
    condition = c(Count_Pc = 0, Count_C7 = 0) # Control: Neutralise observer effort by setting calibration counts to zero
  )
  saveRDS(Camericana_Abundance_Predictions, predictions_rds)
  cat("Saved new Camericana_Abundance_Predictions.rds\n")
}

# Standardize labelling for ecological interpretation
plot_data <- as.data.frame(Camericana_Abundance_Predictions) %>%
  dplyr::mutate(
    facet = factor(facet, 
                   levels = c("-1", "0", "1"), 
                   labels = c("Temperature: Cooler (-1 SD)", 
                              "Temperature: Average (0 SD)", 
                              "Temperature: Warmer (+1 SD)")),
    
    group_name = dplyr::case_when(
      group == "-1" ~ "More Rural (-1 SD)",
      group == "0"  ~ "Average (0 SD)",
      group == "1"  ~ "More Urban (+1 SD)"
    ),
    
    group_name = factor(group_name, levels = c("More Rural (-1 SD)", "Average (0 SD)", "More Urban (+1 SD)"))
  )

cat("\n--- EXTRACTING LOGARITHMIC SLOPE COEFFICIENTS ---\n")
reproducibility_slopes <- plot_data %>%
  dplyr::group_by(facet, group_name) %>%
  dplyr::filter(x == min(x) | x == max(x)) %>%
  dplyr::arrange(x) %>%
  dplyr::summarize(
    Log10_Slope = (log10(dplyr::last(predicted)) - log10(dplyr::first(predicted))) / (dplyr::last(x) - dplyr::first(x)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(facet, group_name)

print(reproducibility_slopes)

# Construct 3-way interaction plot
abundance_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x + 1994, y = predicted, color = group_name, fill = group_name)) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  
  ggplot2::facet_wrap(~ facet) +
  
  ggplot2::scale_color_manual(
    values = c("More Rural (-1 SD)" = "#3182bd", 
               "Average (0 SD)"     = "#756bb1", 
               "More Urban (+1 SD)" = "#de2d26")
  ) +
  ggplot2::scale_fill_manual(
    values = c("More Rural (-1 SD)" = "#3182bd", 
               "Average (0 SD)"     = "#756bb1", 
               "More Urban (+1 SD)" = "#de2d26")
  ) +
  
  ggplot2::scale_x_continuous(breaks = seq(1995, 2020, by = 10)) +
  ggplot2::scale_y_log10() + 
  
  ggplot2::labs(
    title = "Climate and Density-Dependent Expansion of C. americana", 
    subtitle = "Predicted abundance across human population density and historical temperature profiles",
    x = "Year", 
    y = "Predicted Annual Abundance Index (PAAI) (Log10)"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position = "none", 
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    panel.spacing = ggplot2::unit(2, "lines"),
    panel.border = ggplot2::element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

print(abundance_plot)

# ==============================================================================
# PART 6: CLIMATE SIMULATIONS (H3: TEMPERATURE SENSITIVITY)
# ==============================================================================
# Projects future spatial distributions under progressive uniform warming 
# scenarios. This combines the probability of establishment (Presence model) 
# with the expected population size (Abundance model) to map carrying capacity.
# ==============================================================================
cat("\n--- GENERATING TEMPERATURE SENSITIVITY PREDICTIONS ---\n")

# 1. Extract historical baseline statistics to properly scale future scenarios
temp_hist_vector <- sf::st_drop_geometry(uk_spatial_master) %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>% 
  dplyr::pull(Temp_Hist)

temp_mean <- mean(temp_hist_vector, na.rm = TRUE)
temp_sd   <- sd(temp_hist_vector, na.rm = TRUE)

# 2. Create the neutral baseline dataset
base_pred_data <- uk_spatial_master %>%
  sf::st_drop_geometry() %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>%
  dplyr::mutate(
    Year_Centered = 27,   # Locks chronological expansion at year 2021
    Presence_Pc = 0,      # Neutralises woodland observer bias (presence)
    Presence_C7 = 0,      # Neutralises garden observer bias (presence)
    Count_Pc = 0,         # Neutralises woodland observer bias (abundance)
    Count_C7 = 0          # Neutralises garden observer bias (abundance)
  )

future_predictions_table <- base_pred_data %>% dplyr::select(Zone_Code)

# 3. Define the uniform warming increments (in degrees Celsius)
scenarios <- c(0.0, 1.5, 2.0, 2.5, 3.0, 4.0)

# 4. Loop through each scenario, calculate the hurdle model, and store the results
for (temp_inc in scenarios) {
  
  # Format column names dynamically (e.g., Pred_Abund_15C)
  formatted_temp <- sprintf("%.1f", temp_inc) 
  col_name <- ifelse(temp_inc == 0.0, "Pred_Abund_Baseline", paste0("Pred_Abund_", gsub("\\.", "", formatted_temp), "C"))
  
  # Step A: Apply the warming increment to the raw temperature, then re-standardise
  temp_scenario_data <- base_pred_data %>%
    dplyr::mutate(Temp_Hist_Z = ((Temp_Hist + temp_inc) - temp_mean) / temp_sd)
  
  # Step B: Predict the probability of presence (Binomial GLM)
  prob_presence <- predict(
    Camericana_Presence_Model, 
    newdata = temp_scenario_data, 
    type = "response"
  )
  
  # Step C: Predict the conditional abundance (Truncated Negative Binomial GLMM)
  cond_abundance <- predict(
    Camericana_Abundance_Model, 
    newdata = temp_scenario_data, 
    type = "response", 
    re.form = NA 
  )
  
  # Step D: Combine into the final Predicted Annual Abundance Index (PAAI)
  future_predictions_table[[col_name]] <- prob_presence * cond_abundance
}

# ==============================================================================
# PART 7: MASTER UNIFIED SPATIAL EXPORT
# ==============================================================================
cat("\n--- CHECKING GIS EXPORT ---\n")
map_rds <- "Camericana_Future_Projections_Data.rds"
gpkg_file <- "Camericana_Geopackage.gpkg"
gpkg_table <- "Camericana_Geopackage_Table.csv"

if (file.exists(map_rds)) {
  final_master_map <- readRDS(map_rds)
  cat("Loaded existing Camericana_Future_Projections_Data.rds\n")
} else {
  final_master_map <- uk_spatial_master %>%
    dplyr::left_join(zone_gis_summary, by = "Zone_Code") %>%
    dplyr::left_join(future_predictions_table, by = "Zone_Code") %>%
    dplyr::select(Zone_Code, Zone_Name, Total_Pop, Area_SqKm, Pop_Density_SqKm, Pop_Density_Z,
                  Temp_Hist, Temp_Hist_Z, Total_Pc, Total_C7,
                  Total_Ca, Density_Coefficient, 
                  Pred_Abund_Baseline, Pred_Abund_15C, Pred_Abund_20C, Pred_Abund_25C, Pred_Abund_30C, Pred_Abund_40C)
  
  final_master_map$Dummy_Scale <- seq(0, 3.1, length.out = nrow(final_master_map))
  saveRDS(final_master_map, map_rds)
  cat("Saved Camericana_Future_Projections_Data.rds\n")
}

# Ensure Geopackage and CSV table are always written if missing
if (!file.exists(gpkg_file) || !file.exists(gpkg_table)) {
  if (!exists("final_master_map")) final_master_map <- readRDS(map_rds)
  sf::st_write(final_master_map, gpkg_file, layer = "main", append = FALSE)
  write.csv(sf::st_drop_geometry(final_master_map), gpkg_table, row.names = FALSE)
  cat("Success! Combined spatial layer written to Camericana_Geopackage.gpkg (layer 'main') and Camericana_Geopackage_Table.csv.\n")
} else {
  cat("Geopackage and table already exist in directory.\n")
}

# ==============================================================================
# PART 8: EXPORTING PUBLICATION-READY TABLES
# ==============================================================================
cat("\n--- EXPORTING TABLES ---\n")

sjPlot::tab_model(Camericana_Presence_Model, 
                  digits = 3, 
                  file = "Camericana_Table_Presence_Model.html",
                  title = "Presence Model: Binomial GLM")

summary(Camericana_Presence_Model)

sjPlot::tab_model(Camericana_Abundance_Model, 
                  digits = 3,    
                  wrap.labels = 100,             
                  file = "Camericana_Table_Abundance_Model.html",
                  title = "Abundance Model: Truncated Negative Binomial GLMM")

summary(Camericana_Abundance_Model)

cat("\n--- PIPELINE COMPLETE! ---\n")