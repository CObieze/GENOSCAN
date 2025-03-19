# set directory
setwd("./GENOSCAN/")

## Load required libraries
library(vegan)
library(ampvis2)
library(tidyverse)
library(geomtextpath)
library(RColorBrewer)
library(colorRamps)
library(car)

############################# Bacteria Analysis #############################

# Import bacterial data and metadata
df <- amp_load(otutable = "amplicon-data/16s-asv-table.biom")

# Remove overlapping 2021 samples
df <- amp_filter_samples(df, !Year %in% c("RH_2021", "BS_2021", "R_2021"))

# Ensure multiple metadata columns have been read as numeric variables
num_vars <- c("Zn_ppm_", "Ag_ppm_", "Mo_ppm_", "Fe_ppm_", "NO3_ppm_", "N__", "Ni_ppm_", 
              "Pb_ppm_", "S__", "Sn_ppm_", "As_ppm_", "Ba_ppm_", "C__", "Cd_ppm_", "Co_ppm_", 
              "Cr_ppm_", "Cu_ppm_", "K_rs_cm_", "Mn_ppm_", "NH4_ppm_", "P_ppm_", 
              "Altitude", "pH_CaCl2_")

df$metadata[num_vars] <- lapply(df$metadata[num_vars], as.numeric)

# Reorder columns: keep the first column, then numeric, factor, and character columns
df$metadata <- df$metadata %>%
  select(1, where(is.numeric), where(is.factor), where(is.character))

# Subset by region
fermont <- amp_filter_samples(df, Region == "Fermont")
schefferville <- amp_filter_samples(df, Region == "Schefferville")


### Function to Filter, Rarefy, and Normalize Samples ###
process_samples <- function(data, min_reads, niche, rarefy_val) {
  filtered <- amp_filter_samples(data, minreads = min_reads, Niche == niche)
  amp_subset_samples(filtered, rarefy = rarefy_val, normalise = TRUE, removeAbsents = TRUE)
}

# Process Fermont samples
fermont_roots <- process_samples(fermont, 2000, "Root", 2000)
fermont_soil <- process_samples(fermont, 6000, "Bulk soil", 6330)
fermont_rhizosphere <- process_samples(fermont, 6000, "Rhizosphere", 6000)

# Process Schefferville samples
schefferville_roots <- process_samples(schefferville, 2000, "Root", 2000)
schefferville_soil <- process_samples(schefferville, 6000, "Bulk soil", 6330)
schefferville_rhizosphere <- process_samples(schefferville, 6000, "Rhizosphere", 6000)

## Function to normalize metadata for 16S
normalize_metadata <- function(metadata) {
  metadata[, 2:24] <- decostand(metadata[, 2:24], method = "standardize")
  return(metadata)
}

## Function to normalize metadata for ITS
normalize_metadata_fungi <- function(metadata) {
  metadata[, 2:26] <- decostand(metadata[, 2:26], method = "standardize")
  return(metadata)
}

## Function to perform RDA analysis and plotting for Fermont
perform_rda <- function(data, title) {
  rda_result <- amp_ordinate(
    data,
    type = "rda",
    constrain = c("Ecotype", "Site"),
    transform = "hellinger",
    distmeasure = "bray",
    sample_color_by = "Ecotype",
    sample_shape_by = "Site",
    species_plot = FALSE,
    species_nlabels = 5,
    species_label_taxonomy = "Family",
    species_shape = 20,
    species_point_size = 0,
    envfit_numeric = c(2:11,13,14,17:24),
    envfit_signif_level = 0.01,
    sample_colorframe = FALSE,
    sample_label_size = 6,
    sample_point_size = 3,
    detailed_output = FALSE
  ) +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = 0.010, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15:18,7,8,9,10,13)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none", 
          panel.border = element_rect(color = "black", fill = NA))
  
  print(rda_result)
  return(amp_ordinate(data, type = "rda", constrain = c("Ecotype", "Site"), 
                      transform = "hellinger", distmeasure = "bray", 
                      envfit_numeric = c(2:11,13,14,17:24), detailed_output = TRUE))
}

## Function to perform RDA analysis and plotting for schefferville
perform_rda_sch <- function(data, title) {
  rda_result <- amp_ordinate(
    data,
    type = "rda",
    constrain = c("Ecotype", "Site"),
    transform = "hellinger",
    distmeasure = "bray",
    sample_color_by = "Ecotype",
    sample_shape_by = "Site",
    species_plot = FALSE,
    species_nlabels = 5,
    species_label_taxonomy = "Family",
    species_shape = 20,
    species_point_size = 0,
    envfit_numeric = c(2:11,13,14,17:24),
    envfit_signif_level = 0.01,
    sample_colorframe = FALSE,
    sample_label_size = 6,
    sample_point_size = 3,
    detailed_output = FALSE
  ) +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = 0.010, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15:19,7,8,9,10,13)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none", 
          panel.border = element_rect(color = "black", fill = NA))
  
  print(rda_result)
  return(amp_ordinate(data, type = "rda", constrain = c("Ecotype", "Site"), 
                      transform = "hellinger", distmeasure = "bray", 
                      envfit_numeric = c(2:11,13,14,17:24), detailed_output = TRUE))
}

## Function to clean data for multicollinearity
# check for multicolinearity
clean_vif_data <- function(data, vif_threshold = 10) {
  # Ensure it's a data frame
  data <- as.data.frame(data)
  
  # Remove columns with missing values
  data <- data[, colSums(is.na(data)) == 0]
  
  # Create a dummy response variable
  data$dummy_response <- runif(nrow(data))
  
  repeat {
    # Fit linear model
    lm_model <- lm(dummy_response ~ ., data = data)
    
    # Compute VIF values
    vif_values <- vif(lm_model)
    
    # Remove dummy response from VIF check
    vif_values <- vif_values[names(vif_values) != "dummy_response"]
    
    # Find the highest VIF
    max_vif <- max(vif_values)
    
    # Stop if all VIF values are below threshold
    if (max_vif < vif_threshold) {
      break
    }
    
    # Identify the variable with highest VIF
    var_to_remove <- names(which.max(vif_values))
    message(paste("Removing", var_to_remove, "with VIF =", round(max_vif, 2)))
    
    # Drop the highest VIF variable
    data <- select(data, -all_of(var_to_remove))
  }
  
  # Remove dummy response before returning
  data <- select(data, -dummy_response)
  
  return(data)
}

## Function to perform dbRDA analysis
perform_dbrda <- function(hel_data, env_data) {
  dbrda_result <- dbrda(hel_data ~ ., data = env_data, distance = "bray", na.action = na.omit)
  print(dbrda_result)
  print(RsquareAdj(dbrda_result))
  print(anova.cca(dbrda_result, step = 1000))
  print(anova.cca(dbrda_result, step = 1000, by = "term"))
}

## Process bulk soil data (Fermont)
fermont_soil$metadata <- normalize_metadata(fermont_soil$metadata)
fermont_soil <- filter_otus(fermont_soil, filter_otus = 1.0)
stats_BS_rda_fer <- perform_rda(fermont_soil, "Bulk Soil RDA")
stats_BS_rda_fer$model
RsquareAdj(stats_BS_rda_fer$model)
anova.cca(stats_BS_rda_fer$model, by = "term")
BSFermont.hel <- decostand(t(fermont_soil$abund), method = "hellinger")
env.BSFermont <- clean_vif_data(fermont_soil$metadata[, c(2:11,13,14,17:24)])
perform_dbrda(BSFermont.hel, env.BSFermont)

## Process rhizosphere data (Fermont)
fermont_rhizosphere$metadata <- normalize_metadata(fermont_rhizosphere$metadata)
fermont_rhizosphere <- filter_otus(fermont_rhizosphere, filter_otus = 1.0)
stats_RH_rda_fer <- perform_rda(fermont_rhizosphere, "Rhizosphere RDA")
stats_RH_rda_fer$model
RsquareAdj(stats_RH_rda_fer$model)
anova.cca(stats_RH_rda_fer$model, by = "term")
RHFermont.hel <- decostand(t(fermont_rhizosphere$abund), method = "hellinger")
env.RHFermont_clean <- clean_vif_data(fermont_rhizosphere$metadata[, c(2:14, 17:24)])
perform_dbrda(RHFermont.hel, env.RHFermont_clean)

## Process root data (Fermont)
fermont_roots$metadata <- normalize_metadata(fermont_roots$metadata)
fermont_roots <- filter_otus(fermont_roots, filter_otus = 1.0)
stats_roots_rda_fer <- perform_rda(fermont_roots, "Roots RDA")
stats_roots_rda_fer$model
RsquareAdj(stats_roots_rda_fer$model)
anova.cca(stats_roots_rda_fer$model, by = "term")
RootsFermont.hel <- decostand(t(fermont_roots$abund), method = "hellinger")
env.rootsFermont <- clean_vif_data(fermont_roots$metadata[, c(2:14, 17:24)])
perform_dbrda(RootsFermont.hel, env.rootsFermont)

## Process bulk soil data (Schefferville)
schefferville_soil$metadata <- normalize_metadata(schefferville_soil$metadata)
schefferville_soil <- filter_otus(schefferville_soil, filter_otus = 1.0)
stats_BS_rda_sch <- perform_rda_sch(schefferville_soil, "Bulk Soil RDA")
stats_BS_rda_sch$model
RsquareAdj(stats_BS_rda_sch$model)
anova.cca(stats_BS_rda_sch$model, by = "term")
BSSchefferville.hel <- decostand(t(schefferville_soil$abund), method = "hellinger")
env.BSSchefferville <- clean_vif_data(schefferville_soil$metadata[, c(2:11,13,14,17:24)])
perform_dbrda(BSSchefferville.hel, env.BSSchefferville)

## Process rhizosphere data (Schefferville)
schefferville_rhizosphere$metadata <- normalize_metadata(schefferville_rhizosphere$metadata)
schefferville_rhizosphere <- filter_otus(schefferville_rhizosphere, filter_otus = 1.0)
stats_RH_rda_sch <- perform_rda_sch(schefferville_rhizosphere, "Rhizosphere RDA")
stats_RH_rda_sch$model
RsquareAdj(stats_RH_rda_sch$model)
anova.cca(stats_RH_rda_sch$model, by = "term")
RHSchefferville.hel <- decostand(t(schefferville_rhizosphere$abund), method = "hellinger")
env.RHSchefferville_clean <- clean_vif_data(schefferville_rhizosphere$metadata[, c(2:14, 17:24)])
perform_dbrda(RHSchefferville.hel, env.RHSchefferville_clean)

## Process root data (Schefferville)
schefferville_roots$metadata <- normalize_metadata(schefferville_roots$metadata)
schefferville_roots <- filter_otus(schefferville_roots, filter_otus = 1.0)
stats_roots_rda_sch <- perform_rda_sch(schefferville_roots, "Roots RDA")
stats_roots_rda_sch$model
RsquareAdj(stats_roots_rda_sch$model)
anova.cca(stats_roots_rda_sch$model, by = "term")
RootsSchefferville.hel <- decostand(t(schefferville_roots$abund), method = "hellinger")
env.rootsSchefferville <- clean_vif_data(schefferville_roots$metadata[, c(2:14, 17:24)])
perform_dbrda(RootsSchefferville.hel, env.rootsSchefferville)



############################# Fungal Analysis #############################

#import fungal biom and metadata
df_fungi <- amp_load(otutable = "amplicon-data//ITS-asv-table-all-runs-no-taxa.txt",
                     taxonomy = "amplicon-data//ITS-taxonomy-all-runs.txt",
                     metadata = "amplicon-data//ITS-sample-metadata.txt")

# Remove overlapping 2021 samples
df_fungi <- amp_filter_samples(df_fungi, !Year %in% c("2021"))

# Subset by region
fermont_fungi <- amp_filter_samples(df_fungi, Region == "Fermont")
schefferville_fungi <- amp_filter_samples(df_fungi, Region == "Schefferville")

# Process Fermont samples
fermont_fungi_roots <- process_samples(fermont_fungi, 2000, "Root", 2000)
fermont_fungi_soil <- process_samples(fermont_fungi, 6000, "Bulk soil", 6000)
fermont_fungi_rhizosphere <- process_samples(fermont_fungi, 6000, "Rhizosphere", 6000)

# Process Schefferville samples
schefferville_fungi_roots <- process_samples(schefferville_fungi, 2000, "Root", 2000)
schefferville_fungi_soil <- process_samples(schefferville_fungi, 6000, "Bulk soil", 6000)
schefferville_fungi_rhizosphere <- process_samples(schefferville_fungi, 6000, "Rhizosphere", 6000)

## Function to normalize metadata for ITS
normalize_metadata_fungi <- function(metadata) {
  metadata[, c(11,14:31,34,35)] <- decostand(metadata[, c(11,14:31,34,35)], method = "standardize")
  return(metadata)
}

## Function to perform RDA analysis and plotting for Fermont
perform_rda_fungi <- function(data, title) {
  rda_result <- amp_ordinate(
    data,
    type = "rda",
    constrain = c("Ecotype", "Site"),
    transform = "hellinger",
    distmeasure = "bray",
    sample_color_by = "Ecotype",
    sample_shape_by = "Site",
    species_plot = FALSE,
    species_nlabels = 5,
    species_label_taxonomy = "Family",
    species_shape = 20,
    species_point_size = 0,
    envfit_numeric = c(11,14:31,35,35),
    envfit_signif_level = 0.01,
    sample_colorframe = FALSE,
    sample_label_size = 6,
    sample_point_size = 3,
    detailed_output = FALSE
  ) +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = 0.010, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15:18,7,8,9,10,13)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none", 
          panel.border = element_rect(color = "black", fill = NA))
  
  print(rda_result)
  return(amp_ordinate(data, type = "rda", constrain = c("Ecotype", "Site"), 
                      transform = "hellinger", distmeasure = "bray", 
                      envfit_numeric = c(11,14:31,35,35), detailed_output = TRUE))
}

## Function to perform RDA analysis and plotting for schefferville
perform_rda_sch_fungi <- function(data, title) {
  rda_result <- amp_ordinate(
    data,
    type = "rda",
    constrain = c("Ecotype", "Site"),
    transform = "hellinger",
    distmeasure = "bray",
    sample_color_by = "Ecotype",
    sample_shape_by = "Site",
    species_plot = FALSE,
    species_nlabels = 5,
    species_label_taxonomy = "Family",
    species_shape = 20,
    species_point_size = 0,
    envfit_numeric = c(11,14:31,34,35),
    envfit_signif_level = 0.01,
    sample_colorframe = FALSE,
    sample_label_size = 6,
    sample_point_size = 3,
    detailed_output = FALSE
  ) +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = 0.010, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15:19,7,8,9,10,13)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none", 
          panel.border = element_rect(color = "black", fill = NA))
  
  print(rda_result)
  return(amp_ordinate(data, type = "rda", constrain = c("Ecotype", "Site"), 
                      transform = "hellinger", distmeasure = "bray", 
                      envfit_numeric = c(11,14:31,34,35), detailed_output = TRUE))
}

## Process bulk soil data (Fermont)
fermont_fungi_soil$metadata <- normalize_metadata_fungi(fermont_fungi_soil$metadata)
fermont_fungi_soil <- filter_otus(fermont_fungi_soil, filter_otus = 1.0)
stats_BS_rda_fer_fungi <- perform_rda_fungi(fermont_fungi_soil, "Bulk Soil RDA")
stats_BS_rda_fer_fungi$model
RsquareAdj(stats_BS_rda_fer_fungi$model)
anova.cca(stats_BS_rda_fer_fungi$model, by = "term")
BSFermont.hel <- decostand(t(fermont_fungi_soil$abund), method = "hellinger")
env.BSFermont <- clean_vif_data(fermont_fungi_soil$metadata[, c(11,14:31,34,35)])
perform_dbrda(BSFermont.hel, env.BSFermont)

## Process rhizosphere data (Fermont)
fermont_fungi_rhizosphere$metadata <- normalize_metadata_fungi(fermont_fungi_rhizosphere$metadata)
fermont_fungi_rhizosphere <- filter_otus(fermont_fungi_rhizosphere, filter_otus = 1.0)
stats_RH_rda_fer_fungi <- perform_rda_fungi(fermont_fungi_rhizosphere, "Rhizosphere RDA")
stats_RH_rda_fer_fungi$model
RsquareAdj(stats_RH_rda_fer_fungi$model)
anova.cca(stats_RH_rda_fer_fungi$model, by = "term")
RHFermont.hel <- decostand(t(fermont_fungi_rhizosphere$abund), method = "hellinger")
env.RHFermont_clean <- clean_vif_data(fermont_fungi_rhizosphere$metadata[, c(11,14:31,34,35)])
perform_dbrda(RHFermont.hel, env.RHFermont_clean)

## Process root data (Fermont)
fermont_fungi_roots$metadata <- normalize_metadata_fungi(fermont_fungi_roots$metadata)
fermont_fungi_roots <- filter_otus(fermont_fungi_roots, filter_otus = 1.0)
stats_roots_rda_fer_fungi <- perform_rda_fungi(fermont_fungi_roots, "Roots RDA")
stats_roots_rda_fer_fungi$model
RsquareAdj(stats_roots_rda_fer_fungi$model)
anova.cca(stats_roots_rda_fer_fungi$model, by = "term")
RootsFermont.hel <- decostand(t(fermont_fungi_roots$abund), method = "hellinger")
env.rootsFermont <- clean_vif_data(fermont_fungi_roots$metadata[, c(11,14:31,34,35)])
perform_dbrda(RootsFermont.hel, env.rootsFermont)

## Process bulk soil data (Schefferville)
schefferville_fungi_soil$metadata <- normalize_metadata_fungi(schefferville_fungi_soil$metadata)
schefferville_fungi_soil <- filter_otus(schefferville_fungi_soil, filter_otus = 1.0)
stats_BS_rda_sch_fungi <- perform_rda_sch_fungi(schefferville_fungi_soil, "Bulk Soil RDA")
stats_BS_rda_sch_fungi$model
RsquareAdj(stats_BS_rda_sch_fungi$model)
anova.cca(stats_BS_rda_sch_fungi$model, by = "term")
BSSchefferville.hel <- decostand(t(schefferville_fungi_soil$abund), method = "hellinger")
env.BSSchefferville <- clean_vif_data(schefferville_fungi_soil$metadata[, c(11,14:31,34,35)])
perform_dbrda(BSSchefferville.hel, env.BSSchefferville)

## Process rhizosphere data (Schefferville)
schefferville_fungi_rhizosphere$metadata <- normalize_metadata_fungi(schefferville_fungi_rhizosphere$metadata)
schefferville_fungi_rhizosphere <- filter_otus(schefferville_fungi_rhizosphere, filter_otus = 1.0)
stats_RH_rda_sch_fungi <- perform_rda_sch_fungi(schefferville_fungi_rhizosphere, "Rhizosphere RDA")
stats_RH_rda_sch_fungi$model
RsquareAdj(stats_RH_rda_sch_fungi$model)
anova.cca(stats_RH_rda_sch_fungi$model, by = "term")
RHSchefferville.hel <- decostand(t(schefferville_fungi_rhizosphere$abund), method = "hellinger")
env.RHSchefferville_clean <- clean_vif_data(schefferville_fungi_rhizosphere$metadata[, c(11,14:31,34,35)])
perform_dbrda(RHSchefferville.hel, env.RHSchefferville_clean)

## Process root data (Schefferville)
schefferville_fungi_roots$metadata <- normalize_metadata_fungi(schefferville_fungi_roots$metadata)
schefferville_fungi_roots <- filter_otus(schefferville_fungi_roots, filter_otus = 1.0)
stats_roots_rda_sch_fungi <- perform_rda_sch_fungi(schefferville_fungi_roots, "Roots RDA")
stats_roots_rda_sch_fungi$model
RsquareAdj(stats_roots_rda_sch_fungi$model)
anova.cca(stats_roots_rda_sch_fungi$model, by = "term")
RootsSchefferville.hel <- decostand(t(schefferville_fungi_roots$abund), method = "hellinger")
env.rootsSchefferville <- clean_vif_data(schefferville_fungi_roots$metadata[, c(11,14:31,34,35)])
perform_dbrda(RootsSchefferville.hel, env.rootsSchefferville)
