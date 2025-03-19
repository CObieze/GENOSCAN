library(DHARMa)
library(lmerTest)
library(dplyr)

# Set working directory
setwd("./GENOSCAN/")

# ----------- 16S Analysis -----------
# Read and preprocess data
alpha_div_16s <- read.delim("amplicon-data/16s-alpha-diversity.txt", row.names = 1)
env_16s <- read.delim("amplicon-data/16s-sample-metadata.txt", row.names = 1)
df_16s <- merge(alpha_div_16s, env_16s, by = "row.names") %>%
  filter(!Year %in% c("RH_2021", "BS_2021", "R_2021")) %>%
  mutate(observed_features = as.numeric(observed_features)) %>%
  select(1:4) %>%
  column_to_rownames(var = "Row.names")

# Subset and preprocess function
process_niche_data <- function(niche, alpha_div, env_data) {
  df <- merge(alpha_div, env_data, by = "row.names") %>%
    filter(!Year %in% c("RH_2021", "BS_2021", "R_2021")) %>%
    filter(Niche == niche) %>%
    mutate(across(c(shannon_entropy, observed_features, pielou_evenness), as.numeric))
  return(df)
}

# Define niches
niches <- c("Bulk soil", "Rhizosphere", "Root")

# List of physicochemical parameters
env_parameters <- c("Ag_ppm", "As_ppm", "Ba_ppm", "Cd_ppm", "Co_ppm", 
                    "Cr_ppm", "Cu_ppm", "Fe_ppm", "Mn_ppm", "Mo_ppm", 
                    "Ni_ppm", "P_ppm", "Pb_ppm", "Sn_ppm", "Zn_ppm",
                    "C", "N", "S", "Altitude", "pH_CaCl2")

# List of diversity indices
diversity_indices <- c("shannon_entropy", "observed_features", "pielou_evenness")

# run model with site as random effect
for (niche in niches) {
  df_niche <- process_niche_data(niche, alpha_div_16s, env_16s)
  
  # Standardize environmental parameters
  df_niche[env_parameters] <- scale(df_niche[env_parameters])
  
  # Perform PCA
  pca_env <- prcomp(df_niche[, env_parameters], center = TRUE, scale. = TRUE)
  df_niche$env_parameters_PC1 <- pca_env$x[, 1]
  
  for (index in diversity_indices) {
    df_niche[[paste0("log_", index)]] <- log(df_niche[[index]] + 1)
    
    # Fit the linear mixed-effects model
    model <- lmer(as.formula(paste(index, "~ Ecotype + env_parameters_PC1 + (1 | Site)")),
                  data = df_niche)
    
    # Perform ANOVA
    anova_results <- anova(model)
    print(anova_results)
    
    # Model diagnostics
    simulation_output <- simulateResiduals(fittedModel = model)
    plot(simulation_output)
    
    # Save results
    write.csv(anova_results, paste0(niche, "-lmer-", index, "-bacteria.csv"))
  }
}



# -------------- ITS analysis ---------------
# Read and preprocess data
alpha_div_its <- read.delim("amplicon-data/ITS-alpha-diversity.txt", row.names = 1)
env_its <- read.delim("amplicon-data/ITS-sample-metadata.txt", row.names = 1)
alpha_div_its <- merge(alpha_div_its, env_its, by = "row.names") %>%
  filter(Year != "2021") %>%
  mutate(observed_features = as.numeric(observed_features)) %>%
  select(1:4) %>%
  column_to_rownames(var = "Row.names")

# List of physicochemical parameters
env_parameters_fungi <- c("Ag_ppm", "As_ppm", "Ba_ppm", "Cd_ppm", "Co_ppm", 
                    "Cr_ppm", "Cu_ppm", "Fe_ppm", "Mn_ppm", "Mo_ppm", 
                    "Ni_ppm", "P_ppm", "Pb_ppm", "Sn_ppm", "Zn_ppm",
                    "C", "N", "S", "Altitude", "pH_CaCl2")

# run model with site as random effect
for (niche in niches) {
  df_niche <- process_niche_data(niche, alpha_div_its, env_its)
  
  # Standardize environmental parameters
  df_niche[env_parameters_fungi] <- scale(df_niche[env_parameters_fungi])
  
  # Perform PCA
  pca_env <- prcomp(df_niche[, env_parameters_fungi], center = TRUE, scale. = TRUE)
  df_niche$env_parameters_PC1 <- pca_env$x[, 1]
  
  for (index in diversity_indices) {
    df_niche[[paste0("log_", index)]] <- log(df_niche[[index]] + 1)
    
    # Fit the linear mixed-effects model
    model <- lmer(as.formula(paste(index, "~ Ecotype + env_parameters_PC1 + (1 | Site)")),
                  data = df_niche)
    
    # Perform ANOVA
    anova_results <- anova(model)
    print(anova_results)
    
    # Model diagnostics
    simulation_output <- simulateResiduals(fittedModel = model)
    plot(simulation_output)
    
    # Save results
    write.csv(anova_results, paste0(niche, "-lmer-", index, "-fungi.csv"))
  }
}