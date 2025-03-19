# Load required libraries
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set working directory
setwd("./GENOSCAN/")

# Function for statistical analysis
perform_stats <- function(data, y_var, group_var, filename_prefix) {
  stat_result <- compare_means(as.formula(paste(y_var, "~", group_var)), data = data, method = "wilcox.test")
  write.csv(as.data.frame(stat_result), paste0(filename_prefix, "-", y_var, "-", group_var, ".csv"))
}

# Function to create plots
create_plot <- function(data, y_var, y_label, region, is_bacteria = FALSE) {
  axis_text_x <- if (is_bacteria) element_blank() else {
    if (region == "Schefferville") {
      element_text(face = "plain", colour = "black", size = 8, angle = 60, hjust = 1.01)
    } else {
      element_blank()
    }
  }
  
  strip_text <- if (region == "Schefferville") element_blank() else element_text()
  
  ggplot(data, aes(x = Ecotype, y = !!sym(y_var), color = Ecotype)) +
    geom_boxplot() +
    scale_color_manual(values = c("tomato", "cornflowerblue")) +
    stat_compare_means(method = "kruskal.test", label = "p.signif", hide.ns = TRUE, size = 3) +
    facet_wrap(.~Niche, scales = "free", ncol = 3) +
    stat_compare_means(label.y = NULL, size = 2) +
    stat_summary(fun = mean, colour = 1, geom = "point", shape = 5, size = 2, show.legend = FALSE) +
    scale_y_continuous(name = y_label) +
    theme(
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8, face = "plain"),
      legend.position = "none",
      axis.text.x = axis_text_x,
      axis.ticks.x = if (region == "Fermont") element_blank() else element_line(),
      axis.title.x = element_blank(),
      strip.text = strip_text,
      axis.text.y = element_text(face = "plain", size = 8, colour = "black"),
      axis.title.y = element_text(face = "plain", size = 8, colour = "black"),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      legend.key = element_blank()
    )
}

# --- ITS Analysis ---
process_ITS_data <- function() {
  # Read and preprocess data
  alpha_div_its <- read.delim("amplicon-data/ITS-alpha-diversity.txt", row.names = 1)
  env_its <- read.delim("amplicon-data/ITS-sample-metadata.txt", row.names = 1)
  df_its <- merge(alpha_div_its, env_its, by = "row.names") %>%
    filter(Year != "2021") %>%
    mutate(observed_features = as.numeric(observed_features))
  
  # Regions and variables for statistical tests
  regions <- c("Fermont", "Schefferville")
  variables <- c("observed_features", "shannon_entropy", "pielou_evenness")
  group_vars <- c("Niche", "Ecotype")
  
  # Plots and stats function
  generate_plots_and_stats <- function(df, create_plot_func, prefix = "") {
    plots <- list()
    for (region in regions) {
      data_region <- df %>% filter(Region == region)
      for (i in 1:length(variables)) {
        y_var <- variables[i]
        y_label <- switch(y_var,
                          observed_features = "Observed ASVs",
                          shannon_entropy = "Shannon-entropy",
                          pielou_evenness = "Pielou's evenness index")
        plot_name <- paste0("plot_", tolower(region), "_", i)
        plots[[plot_name]] <- create_plot_func(data_region, y_var, y_label, region)
      }
      
      # statistical analysis
      for (var in variables) {
        for (group in group_vars) {
          perform_stats(data_region, var, group, paste0(prefix, region))
        }
        for (niche in c("Root", "Bulk soil", "Rhizosphere")) {
          niche_data <- data_region %>% filter(Niche == niche)
          perform_stats(niche_data, var, "Ecotype", paste0(prefix, region, "-", tolower(niche)))
        }
      }
    }
    return(plots)
  }
  
  its_plots <- generate_plots_and_stats(df_its, create_plot)
  
  # Combine and display plots
  combined_plot <- (its_plots$plot_fermont_1 + its_plots$plot_fermont_2 + its_plots$plot_fermont_3) /
    (its_plots$plot_schefferville_1 + its_plots$plot_schefferville_2 + its_plots$plot_schefferville_3)
  print(combined_plot)
}

# --- 16S Analysis ---
process_16S_data <- function() {
  # Read and preprocess data
  alpha_div_16s <- read.delim("amplicon-data/16s-alpha-diversity.txt", row.names = 1)
  env_16s <- read.delim("amplicon-data/16s-sample-metadata.txt", row.names = 1)
  df_16s <- merge(alpha_div_16s, env_16s, by = "row.names") %>%
    filter(!Year %in% c("RH_2021", "BS_2021", "R_2021")) %>%
    mutate(observed_features = as.numeric(observed_features))
  
  # Regions and variables for statistical tests
  regions <- c("Fermont", "Schefferville")
  variables <- c("observed_features", "shannon_entropy", "pielou_evenness")
  group_vars <- c("Niche", "Ecotype")
  
  # Plots and stats function
  generate_plots_and_stats <- function(df, create_plot_func, prefix = "") {
    plots <- list()
    for (region in regions) {
      data_region <- df %>% filter(Region == region)
      for (i in 1:length(variables)) {
        y_var <- variables[i]
        y_label <- switch(y_var,
                          observed_features = "Observed ASVs",
                          shannon_entropy = "Shannon-entropy",
                          pielou_evenness = "Pielou's evenness index")
        plot_name <- paste0("plot_", tolower(region), "_", i)
        plots[[plot_name]] <- create_plot_func(data_region, y_var, y_label, region)
      }
      
      # statistical analysis
      for (var in variables) {
        for (group in group_vars) {
          perform_stats(data_region, var, group, paste0(prefix, region))
        }
        for (niche in c("Root", "Bulk soil", "Rhizosphere")) {
          niche_data <- data_region %>% filter(Niche == niche)
          perform_stats(niche_data, var, "Ecotype", paste0(prefix, region, "-", tolower(niche)))
        }
      }
    }
    return(plots)
  }
  
  # Using the streamlined function for 16S data
  its_plots <- generate_plots_and_stats(df_16s, function(data, y_var, y_label, region) create_plot(data, y_var, y_label, region, is_bacteria = TRUE))
  
  # Combine and display plots
  combined_plot <- (its_plots$plot_fermont_1 + its_plots$plot_fermont_2 + its_plots$plot_fermont_3) /
    (its_plots$plot_schefferville_1 + its_plots$plot_schefferville_2 + its_plots$plot_schefferville_3)
  print(combined_plot)
}

# Run the analysis for both datasets
process_16S_data()
process_ITS_data()
