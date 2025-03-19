# set directory
setwd("./GENOSCAN/")

# Load required packages efficiently
packages <- c("ampvis2", "tidyverse", "RColorBrewer", "colorRamps", 
              "geomtextpath", "ggforce", "patchwork")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

lapply(packages, install_if_missing)
lapply(packages, library, character.only = TRUE)

#-------------- ITS analysis ----------------
# Import fungal data
df_fungi <- amp_load(
  otutable = "amplicon-data/ITS-asv-table-all-runs-no-taxa.txt",
  taxonomy = "amplicon-data/ITS-taxonomy-all-runs.txt",
  metadata = "amplicon-data/ITS-sample-metadata.txt"
)

# View unique taxa (Phylum level)
unique(df_fungi$tax$Phylum)

# Filter by Region (Using dplyr for efficiency)
fermont_fungi <- df_fungi %>% amp_filter_samples(Region == "Fermont")
schefferville_fungi <- df_fungi %>% amp_filter_samples(Region == "Schefferville")

# Function to generate genus heatmaps
plot_genus_heatmap <- function(data, x_angle = 0, strip_text = TRUE) {
  theme_strip_text <- if (strip_text) element_text() else element_blank()
  
  amp_heatmap(
    data,
    group_by = "Ecotype",
    facet_by = "Niche",
    tax_aggregate = "Genus",
    tax_show = 30,
    color_vector = c("white", "firebrick4"),
    plot_colorscale = "log10",
    plot_values = FALSE
  ) +
    theme(
      axis.text.x = element_text(angle = x_angle, size = 8, vjust = 1),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      strip.text = theme_strip_text  # Fixing strip.text assignment
    )
}

# Generate heatmaps
genus_plot_fungi_F <- plot_genus_heatmap(fermont_fungi, x_angle = 0, strip_text = TRUE)
genus_plot_sch_F <- plot_genus_heatmap(schefferville_fungi, x_angle = 45, strip_text = FALSE)

# Combine and display plots
genus_plot_fungi_F + genus_plot_sch_F + plot_layout(ncol = 1)


#------------- 16S Analysis --------------------

# Import data
df_bacteria <- amp_load(otutable = "amplicon-data/16s-asv-table.biom")

# Rename taxa for plotting aesthetics
rename_taxa <- function(tax_table, old_names, new_names) {
  for (i in seq_along(old_names)) {
    tax_table$Genus[tax_table$Genus == old_names[i]] <- new_names[i]
  }
  return(tax_table)
}

df_bacteria$tax <- rename_taxa(
  df_bacteria$tax,
  c("Burkholderia-Caballeronia-Paraburkholderia", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"),
  c("Burkholderia", "Allorhizobium")
)

# Filter out unwanted taxa
df_bacteria <- amp_filter_taxa(df_bacteria, tax_vector = c("uncultured", "Subgroup_2"), remove = TRUE)

# Remove overlapping 2021 samples
df_bacteria <- amp_filter_samples(df_bacteria, !Year %in% c("BS_2021", "R_2021", "RH_2021"))

# Subset data by region
fermont_bact <- amp_filter_samples(df_bacteria, Region == "Fermont")
schefferville_bact <- amp_filter_samples(df_bacteria, Region == "Schefferville")

# Function to generate genus-level heatmap
plot_genus_heatmap <- function(data, x_angle = 0) {
  amp_heatmap(data,
              group_by = "Ecotype",
              facet_by = "Niche",
              tax_aggregate = "Genus",
              tax_show = 30,
              tax_empty = "none",
              normalise = TRUE,
              color_vector = c("white", "firebrick4"),
              plot_colorscale = "log10",
              plot_values = FALSE) +
    theme(axis.text.x = element_text(angle = x_angle, size = 8, vjust = 1),
          axis.text.y = element_text(size = 8),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "right")
}

# Generate heatmaps
genus_plot_bact_F <- plot_genus_heatmap(fermont_bact, x_angle = 0) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

genus_plot_sch_B <- plot_genus_heatmap(schefferville_bact, x_angle = 45) +
  theme(strip.text = element_blank())

# Combine plots
combined_plot <- genus_plot_bact_F + genus_plot_sch_B + plot_layout(ncol = 1)
print(combined_plot)
