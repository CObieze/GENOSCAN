# set directory
setwd("./GENOSCAN/")

# Load necessary libraries
library(ampvis2)
library(tidyverse)
library(geomtextpath)
library(RColorBrewer)
library(colorRamps)

############################# Bacteria Analysis #############################

# Import bacterial data and metadata
df <- amp_load(otutable = "amplicon-data/16s-asv-table.biom")

# Remove overlapping 2021 samples
df <- amp_filter_samples(df, !Year %in% c("RH_2021", "BS_2021", "R_2021"))

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

### Function for PCoA Plot ###
create_pcoa_plot <- function(data, hjust_val) {
  data <- filter_otus(data, filter_otus = 1.0)
  amp_ordinate(data, 
               type = "pcoa",
               constrain = c("Ecotype", "Site"),
               transform = "none",
               distmeasure = "bray",
               sample_color_by = "Ecotype",
               sample_shape_by = "Site",
               species_plot = FALSE,
               species_point_size = 0,
               sample_colorframe = FALSE,
               sample_label_size = 6,
               sample_point_size = 3,
               detailed_output = FALSE) +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = hjust_val, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15:19, 7, 8, 9, 10, 13)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA))
}

# Generate PCoA plots
BS_pcoa_fer <- create_pcoa_plot(fermont_soil, 0.40)
RH_pcoa_fer <- create_pcoa_plot(fermont_rhizosphere, 0.40)
Root_pcoa_fer <- create_pcoa_plot(fermont_roots, 0.50)

BS_pcoa_sch <- create_pcoa_plot(schefferville_soil, 0.30)
RH_pcoa_sch <- create_pcoa_plot(schefferville_rhizosphere, 0.55)
Root_pcoa_sch <- create_pcoa_plot(schefferville_roots, 0.20)

# Display plots
BS_pcoa_fer
RH_pcoa_fer
Root_pcoa_fer
BS_pcoa_sch
RH_pcoa_sch
Root_pcoa_sch


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

# Generate PCoA plots
Fungi_BS_pcoa_fer <- create_pcoa_plot(fermont_fungi_soil, 0.40)
Fungi_RH_pcoa_fer <- create_pcoa_plot(fermont_fungi_rhizosphere, 0.40)
Fungi_Root_pcoa_fer <- create_pcoa_plot(fermont_fungi_roots, 0.40)

Fungi_BS_pcoa_sch <- create_pcoa_plot(schefferville_fungi_soil, 0.30)
Fungi_RH_pcoa_sch <- create_pcoa_plot(schefferville_fungi_rhizosphere, 0.55)
Fungi_Root_pcoa_sch <- create_pcoa_plot(schefferville_fungi_roots, 0.20)

# Display plots
Fungi_BS_pcoa_fer
Fungi_RH_pcoa_fer
Fungi_Root_pcoa_fer
Fungi_BS_pcoa_sch
Fungi_RH_pcoa_sch
Fungi_Root_pcoa_sch

