# Load necessary libraries
library(tidyverse)
library(ggpubr)
library(phyloseq)
library(ampvis2)
library(vegan)
library(here)

# Set working directory
data_dir <- here("./GENOSCAN/")

# Import data
data <- read_tsv(file.path(data_dir, "KEGG-KO-GENE-abundance-GPM-short.txt"))
taxonomy <- read.delim(file.path(data_dir, "annotation-KOfam.txt"))
smd <- read.delim(file.path(data_dir, "sample-metadata.txt"))

# Create an abundance matrix and remove duplicates, if any
df_abundance <- data %>%
  select(1, 3:22) %>%
  distinct(across(1), .keep_all = TRUE) %>%
  column_to_rownames(var = names(.)[1])

# Convert taxonomy to matrix format
taxonomy <- as.matrix(taxonomy)

# Load data into ampvis2 object
dfc <- amp_load(otutable = df_abundance, taxonomy = taxonomy, metadata = smd)

# Normalize numeric metadata columns
numeric_cols <- names(dfc$metadata)[sapply(dfc$metadata, is.numeric)]
dfc$metadata[numeric_cols] <- decostand(dfc$metadata[numeric_cols], method = "standardize")

# Display mean and standard deviation
round(colMeans(dfc$metadata[numeric_cols], na.rm = TRUE), 1)
apply(dfc$metadata[numeric_cols], 2, sd, na.rm = TRUE)

# Function to customize plots
customize_plot <- function(plot) {
  plot +
    stat_ellipse(aes(label = Ecotype, group = Ecotype), color = "gray50", 
                 geom = "textpath", hjust = 0.7, vjust = 1.2, linetype = 2) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    scale_color_manual(values = c("coral", "cornflowerblue"), name = "Ecotype") +
    scale_shape_manual("Site", values = c(15, 16)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA))
}

# Perform RDA (Redundancy Analysis)
rda_result <- amp_ordinate(
  dfc, type = "rda", transform = "hellinger", distmeasure = "bray", 
  constrain = "Ecotype", envfit_numeric = c(11, 14:31, 34, 35), 
  envfit_signif_level = 0.05, detailed_output = TRUE
)

# Display RDA results
rda_result$model
RsquareAdj(rda_result$model)
anova.cca(rda_result$model)
rda_result$evf_numeric_model

# Plot RDA
amp_ordinate(
  dfc, type = "rda", transform = "hellinger", distmeasure = "bray", 
  sample_color_by = "Ecotype", sample_shape_by = "Region",
  constrain = "Ecotype", envfit_numeric = c(11, 14:31, 34, 35), 
  envfit_signif_level = 0.05
) %>%
  customize_plot()

# Plot PCoA
amp_ordinate(
  dfc, type = "pcoa", transform = "none", distmeasure = "bray", 
  sample_color_by = "Ecotype", sample_shape_by = "Region"
) %>%
  customize_plot()

# Perform PERMANOVA (adonis2)
permanova_results <- list(
  Ecotype = adonis2(decostand(t(dfc$abund), method = "hellinger") ~ Ecotype, data = dfc$metadata),
  Region = adonis2(decostand(t(dfc$abund), method = "hellinger") ~ Region, data = dfc$metadata),
  Ecotype_Region = adonis2(decostand(t(dfc$abund), method = "hellinger") ~ Ecotype + Region, data = dfc$metadata)
)

# Print PERMANOVA results
permanova_results
