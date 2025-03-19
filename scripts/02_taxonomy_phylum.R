# Load packages
library(ampvis2)
library(tidyverse)
library(RColorBrewer)
library(colorRamps)
library(phyloseq)
library(patchwork)
library(scales)

# Set directory
setwd("./GENOSCAN/")

# ---------------------- AMPLICON DATA -------------------
# Function to process and plot phylum data
process_and_plot_phylum <- function(phyloseq_obj, data_type) {
  # Subset and process data
  phyloseq_obj <- subset_samples(phyloseq_obj, !Year=="2021")
  phylum_data <- tax_glom(phyloseq_obj, taxrank = 'Phylum') %>%
    psmelt() %>%
    group_by(Region, Ecotype, Niche) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    group_by(Region, Ecotype, Niche, Site, Phylum) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop')
  
  # Rename rare phyla and reorder
  phylum_data$Phylum[phylum_data$Abundance < 0.01] <- "Others"
  phylum_data <- phylum_data %>%
    mutate(Phylum = fct_reorder(Phylum, Abundance, .desc = T))
  
  # Create color palette
  all_phylum <- unique(phylum_data$Phylum)
  colourCount <- length(all_phylum)
  color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)
  
  # Plot
  plot <- phylum_data %>% 
    count(Niche, Region, Phylum, Site, Ecotype, wt=Abundance, name = "Abundance") %>%
    ggplot() + 
    geom_bar(aes(x=Ecotype, y=Abundance, fill=Phylum), stat="identity", position="stack") +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(name = "Mean relative abundance", labels = scales::percent, breaks = seq(0, 1, by = 0.1)) +
    scale_x_discrete(name = "") +
    theme_minimal(base_line_size = 0.05, base_rect_size = 0.05) +
    facet_grid(Region~Niche, scales = "free", space = "free") +
    theme(
      legend.title = element_text(size = 8, face = "bold", angle = 90),
      legend.text = element_text(size = 8, face = "plain"),
      legend.position = "bottom",
      axis.text.x = element_text(colour = "black", size = 8, face = "plain", angle = 60, vjust = 1, hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_text(colour = "black", size = 8, face = "plain", angle = 0, vjust = 0.3, hjust = 1),
      axis.title.y = element_text(face = "bold", size = 8, colour = "black"),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      legend.key = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 3))
  
  # Adjust y-axis for fungi plot
  if (data_type == "fungi") {
    plot <- plot +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }
  
  return(plot)
}

# Process and plot 16S rRNA data
biom_file <- "amplicon-data/16s-asv-table.biom"
cc <- import_biom(biom_file)
colnames(tax_table(cc)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
bacteria_phylum_plot <- process_and_plot_phylum(cc, "bacteria")

# Process and plot ITS data
df_fungi <- amp_load(
  otutable = "amplicon-data//ITS-asv-table-all-runs-no-taxa.txt",
  taxonomy = "amplicon-data//ITS-taxonomy-all-runs.txt",
  metadata = "amplicon-data//ITS-sample-metadata.txt"
)
cc_fungi <- phyloseq(
  otu_table(df_fungi$abund, taxa_are_rows = TRUE),
  tax_table(as.matrix(df_fungi$tax)),
  sample_data(df_fungi$metadata)
)
fungi_phylum_plot <- process_and_plot_phylum(cc_fungi, "fungi")

# Combine plots
combined_plot <- bacteria_phylum_plot + fungi_phylum_plot +
  plot_layout(ncol = 2, nrow = 1)

# Display combined plot
print(combined_plot)

# ------------------- SHOTGUN analysis -------------------
# Read data
abund_table <- read.delim("shotgun-data/kaiju_summary_kingdom.tsv")
smd <- read.delim("shotgun-data/sample-metadata.txt")

# Merge datasets
df_main <- merge(abund_table, smd, by = "sampleid")

# Pre-filtering: Keep only Bacteria, Archaea, and Fungi (ending in "mycota")
df_main <- df_main %>% filter(Kingdom %in% c("Bacteria", "Archaea") | grepl("mycota$", Phylum))

# Function to process each kingdom
process_phylum_data <- function(df, kingdom_filter, abundance_threshold, is_fungi = FALSE) {
  if (is_fungi) {
    df <- df[grep("mycota$", df$Phylum), ]  # Proper fungi filtering
  } else {
    df <- df %>% filter(Kingdom == kingdom_filter)
  }
  
  if (nrow(df) == 0) {
    message(paste("Warning: No data found for", kingdom_filter))
    return(list(df = data.frame(Region = character(), Ecotype = character(), Phylum = character(), Abundance = numeric()), palette = c()))
  }
  
  # Calculate relative abundance
  df <- df %>%
    group_by(Region, Ecotype) %>%
    mutate(Abundance = reads / sum(reads)) %>%
    ungroup()
  
  # Sum ASVs according to phylum
  df <- df %>%
    group_by(Region, Ecotype, Phylum) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop')
  
  # Rename rare phyla with <1% abundance
  df$Phylum[df$Abundance < abundance_threshold] <- "Others"
  
  # Reorder Phylum factor based on Abundance
  df <- df %>%
    mutate(Phylum = fct_reorder(Phylum, Abundance, .desc = TRUE))
  
  # Get unique phyla
  all_phylum <- unique(df$Phylum)
  colourCount <- length(all_phylum)
  
  # Create color palette
  color_palette <- setNames(
    colorRampPalette(brewer.pal(12, "Set3"))(max(colourCount, 1)), 
    all_phylum
  )
  
  return(list(df = df, palette = color_palette))
}


# Process each group
df_bact_result <- process_phylum_data(df_main, "Bacteria", 0.01)
df_bact <- df_bact_result$df
palette_bact <- df_bact_result$palette

df_fungi_result <- process_phylum_data(df_main, "Eukaryota", 0.01, is_fungi = TRUE)  # ✅ Correct fungi subsetting
df_fungi <- df_fungi_result$df
palette_fungi <- df_fungi_result$palette

df_archaea_result <- process_phylum_data(df_main, "Archaea", 0.01)
df_archaea <- df_archaea_result$df
palette_archaea <- df_archaea_result$palette

# Function to generate plots
generate_phylum_plot <- function(df, palette, legend_cols = 3, legend_pos = "bottom") {
  ggplot(df, aes(x = Ecotype, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = palette) +
    scale_y_continuous(name = "Mean relative abundance", labels = percent, breaks = seq(0, 1, by = 0.1)) +
    scale_x_discrete(name = "") +
    facet_grid(. ~ Region, scales = "free", space = "free") +
    theme_minimal(base_line_size = 0.05, base_rect_size = 0.05) +
    theme(
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8),
      legend.position = legend_pos,
      axis.text.x = element_text(colour = "black", size = 8, angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 8),
      axis.title.y = element_text(face = "bold", size = 8, colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black"),
      legend.key = element_blank()
    ) +
    guides(fill = guide_legend(ncol = legend_cols))
}

# Generate and display plots
df_bact_plot <- generate_phylum_plot(df_bact, palette_bact)
df_fungi_plot <- generate_phylum_plot(df_fungi, palette_fungi)
df_archaea_plot <- generate_phylum_plot(df_archaea, palette_archaea, legend_cols = 1, legend_pos = "right")

df_bact_plot
df_fungi_plot
df_archaea_plot
