# Load packages
library(dplyr)
library(stringr)
library(pacman)
library(IRanges)
p_load(phyloseq, tidyverse, magrittr, ComplexHeatmap, colorRamp2, circlize,
       purrr, patchwork, kableExtra)

# Set working directory
setwd("./GENOSCAN/")

# Load module completeness data
df <- read.delim("metabolism_MAGS101_no_coverage_modules.txt") 

# Define non-prokaryotic pathways to remove
non_prokaryotic_pathways <- c(
  "CAM \\(Crassulacean acid metabolism\\), dark",
  "CAM \\(Crassulacean acid metabolism\\), light",
  "C4−dicarboxylic acid cycle, NADP − malic enzyme type",
  "C4−dicarboxylic acid cycle, NAD − malic enzyme type",
  "C4−dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type",
  "Reductive pentose phosphate cycle, ribulose−5P",
  "Reductive pentose phosphate cycle \\(Calvin cycle\\)",
  "Reductive pentose phosphate cycle, glyceraldehyde−3P"
)

# Clean and filter module names
cleaned_df <- df %>%
  mutate(
    module_name = str_trim(module_name),
    module_name = str_replace_all(module_name, "[-−‐‑–—]", "-")
  ) %>%
  filter(!str_detect(module_name, paste(non_prokaryotic_pathways, collapse = "|")))

# Reshape data for heatmap
reshaped_data <- cleaned_df %>%
  select(module, bin_name, pathwise_module_completeness) %>%
  pivot_wider(names_from = bin_name, values_from = pathwise_module_completeness)

# Extract module metadata
metrics <- c('module_name', 'module_subcategory')
module_metrics <- cleaned_df %>%
  select(module, all_of(metrics)) %>%
  distinct()

# Merge module completeness with metadata
merged_data <- module_metrics %>%
  left_join(reshaped_data, by = "module")

# Load species taxonomy data
species_taxonomy <- read.delim("taxonomy-main-2024-11-15.txt")

# Define target metabolic pathways of interest
pathway_groups_of_interest <- c("Symbiosis",
                                "Photosynthesis",
                                "Nitrogen metabolism",
                                "Methane metabolism", 
                                "Carbon fixation"
)
modules_of_interest <- c("M00173", "M00376", "M00375", "M00374",
                         "M00377", "M00579", "M00260","M00664",
                         "M00672")

# Filter pathways of interest
filtered_pathways <- merged_data %>% 
  filter(module_subcategory %in% pathway_groups_of_interest | module %in% modules_of_interest) %>%
  mutate(across(module_subcategory, as_factor)) %>% 
  mutate(module_subcategory = fct_relevel(module_subcategory, c("Symbiosis",
                                                                "Photosynthesis",
                                                                "Nitrogen metabolism",
                                                                "Methane metabolism", 
                                                                "Carbon fixation"))) %>% 
  arrange(module_subcategory)

write.csv(filtered_pathways, "SUBSET-of-KEGG-MODULES-of INTEREST.csv")
# Convert NA to 0.00
filtered_pathways <- filtered_pathways %>% mutate(across(everything(), ~ replace_na(., 0)))

# Transpose data to match heatmap format
transposed_pathways <- filtered_pathways %>% 
  dplyr::select(-module, -module_subcategory, -module_name) %>% 
  t() %>% 
  data.frame() %>% 
  setNames(filtered_pathways$module) %>% 
  rownames_to_column(var = "id") 
# Merge with species data
merged_pathways_species <- left_join(transposed_pathways, species_taxonomy, by = "id")

# Sort by taxonomy
sorted_data <- merged_pathways_species %>% 
  arrange(Phylum) %>% 
  arrange(desc(Ecotype))

# Format data for heatmap
formatted_pathways <- sorted_data %>% 
  tibble::column_to_rownames(var = "Species") %>% 
  select(starts_with("M0")) 
pathway_groups <- filtered_pathways %>% select(module, module_subcategory)
pathway_descriptions <- filtered_pathways %>% select(module, module_name)

# Select only columns that total more than zero 
filtered_data <- formatted_pathways %>% 
  select_if(colSums(.) > 0) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "module")

# Make grouping variable
grouped_pathways <- left_join(filtered_data, pathway_groups) %>% select(module, module_subcategory)
pathway_group_vector <- grouped_pathways %>% column_to_rownames(var = "module") 
group_names <- rownames(pathway_group_vector) 
unlisted_groups <- unlist(pathway_group_vector)
grouping_structure <- structure(unlisted_groups, names = group_names)

# Make module names more meaningful
pathway_descriptions_cleaned <- left_join(filtered_data, pathway_descriptions) %>% select(module, module_name)
pathway_descriptions_cleaned$module_name <- sub(" =>.*", "", pathway_descriptions_cleaned$module_name) 
description_vector <- pathway_descriptions_cleaned %>% column_to_rownames(var = "module") 
description_names <- rownames(description_vector) 
unlisted_descriptions <- unlist(description_vector)
description_labels <- structure(unlisted_descriptions, names = description_names)
heatmap_data <- filtered_data %>% tibble::column_to_rownames(var="module")

# Order category 
order_vector <- sorted_data %>% 
  select(Species, Phylum) %>%
  column_to_rownames(var = "Species")
order_names <- rownames(order_vector) 
unlisted_order <- unlist(order_vector)
order_grouping <- structure(unlisted_groups, names = group_names)

# Define color scale for heatmap
color_function <- colorRamp2(c(0, 0.2, 0.6, 0.8, 1), c("white", "#BA8A65", "#B37A5F", "#AB6A58", "#A45A51"))

# Generate a color palette for Phylum
unique_phyla <- unique(sorted_data$Phylum)
phylum_colors <- setNames(
  colorRampPalette(c("#2986ccff","#6a329fff","#8fce00ff","#c90076ff","#c9a39bff",
                     "#f44336ff","#744700ff","#ce7e00ff"))(length(unique_phyla)), 
  unique_phyla
)

# Create the Heatmap
Heatmap(as.matrix(heatmap_data), 
        name = "Module Completeness",  # Name for the main heatmap legend
        col = color_function,
        column_split = sorted_data$Ecotype,
        row_split = grouping_structure, 
        column_names_gp = gpar(fontsize = 24),
        column_title_gp = gpar(fontsize = 34, fontface = "bold"),
        column_names_side = "bottom",
        column_names_rot = 36,
        column_gap = unit(5, "mm"),
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        row_names_max_width = max_text_width(description_labels),
        show_row_dend = FALSE,
        row_labels = description_labels,
        row_gap = unit(5, "mm"),
        row_names_gp = gpar(fontsize = 26),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 32),
        cluster_rows = FALSE, 
        cluster_row_slices = FALSE, 
        bottom_annotation = HeatmapAnnotation(
          Phylum = sorted_data$Phylum,
          col = list(Phylum = phylum_colors),
          annotation_legend_param = list(Phylum = list(nrow = 3, title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20))),
          height = unit(8, "cm"),
          show_annotation_name = TRUE,
          annotation_name_gp = gpar(fontsize = 24),
          border = TRUE,
          simple_anno_size = unit(2, "cm")
        ),
        top_annotation = HeatmapAnnotation(
          Compartment = sorted_data$Ecotype,
          col = list(Compartment = c("NT" = "cornflowerblue", "DT" = "coral")),
          show_legend = TRUE,
          annotation_legend_param = list(Compartment = list(title_gp = gpar(fontsize = 24), labels_gp = gpar(fontsize = 20))),
          show_annotation_name = TRUE,
          annotation_name_gp = gpar(fontsize = 24),
          simple_anno_size = unit(2, "cm")
        ),
        width = unit(90, "cm"),
        height = unit(55, "cm"),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 24, fontface = "bold"),
          labels_gp = gpar(fontsize = 20),
          legend_height = unit(20, "cm"),
          grid_width = unit(1, "cm")
        ),
        border = TRUE
)

# Legends
legend_list <- list(Order = c(
  "Acidobacteriota" = "#2986ccff", "Actinobacteriota" = "#6a329fff",
  "Actinomycetota" = "#c9a39bff", "Bacteroidota" = "#f44336ff",
  "Myxococcota" = "#744700ff", "Planctomycetota" = "#8fce00ff",
  "Proteobacteria" = "#c90076ff", "Verrucomicrobiota" = "#ce7e00ff" 
))

legend_colors <- unlist(legend_list)
taxa_names <- names(legend_list[[1]]) 

phylum_legend <- Legend(labels = taxa_names, title = "Phylum", 
                        legend_gp = gpar(fill = legend_colors))

completeness_legend <- Legend(col_fun = color_function, title = "Pathway \nCompleteness (%)", border = "black")

draw(completeness_legend, just = c("right"))
draw(phylum_legend, just = c("left"))
