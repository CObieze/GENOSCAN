# Load necessary libraries
library(ampvis2)
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(mia)
library(patchwork)

#-------------------- Bacteria --------------------

# Set working directory
setwd("C:/Users/Admin/Downloads/GENOSCAN-LAB-FILES/shotgun-sequences/Sequence_processing_output/manuscript/")

# Import ASV table
df <- amp_load(otutable = "data-2025-03-04/amplicon-data/16S-asv-table.biom")
# Make the replacement to shorten name for plotting aestetics
df$tax$Genus[df$tax$Genus=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium"

####### subset the table to samples from Fermont
fermont <- amp_filter_samples(
  df, Region %in% c("Fermont")
)

# subset according to soil compartment (Root) and filter low abundance samples for consistency
root_fermont <- amp_filter_samples(
  fermont, minreads = 2000, Niche %in% c("Root"))
# remove overlapping 2021 samples
root_fermont <- amp_filter_samples(root_fermont, !Year %in% c("R_2021"))

# subset according to soil compartment (bulk soil) and filter low abundance samples for consistency
soil_fermont <- amp_filter_samples(
  fermont, minreads = 6000, Niche %in% c("Bulk soil"))
# remove overlapping 2021 samples
soil_fermont <- amp_filter_samples(soil_fermont, !Year %in% c("BS_2021"))

# subset according to soil compartment (Rhizosphere) and filter low abundance samples
rhizosphere_fermont <- amp_filter_samples(
  fermont, minreads = 6000, Niche %in% c("Rhizosphere"))
# remove overlapping 2021 samples
rhizosphere_fermont <- amp_filter_samples(rhizosphere_fermont, !Year %in% c("RH_2021"))

####### subset the table to samples from Schefferville (no 2021 overlapping samples)
schefferville <- amp_filter_samples(
  df, Region %in% c("Schefferville")
)

# subset according to soil compartment (Root) and filter low abundance samples
root_sch <- amp_filter_samples(
  schefferville, minreads = 2000, Niche %in% c("Root"))

# subset according to soil compartment (bulk soil) and filter low abundance samples
soil_sch <- amp_filter_samples(
  schefferville, minreads = 6000, Niche %in% c("Bulk soil"))

# subset according to soil compartment (Rhizosphere) and filter low abundance samples
rhizosphere_sch <- amp_filter_samples(
  schefferville, minreads = 6000, Niche %in% c("Rhizosphere"))

# Define ANCOMBC2 parameters
ancombc2_params <- list(
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Ecotype",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo = TRUE,
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE
)

# Function to create phyloseq objects
create_phyloseq <- function(ampvis_object, metadata_columns) {
  phyloseq(
    otu_table(ampvis_object$abund, taxa_are_rows = TRUE),
    tax_table(as.matrix(ampvis_object$tax)),
    sample_data(ampvis_object$metadata[metadata_columns])
  )
}

# Function to run ANCOMBC2 and save results
run_ancombc2 <- function(tse_data, output_filename) {
  # Ensure the results directory exists
  results_dir <- dirname(output_filename)
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  set.seed(123)
  
  tryCatch({
    output <- do.call(ancombc2, c(list(data = tse_data), ancombc2_params))
    write.csv(output$res, output_filename)
    return(output$res)
  }, error = function(e) {
    zero_var_taxa <- str_extract_all(conditionMessage(e), "[a-f0-9]{32}")[[1]]
    
    if (length(zero_var_taxa) > 0) {
      message("Removing zero-variance taxa and retrying ANCOM-BC2...")
      tse_data <- tse_data[!rownames(tse_data) %in% zero_var_taxa, ]
      output <- do.call(ancombc2, c(list(data = tse_data), ancombc2_params))
      write.csv(output$res, output_filename)
      return(output$res)
    } else {
      stop("Unexpected ANCOM-BC2 error: ", conditionMessage(e))
    }
  })
}

# Process each region and compartment
compartments_order <- c("soil_Fermont", "soil_Schefferville", "rhizosphere_Fermont", "rhizosphere_Schefferville", "root_Fermont", "root_Schefferville")
index <- 1

# make nested list of sites 
regions <- list(
  "Fermont" = list(
    soil = soil_fermont, 
    root = root_fermont, 
    rhizosphere = rhizosphere_fermont
  ),
  "Schefferville" = list(
    soil = soil_sch, 
    root = root_sch, 
    rhizosphere = rhizosphere_sch
  )
)

# plot function
plot_taxa <- function(df, title_text = NULL) {
  ggplot(df, aes(x = taxon, y = lfc_EcotypeDisturbed, fill = direct)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_EcotypeDisturbed - se_EcotypeDisturbed, 
                      ymax = lfc_EcotypeDisturbed + se_EcotypeDisturbed), 
                  width = 0.2, position = position_dodge(0.05), color = "black") +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(name = NULL, values = c("coral", "cornflowerblue")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(color = df$color),
          legend.position = "none") +
    coord_flip()
}

# Process each region and compartment
for (region in names(regions)) {
  for (compartment in names(regions[[region]])) {
    phyloseq_obj <- create_phyloseq(regions[[region]][[compartment]], c(1, 12, 30, 34, 36))
    tse_obj <- mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_obj)
    tse_obj$Ecotype <- factor(tse_obj$Ecotype, levels = c("Natural", "Disturbed"))
    
    output_filename <- paste0("results/ANCOMBC2_DAA_", compartment, "_", region, ".csv")
    ancom_results <- run_ancombc2(tse_obj, output_filename)
    
    # Subset to top 15 differentially abundant taxa and split into positive and negative LFC
    df_fig_disturbed_negative <- ancom_results %>%
      dplyr::select(taxon, ends_with("Disturbed")) %>%
      dplyr::mutate(taxon = stringr::str_replace(taxon, "^Genus:", "")) %>%
      dplyr::filter(!stringr::str_detect(taxon, "uncultured")) %>%
      dplyr::filter(diff_EcotypeDisturbed == 1, lfc_EcotypeDisturbed < 0, q_EcotypeDisturbed < 0.05) %>%
      dplyr::arrange(lfc_EcotypeDisturbed) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(direct = "Natural", color = "black") %>%
      dplyr::arrange(desc(lfc_EcotypeDisturbed))
    
    df_fig_disturbed_positive <- ancom_results %>%
      dplyr::select(taxon, ends_with("Disturbed")) %>%
      dplyr::mutate(taxon = stringr::str_replace(taxon, "^Genus:", "")) %>%
      dplyr::filter(!stringr::str_detect(taxon, "uncultured")) %>%
      dplyr::filter(diff_EcotypeDisturbed == 1, lfc_EcotypeDisturbed > 0, q_EcotypeDisturbed < 0.05) %>%
      dplyr::arrange(desc(lfc_EcotypeDisturbed)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(direct = "Disturbed", color = "black")
    
    df_fig_disturbed <- bind_rows(df_fig_disturbed_positive, df_fig_disturbed_negative)
    
    # Reorder factors for plotting
    df_fig_disturbed$taxon <- factor(df_fig_disturbed$taxon, levels = df_fig_disturbed$taxon)
    df_fig_disturbed$direct <- factor(df_fig_disturbed$direct, levels = c("Disturbed", "Natural"))
    
    assign(paste0("plot_", compartment, "_", region), plot_taxa(df_fig_disturbed, NULL), envir = .GlobalEnv)
  }
}

plot_soil_Fermont<-plot_soil_Fermont +
  labs(title = "Fermont: Bacteria") +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
plot_soil_Schefferville<-plot_soil_Schefferville+
  labs(title = "Schefferville: Bacteria") 

# merge plots
plot_soil_Fermont+plot_soil_Schefferville+plot_rhizosphere_Fermont+plot_rhizosphere_Schefferville+plot_root_Fermont+plot_root_Schefferville+
  plot_layout(ncol = 2, nrow = 3)


######################################################################################################
############################## Fungi ######################################
######################################################################################################
# Set working directory for Fungi
setwd("C:/Users/Admin/Downloads/GENOSCAN-LAB-FILES/Sequence-runs/concanated-seqs-run/run1_3/ITS/Decontamination/no-decontamination-tables/")

# Import fungal biom and metadata
df_fungi <- amp_load(otutable = "data/ITS-asv-table-all-runs-no-taxa.txt",
                     taxonomy = "data/ITS-taxonomy-all-runs.txt",
                     metadata = "data/ITS-sample-metadata.txt")

# Subset the table to samples from Fermont
fermont_fungi <- amp_filter_samples(
  df_fungi, Region %in% c("Fermont")
)

# subset according to soil compartment (Root) and filter low abundance samples
root_fermont_fungi <- amp_filter_samples(
  fermont_fungi, minreads = 2000, Niche %in% c("Root"))
# remove overlapping 2021 samples
root_fermont_fungi <- amp_filter_samples(root_fermont_fungi, !Year %in% c("2021"))

# subset according to soil compartment (bulk soil) and filter low abundance samples
soil_fermont_fungi <- amp_filter_samples(
  fermont_fungi, minreads = 6000, Niche %in% c("Bulk soil"))
# remove overlapping 2021 samples
soil_fermont_fungi <- amp_filter_samples(soil_fermont_fungi, !Year %in% c("2021"))

# subset according to soil compartment (Rhizosphere) and filter low abundance samples
rhizosphere_fermont_fungi <- amp_filter_samples(
  fermont_fungi, minreads = 6000, Niche %in% c("Rhizosphere"))
# remove overlapping 2021 samples
rhizosphere_fermont_fungi <- amp_filter_samples(rhizosphere_fermont_fungi, !Year %in% c("2021"))

# subset table to schefferville (no overlapping 2021 samples)
schefferville_fungi <- amp_filter_samples(
  df_fungi, Region %in% c("Schefferville")
)

# subset according to soil compartment (root)
root_sch_fungi <- amp_filter_samples(
  schefferville_fungi, minreads = 2000, Niche %in% c("Root"))

# subset according to compartment (Bulk soil)
soil_sch_fungi <- amp_filter_samples(
  schefferville_fungi, minreads = 6000, Niche %in% c("Bulk soil"))

# subset according to compartment (rhizosphere)
rhizosphere_sch_fungi <- amp_filter_samples(
  schefferville_fungi, minreads = 6000, Niche %in% c("Rhizosphere"))

# Process each region and compartment
compartments_order <- c("soil_fermont_fungi", "soil_sch_fungi", "rhizosphere_fermont_fungi", "rhizosphere_sch_fungi", "root_fermont_fungi", "root_sch_fungi")
index <- 1
regions <- list(
  "Fermont" = list(
    soil = soil_fermont_fungi, 
    root = root_fermont_fungi, 
    rhizosphere = rhizosphere_fermont_fungi
  ),
  "Schefferville" = list(
    soil = soil_sch_fungi, 
    root = root_sch_fungi, 
    rhizosphere = rhizosphere_sch_fungi
  )
)

# Process each region and compartment for fungi
for (region in names(regions)) {
  for (compartment in names(regions[[region]])) {
    phyloseq_obj <- create_phyloseq(regions[[region]][[compartment]], c(1:8))
    tse_obj <- mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq_obj)
    tse_obj$Ecotype <- factor(tse_obj$Ecotype, levels = c("Natural", "Disturbed"))
    
    output_filename <- paste0("results_fungi/ANCOMBC2_DAA_", compartment, "_", region, ".csv")
    ancom_results <- run_ancombc2(tse_obj, output_filename)
    
    # Check if ANCOM-BC2 returned results
    if (is.null(ancom_results) || nrow(ancom_results) == 0) {
      message("No differentially abundant taxa found for ", compartment, " in ", region, ". Skipping...")
      next
    }
    
    # Print column names to debug
    print(colnames(ancom_results))
    
    # Adjust selection if needed
    if (!any(grepl("Disturbed", colnames(ancom_results)))) {
      stop("Expected column pattern 'Disturbed' not found. Check column names!")
    }
    
    # Subset to top 15 differentially abundant taxa and split into positive and negative LFC
    df_fig_disturbed_negative <- ancom_results %>%
      dplyr::select(taxon, ends_with("Disturbed")) %>%
      dplyr::mutate(taxon = stringr::str_replace(taxon, "^Genus:", "")) %>%
      dplyr::filter(!stringr::str_detect(taxon, "uncultured")) %>%
      dplyr::filter(diff_EcotypeDisturbed == 1, lfc_EcotypeDisturbed < 0, q_EcotypeDisturbed < 0.05) %>%
      dplyr::arrange(lfc_EcotypeDisturbed) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(direct = "Natural", color = "black") %>%
      dplyr::arrange(desc(lfc_EcotypeDisturbed))
    
    df_fig_disturbed_positive <- ancom_results %>%
      dplyr::select(taxon, ends_with("Disturbed")) %>%
      dplyr::mutate(taxon = stringr::str_replace(taxon, "^Genus:", "")) %>%
      dplyr::filter(!stringr::str_detect(taxon, "uncultured")) %>%
      dplyr::filter(diff_EcotypeDisturbed == 1, lfc_EcotypeDisturbed > 0, q_EcotypeDisturbed < 0.05) %>%
      dplyr::arrange(desc(lfc_EcotypeDisturbed)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(direct = "Disturbed", color = "black")
    
    df_fig_disturbed <- bind_rows(df_fig_disturbed_positive, df_fig_disturbed_negative)
    
    # Reorder factors for plotting
    df_fig_disturbed$taxon <- factor(df_fig_disturbed$taxon, levels = df_fig_disturbed$taxon)
    df_fig_disturbed$direct <- factor(df_fig_disturbed$direct, levels = c("Disturbed", "Natural"))
    
    assign(paste0("plot_fungi_", compartment, "_", region), plot_taxa(df_fig_disturbed, NULL), envir = .GlobalEnv)
  }
}

# slightly customize plots further
plot_fungi_soil_Fermont<-plot_fungi_soil_Fermont+
  labs(title = "Fermont: Fungi") 
plot_fungi_soil_Schefferville<-plot_fungi_soil_Schefferville+
  labs(title = "Schefferville: Fungi") 
plot_fungi_root_Fermont<-plot_fungi_root_Fermont+
  labs(x = NULL, y = "Log fold change")
plot_fungi_root_Schefferville<-plot_fungi_root_Schefferville+
  labs(x = NULL, y = "Log fold change")

# Merge plots
plot_fungi_soil_Fermont+plot_fungi_soil_Schefferville+plot_fungi_rhizosphere_Fermont+plot_fungi_rhizosphere_Schefferville+plot_fungi_root_Fermont+plot_fungi_root_Schefferville+
  plot_layout(ncol = 2, nrow = 3)

# Merge bacteria and fungi
plot_soil_Fermont+plot_soil_Schefferville+plot_rhizosphere_Fermont+plot_rhizosphere_Schefferville+plot_root_Fermont+plot_root_Schefferville+plot_fungi_soil_Fermont+plot_fungi_soil_Schefferville+plot_fungi_rhizosphere_Fermont+plot_fungi_rhizosphere_Schefferville+plot_fungi_root_Fermont+plot_fungi_root_Schefferville+
  plot_layout(ncol = 2, nrow = 6)
