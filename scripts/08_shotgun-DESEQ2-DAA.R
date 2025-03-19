library(dplyr)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(patchwork)
library(readr)
library(tibble)

# set the working directory
setwd("./GENOSCAN/")

# import data
abundance_file <- "KEGG-KO-GENE-abundance-GPM-short.txt"
class_file <- "KEGG-annotation-for-metabolism.txt"
sample_metadata <- "sample-metadata.txt"
KEGG_GPM <- read_tsv(abundance_file)
KEGG_class <- read.delim(class_file)
smd <- read.delim(sample_metadata)
smd <- column_to_rownames(smd, var = colnames(smd)[1])

# data cleaning, merging and summing
KEGG_class[KEGG_class == ""] <- NA # Convert empty strings to NA
# remove rows with NA's
KEGG_class = subset(KEGG_class, !is.na(accession))
df_main <- merge(KEGG_class, KEGG_GPM, by = "accession")
count_table <- df_main[c(1,7:26)]
count_table <- count_table %>%
  group_by_at(1) %>%  # Group by the first column
  summarise(across(everything(), sum, na.rm = TRUE)) %>%  # Sum the values
  column_to_rownames(var = colnames(count_table)[1])
KEGG_class <- df_main[c(1:6)]
KEGG_class <- KEGG_class %>%
  group_by_at(1) %>%  # Group by the first column
  summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%  # Concatenate the values
  column_to_rownames(var = colnames(KEGG_class)[1])
KEGG_class <- as.matrix(KEGG_class)

# create phyloseq object
gene_counts = otu_table(count_table, taxa_are_rows = TRUE)
gene_class = tax_table(KEGG_class)
SMD = sample_data(smd)

KEGG_table = smdKEGG_table = phyloseq(gene_counts, gene_class, SMD)

# view rank names
rank_names(KEGG_table)

# Function to process and plot metabolism data
process_and_plot <- function(phyloseq_obj, title, show_legend = FALSE) {
  dds.data <- phyloseq_to_deseq2(phyloseq_obj, ~ Ecotype)
  dds <- estimateSizeFactors(dds.data, type = 'poscounts')
  dds$Ecotype <- relevel(dds$Ecotype, ref = "Natural")
  dds <- DESeq(dds)
  
  res <- results(dds)[order(results(dds)$padj, na.last=NA), ]
  significant <- res[res$padj < 0.05, ]
  significant <- cbind(as(significant, "data.frame"), as(tax_table(phyloseq_obj)[rownames(significant), ], "matrix"))
  significant <- subset(significant, !is.na(function_def))
  
  significant_negative <- significant %>%
    filter(log2FoldChange < 0) %>% 
    arrange(desc(log2FoldChange)) %>%  
    slice_head(n = 15) %>%
    mutate(direct = "Natural")
  
  significant_positive <- significant %>%
    filter(log2FoldChange > 0) %>% 
    arrange(desc(log2FoldChange)) %>%  
    slice_head(n = 15) %>%
    mutate(direct = "Disturbed")
  
  combined_significant <- bind_rows(significant_positive, significant_negative)
  combined_significant$gene.code <- factor(combined_significant$gene.code, levels = combined_significant$gene.code)
  combined_significant$direct <- factor(combined_significant$direct, levels = c("Disturbed", "Natural"))
  
  plot <- ggplot(combined_significant, aes(x = gene.code, y = log2FoldChange, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2, position = position_dodge(0.05), color = "black") +
    labs(x = NULL, y = "Log fold change", title = title) +
    scale_fill_manual(name = "Ecotype", values = c("coral","cornflowerblue")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = ifelse(show_legend, "bottom", "none")) +
    coord_flip()
  
  return(plot)
}

# Create plots for each metabolism type
carb_plot <- process_and_plot(subset_taxa(KEGG_table, Level2 == "Carbohydrate metabolism"), "Carbohydrate metabolism")
energy_plot <- process_and_plot(subset_taxa(KEGG_table, Level2 == "Energy metabolism"), "Energy metabolism")
xenobiotics_plot <- process_and_plot(subset_taxa(KEGG_table, Level2 == "Xenobiotics biodegradation and metabolism"), "Xenobiotics biodegradation and metabolism")

# Select PGPGs
PGPG_genes <- paste(c("nifA", "nifB", "nifD", "nifE", "nifF", "nifH", "nifHD1", "nifHD2", "nifJ", "nifK", "nifM", "nifN", "nifQ", "nifS", "nifT", "nifU", "nifV", "nifW", "nifX", "nifZ", "epsE", "epsD", "epsF", "epsH", "epsI", "epsJ", "epsL", "epsM", "epsN", "epsO", "nodA", "nodB", "nodC", "node", "nodF", "nodI", "nodJ", "nodU", "nod", "nodT", "nolM", "noeA", "noeB", "noeC", "noeD", "noeE", "nodN", "nodN_like", "nodO", "nodP", "nodS", "nodS_like", "nodY", "nodZ", "nodV", "nodV_like", "nodW", "nodX", "nodX", "nodO", "sodN", "sodC", "sod3", "lipA", "lipB", "lipL", "lipL2", "lipM", "lplA", "kdpA", "kdpB", "kdpC", "kdpD", "kdpE", "kdpF", "puuA", "puuB", "puuC", "puuD", "puuE", "pup", "trpA", "trpB", "trpC", "trpCF", "trpD", "trpE", "trpEG", "trpG", "trpDG", "trpF", "trpS", "trpR", "iaaM", "entC", "acdS", "sodA", "katG", "gpx", "mycC", "gloA", "arsC", "czcA", "pbpA", "ctrA", "zntA", "merA", "arsC", "arsB", "acr3", "cadA", "czcA", "ctrA", "copB", "pbpA", "merA", "merT", "zntA", "czcD", "pelA", "cdiG", "pgaC"), collapse = "|")

# Subset to PGPGs
PGPGs <- phyloseq::subset_taxa(KEGG_table, grepl(PGPG_genes, gene.code, ignore.case = TRUE))
taxa_names(PGPGs)

# Create PGPG plot with legend
PGPG_plot <- process_and_plot(PGPGs, "PGPG / heavy metals", show_legend = TRUE)

# Plot all using patchwork
final_plot <- carb_plot + energy_plot + xenobiotics_plot + PGPG_plot + plot_layout(ncol = 2)

# Display the final plot
print(final_plot)
