# Load required libraries
library(tidyverse)

setwd("./GENOSCAN/")

#---------------- Plot MAGS enrichment ----------------------

# Read data
COG20_functions_mod <- read.csv("data/significatly-enriched-COG20_funtion.csv")

# Create a custom ordering based on associated_groups and ascending enrichment_score
COG20_functions_mod <- COG20_functions_mod %>%
  arrange(associated_groups, enrichment_score) %>%
  mutate(Annotation = factor(Annotation, levels = unique(Annotation)))

# Calculate the maximum enrichment score for positioning
max_score <- max(COG20_functions_mod$enrichment_score)

p2 <- COG20_functions_mod %>%
  ggplot(aes(x = Annotation, y = enrichment_score, fill = associated_groups)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("coral", "cornflowerblue")) +
  scale_y_continuous(name = "Enrichment score", 
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(name = "Function") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.box.background = element_blank(),
    legend.box.margin = margin(6, 6, 6, 6),
    axis.text.y.right = element_text(size = 8),
    axis.title.y.right = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1, title = "Groups")) +
  # Add q-values as text
  geom_text(aes(label = sprintf("%.3f", adjusted_q_value), y = max_score * 1.05),
            hjust = 0, size = 3, angle = 0)

# Display the plot
print(p2)
