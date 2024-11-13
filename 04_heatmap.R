# Heatmaps

# Libraries 
library(ComplexHeatmap)
library(dplyr)
library(here)

# Load the data
load(here("data", "data-preprocessed.RData"))

# Select 100 most variable genes
tmp <- subset %>%
  select(GENEID, CTRL_1:PRK3_5) %>%
  rowwise() %>%
  mutate(most_variable = sd(c_across(CTRL_1:PRK3_5))) %>%
  arrange(desc(most_variable)) %>%
  ungroup() %>%
  slice_head(n = 100)

tmp.mat <- tmp %>%
  select(CTRL_1:PRK3_5)

# Create a matrix, with genenames as row IDs
tmp.mat <- as.matrix(tmp.mat)
rownames(tmp.mat) <- tmp$GENEID

# Scale the matrix and center
tmp.scaled <- t(scale(t(tmp.mat))) 

# Plot the heatmap
Heatmap(tmp.scaled,  row_names_gp = grid::gpar(fontsize = 4))
