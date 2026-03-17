# Core microbiome analysis of female and male MS patients 

# Load necessary libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(BiocManager)

# Load phyloseq object 
ms <- readRDS("phyloseq_object_filtered.rds")

# Subset to only include MS patients
ms_only <- subset_samples(ms, disease == "MS")

# Convert relative abundance 
ms_RA <- transform_sample_counts(ms_only, fun=function(x) x/sum(x))

# Subset dataset into treatment and control groups
female_stat <- subset_samples(ms_RA, `sex`== "F")
male_stat <- subset_samples(ms_RA, `sex`== "M")                         

# Set the prevalence threshold and abundance thresholds 
female_ASVs <- core_members(female_stat, detection=0.001, prevalence = 0.8)
male_ASVs <- core_members(male_stat, detection=0.001, prevalence = 0.8)

# Make a Venn-diagram
venn <- ggVennDiagram(x=list(Female= female_ASVs, Male = male_ASVs))

# Save venn diagram
ggsave("venn_ms_sex.png", plot = venn, width = 16, height = 8)

# Get taxonomy table
tax <- as.data.frame(tax_table(ms_only))

library(patchwork)
library(gridExtra)

# Clean the genus names (remove g__ prefix if present)
clean_names <- function(x) gsub("^g__", "", x)

female_only_clean <- clean_names(female_only)
shared_clean      <- clean_names(shared)
male_only_clean   <- clean_names(male_only)

# Plain Venn with counts only - no annotations
venn_plain <- ggVennDiagram(
  x = list(Female = female_labels, Male = male_labels),
  label = "count",
  set_size = 5,
  label_size = 5
) +
  scale_fill_gradient(low = "#cce5ff", high = "#4a90d9") +
  labs(title = "Core Microbiome: Female vs Male MS Patients",
       subtitle = "detection = 0.001 | prevalence = 80%") +
  theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "none")

# Build summary table with cleaned names
max_len <- max(length(female_only_clean), length(shared_clean), length(male_only_clean))
pad <- function(x, n) c(x, rep("", n - length(x)))

table_df <- data.frame(
  `Female only (6)` = pad(female_only_clean, max_len),
  `Shared (4)`      = pad(shared_clean,      max_len),
  `Male only (1)`   = pad(male_only_clean,   max_len),
  check.names       = FALSE
)

table_plot <- tableGrob(
  table_df,
  rows = NULL,
  theme = ttheme_minimal(
    base_size     = 12,
    core = list(fg_params = list(fontface = "italic")),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 12))
  )
)

# Combine Venn + table
combined <- venn_plain + wrap_elements(table_plot) +
  plot_layout(widths = c(2, 1.2))

ggsave("venn_ms_sex_taxonomy.png", plot = combined,
       width = 16, height = 8, dpi = 300)
