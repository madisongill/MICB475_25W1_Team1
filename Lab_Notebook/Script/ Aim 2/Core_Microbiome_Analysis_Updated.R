# Core microbiome analysis of female and male MS patients 

# Load necessary libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(BiocManager)
library(ggvenn)

# Load phyloseq object 
ms <- readRDS("phyloseq_object_filtered.rds")

# Subset to only include MS patients
ms_only <- subset_samples(ms_filt_nolow_samps, disease == "MS")

# Convert relative abundance 
ms_RA <- transform_sample_counts(ms_only, fun=function(x) x/sum(x))

# Subset dataset into treatment and control groups
female_stat <- subset_samples(ms_RA, `sex`== "F")
male_stat <- subset_samples(ms_RA, `sex`== "M")                         

# Set the prevalence threshold and abundance thresholds 
female_ASVs <- core_members(female_stat, detection=0.001, prevalence = 0.8)
male_ASVs <- core_members(male_stat, detection=0.001, prevalence = 0.8)

#Convert ASV IDs 
tax_df <- as.data.frame(tax_table(ms_only))

#Remove duplicate genera and NAs
female_genera <- unique(na.omit(tax_df[female_ASVs, "Genus"]))
male_genera   <- unique(na.omit(tax_df[male_ASVs,   "Genus"]))

# Make venn diagram
venn <- ggVennDiagram(x = list(Female = female_genera, Male = male_genera))

ggsave("venn_ms_sex_filtered.png", plot = venn, width = 16, height = 8)

# Create a named list
venn_list <- list(
  Female = female_genera,
  Male   = male_genera)

# Make venn diagram with genus names
venn_genus <- ggvenn(
  venn_list,
  show_elements    = TRUE,     
  label_sep        = "\n",     # one genus per line
  fill_color       = c("#1B4F72", "#56B4E9"),
  stroke_size      = 0.8,
  set_name_size    = 5,
  text_size        = 2.5,      # reduce if names overlap
  show_percentage  = TRUE)

ggsave("venn_ms_sex_genera_named.png", plot = venn, width = 16, height = 12)
