# Core microbiome analysis of female and male control patients 

# Load necessary libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(BiocManager)
library(ggvenn)

# Load phyloseq object 
ms <- readRDS("phyloseq_object_filtered.rds")

# Subset to only include helathy patients
control_only <- subset_samples(ms, disease == "Control")

# Convert relative abundance 
control_RA <- transform_sample_counts(control_only, fun=function(x) x/sum(x))

# Subset dataset into treatment and control groups
f_control_stat <- subset_samples(control_RA, `sex`== "F")
m_control_stat <- subset_samples(control_RA, `sex`== "M")                         

# Set the prevalence threshold and abundance thresholds 
f_control_ASVs <- core_members(f_control_stat, detection=0.001, prevalence = 0.8)
m_control_ASVs <- core_members(m_control_stat, detection=0.001, prevalence = 0.8)

#Convert ASV IDs 
tax_df <- as.data.frame(tax_table(ms_only))

#Remove duplicate genera and NAs
female_genera <- unique(na.omit(tax_df[f_control_ASVs, "Genus"]))
male_genera   <- unique(na.omit(tax_df[m_control_ASVs,   "Genus"]))

# Make venn diagram
venn_control <- ggVennDiagram(x = list(Female = female_genera, Male = male_genera))

ggsave("venn_control_sex_filtered.png", plot = venn_control, width = 16, height = 8)

# Create a named list
venn_list <- list(
  Female = female_genera,
  Male   = male_genera)

print(venn_list)

# Make venn diagram with genus names
venn_genus_control <- ggvenn(
  venn_list,
  show_elements    = TRUE,     
  label_sep        = "\n",     # one genus per line
  fill_color       = c("#1B4F72", "#56B4E9"),
  stroke_size      = 0.6,
  set_name_size    = 5,
  text_size        = 2.5,      # reduce if names overlap
  show_percentage  = TRUE)

venn_genus_control

ggsave("venn_control_sex_genera_named.png", plot = venn_genus_control, width = 16, height = 12)

