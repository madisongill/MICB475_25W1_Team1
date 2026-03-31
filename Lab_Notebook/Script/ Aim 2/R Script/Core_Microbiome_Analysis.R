# Core microbiome analysis of female and male MS patients 

# Load necessary libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(BiocManager)

# Load phyloseq object 
ms <- readRDS("phyloseq_object.rds")

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
