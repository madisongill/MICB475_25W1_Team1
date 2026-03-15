# Core microbiome analysis of female and male MS patients 

# Load necessary libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)
library(BiocManager)

# Load phyloseq object 
ms <- readRDS("phyloseq_object.rds")

# Convert relative abundance 
ms_RA <- transform_sample_counts(ms, fun=function(x) x/sum(x))

# Subset dataset into treatment and control groups
ms_stat <- subset_samples(ms_RA, `disease`=="MS")
healthy_stat <- subset_samples(ms_RA, `disease`=="Control")                         

# Set the prevalence threshold and abundance thresholds 
ms_ASVs <- core_members(ms_stat, detection=0.001, prevalence = 0.8)
healthy_ASVs <- core_members(healthy_stat, detection=0.001, prevalence = 0.8)

# Make a Venn-diagram
venn <- ggVennDiagram(x=list(MS = ms_ASVs, Control = healthy_ASVs))

# Save venn diagram
ggsave("venn_ms_healthy.png", venn)

# Make a core microbiome heatmap
ms_core <- core(ms_RA,
                detection = 0.001,
                prevalence = 0.8)

plot_heatmap(ms_core,
             method = "NMDS",
             distance = "bray",
             sample.label = "disease_course",
             taxa.label = "Genus")

