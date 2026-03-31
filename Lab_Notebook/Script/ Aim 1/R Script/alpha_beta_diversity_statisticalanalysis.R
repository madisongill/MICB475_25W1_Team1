############################
# Load packages
############################

library(phyloseq)
library(tidyverse)
library(picante)
library(vegan)


############################
# Shannon Alpha Diversity
############################

# Calculate Shannon diversity
shannon <- estimate_richness(ps, measures = "Shannon")

# Extract sample IDs
shannon$SampleID <- rownames(shannon)

# Assign MS subtype from sample name
shannon$Subtype <- case_when(
  grepl("PD", shannon$SampleID) ~ "PPMS",
  grepl("ASO", shannon$SampleID) ~ "SPMS",
  grepl("WT", shannon$SampleID) ~ "RRMS"
)

# Pairwise Wilcoxon tests
rrms_vs_spms <- wilcox.test(
  Shannon ~ Subtype,
  data = subset(shannon, Subtype %in% c("RRMS","SPMS"))
)

rrms_vs_ppms <- wilcox.test(
  Shannon ~ Subtype,
  data = subset(shannon, Subtype %in% c("RRMS","PPMS"))
)

spms_vs_ppms <- wilcox.test(
  Shannon ~ Subtype,
  data = subset(shannon, Subtype %in% c("SPMS","PPMS"))
)

# Benjamini-Hochberg correction
pvals <- c(
  rrms_vs_spms$p.value,
  rrms_vs_ppms$p.value,
  spms_vs_ppms$p.value
)

shannon_adj_pvals <- p.adjust(pvals, method = "BH")

# View results
rrms_vs_spms
rrms_vs_ppms
spms_vs_ppms
shannon_adj_pvals


############################
# Faith's Phylogenetic Diversity
############################

# Extract OTU table
otu_mat <- t(as(otu_table(ps), "matrix"))

# Extract phylogenetic tree
tree <- phy_tree(ps)

# Calculate Faith's PD
faith_pd <- pd(otu_mat, tree, include.root = FALSE)

# Convert to dataframe
faith_pd <- as.data.frame(faith_pd)

# Add sample IDs
faith_pd$SampleID <- rownames(faith_pd)

# Assign MS subtype
faith_pd$Subtype <- case_when(
  grepl("PD", faith_pd$SampleID) ~ "PPMS",
  grepl("ASO", faith_pd$SampleID) ~ "SPMS",
  grepl("WT", faith_pd$SampleID) ~ "RRMS"
)

# Wilcoxon tests
rrms_vs_spms_pd <- wilcox.test(
  PD ~ Subtype,
  data = subset(faith_pd, Subtype %in% c("RRMS","SPMS"))
)

rrms_vs_ppms_pd <- wilcox.test(
  PD ~ Subtype,
  data = subset(faith_pd, Subtype %in% c("RRMS","PPMS"))
)

spms_vs_ppms_pd <- wilcox.test(
  PD ~ Subtype,
  data = subset(faith_pd, Subtype %in% c("SPMS","PPMS"))
)

# Benjamini-Hochberg correction
pd_pvals <- c(
  rrms_vs_spms_pd$p.value,
  rrms_vs_ppms_pd$p.value,
  spms_vs_ppms_pd$p.value
)

faith_adj_pvals <- p.adjust(pd_pvals, method = "BH")

# View results
rrms_vs_spms_pd
rrms_vs_ppms_pd
spms_vs_ppms_pd
faith_adj_pvals


############################
# Beta Diversity (Weighted UniFrac)
############################

# Calculate weighted UniFrac distance matrix
unifrac_dist <- UniFrac(ps, weighted = TRUE)

# PERMANOVA test
adonis_result <- adonis2(
  unifrac_dist ~ Subtype,
  data = shannon
)

# View PERMANOVA results
adonis_result
