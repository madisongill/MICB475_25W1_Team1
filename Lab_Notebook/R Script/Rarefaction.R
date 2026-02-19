# Load packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(phyloseq)

#Load otu table
otu <- read.table("feature-table.txt",
                        header = TRUE,
                        sep = "\t",
                        comment.char = "",
                        skip = 1,
                        row.names = 1,
                        check.names = FALSE)

# Make otu table into phyloseq object --> IF DOING IT AS PHYLOSEQ
otu_mat <- as.matrix(otu[,-1])
head(otu_mat)
rownames(otu_mat) <- otu$X.OTU.ID
otu_mat
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)
 
# Transpose OTU table
otu_table <- t(otu)

# Remove NA values if present
sum(is.na(otu_table))

#Ensure data is numeric
otu_table <- as.data.frame(otu_table)
otu_table <- otu_table[rowSums(otu_table) > 0, ]

rarecurve(otu_table,
          step = 1000,
          label = FALSE)
