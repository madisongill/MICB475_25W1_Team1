# Load packages
library(tidyverse)
library(vegan)

#Load otu table
otu <- read.table("feature-table.txt",
                        header = TRUE,
                        sep = "\t",
                        comment.char = "",
                        skip = 1,
                        row.names = 1,
                        check.names = FALSE)
 
# Transpose OTU table
otu_table <- t(otu)

# Remove NA values if present
sum(is.na(otu_table))

#Ensure data is numeric
otu_table <- as.data.frame(otu_table)
otu_table <- otu_table[rowSums(otu_table) > 0, ]

#Generate rarefaction curve
rarecurve(otu_table,
          step = 1000,
          label = FALSE,
          xlab = "Sequencing depth",
          ylab = "Observed features")
