
#Load necessary packages
library(ggpicrust2)
library(tidyverse)
library(DESeq2)
library(dplyr)

#Load pathway file
abundance <- read.table("path_abun_unstrat.tsv.gz",
                        header = TRUE,
                        sep = "\t",
                        row.names = 1,
                        check.names = FALSE)

#Convert to data frame
abundance_data <-  as.data.frame(abundance, check.names = FALSE)

#Load metadata
meta <- read.delim("metadata_usa.tsv",sep = "\t")

#Filter metadata to only include healthy patients
meta_control <- filter(meta, disease == "Control")

#Filter abundance table to only include samples that are in the filtered metadata
common_samples <- intersect(colnames(abundance_data), meta_control$sampleid)
abundance_control <- abundance_data[, common_samples]
meta_control <- meta_control[meta_control$sampleid %in% common_samples, ]
meta_control <- meta_control[match(common_samples, meta_control$sampleid), ]


#Remove zero only samples
abundance_control <- abundance_control[, colSums(abundance_control != 0) > 0]

#Check if metadata matches samples in abundance data
all(colnames(abundance_control) == meta_control$sampleid)

## Run DESEq ##
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_control, 
                                        metadata = meta_control, 
                                        group = "sex", 
                                        daa_method = "DESeq2")

feature_with_p_0.005 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Make PCA Plot to show no significant difference

pca_plot_control <- pathway_pca(
  abundance = abundance_control,
  metadata = meta_control,
  group = "sex" 
)
