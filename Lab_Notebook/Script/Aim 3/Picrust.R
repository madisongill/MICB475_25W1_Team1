
#Load necessary packages
install.packages("ggpicrust2")
library(ggpicrust2)
install.packages("devtools")
library(devtools)
library(tidyverse)
BiocManager::install('DESeq2')
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

#Filter metadata to only include ms patients
meta_ms <- filter(meta, disease == "MS")

#Filter abundance table to only include samples that are in the filtered metadata
common_samples <- intersect(colnames(abundance_data), meta_ms$sampleid)
abundance_ms <- abundance_data[, common_samples]
meta_ms <- meta_ms[meta_ms$sampleid %in% common_samples, ]
meta_ms <- meta_ms[match(common_samples, meta_ms$sampleid), ]


#Remove zero only samples
abundance_ms <- abundance_ms[, colSums(abundance_ms != 0) > 0]

#Check if metadata matches samples in abundance data
all(colnames(abundance_ms) == meta_ms$sampleid)

## Run DESEq ##
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_ms, 
                                        metadata = meta_ms, 
                                        group = "sex", 
                                        daa_method = "DESeq2")

feature_with_p_0.005 <- abundance_daa_results_df %>% filter(p_values < 0.005)

# Make PCA Plot For MS 

ms_pca_plot <- pathway_pca(
  abundance = abundance_ms,
  metadata = meta_ms,
  group = "sex" 
)
