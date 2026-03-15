#Load in libraries 
library(phyloseq)
library(tidyverse)
library(indicspecies)

#Load phyloseq object
ms <-readRDS("phyloseq_object.rds")

#subset data into treatment and control groups
ms_stat <- subset_samples(ms, disease=="MS")
healthy_stat<-subset_samples(ms, disease=="Control")

# Run indicator species analysis using the sex variable
isa_ms<- multipatt(t(otu_table(ms_stat)), cluster = sample_data(ms_stat)$'sex')

summary (isa_ms)

#Extract taxonomy table 
taxtable <- tax_table(ms_stat) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#merge taxonomy table with phyloseq object and filter by significant p-value and stat value
isa_results<- isa_ms$sign %>% 
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05 & stat >0.7)

#View results, but results return nothing. Stat value may be too high. 
View (isa_results)