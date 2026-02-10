library(tidyverse)


#Load the metadata in
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Remove Edinburgh, Buenos Aires, and San Sebastian
str(metadata)
colnames(metadata)
unique(metadata$site_x)
table(metadata$site_x)

remove_site <- c("Edinburgh", "Buenos Aires", "San Sebastian")
metadata_usa <- metadata[ !metadata$site_x %in% remove_site, ]

unique(metadata_usa$site_x)

#Eczema has previously been found to be a possible confounding variable influencing the gut microbiome, checking if we can remove if cases are smallex
eczema_usa<- metadata_usa %>%
  mutate(eczema = ifelse(eczema == 1, "Yes", "No"))
eczema_usa
table(eczema_usa$eczema)
eczema_yes <- eczema_usa %>%
  filter(eczema == "Yes") %>%
  select(sample.id, site_x)
eczema_yes
table(eczema_yes$site_x)

remove_eczema <- c("1")
metadata_filtered <- metadata_usa[ !metadata_usa$eczema %in% remove_eczema, ]

unique(metadata_filtered$eczema)

#Check the different types of MS, see the counts for each group
unique(metadata_filtered$disease_course)

disease_type <- (metadata_filtered %>% count(disease_course))
print(disease_type)

disease_sex <- table(metadata_filtered$disease_course, metadata_filtered$sex)
print(disease_sex)

#Export edited metadata tsv file

library(readr)

write_tsv(metadata_filtered, "metadatafiltered_usa.tsv")

#Edit the manifest data to match
keep_samples <- metadata_filtered$sample.id

manifest <- read.table("manifest.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

manifest_filtered <- manifest[manifest$`sample.id` %in% keep_samples, ]

write_tsv(manifest_filtered, "manifestfiltered_usa.tsv")

#Check if the manifest and the metadata match :)
metadata_ids <- metadata_filtered$sample.id
manifest_ids <- manifest_filtered$sample.id

all(manifest_ids %in% metadata_ids)
