#February 10th: Team Meeting concluded that filtering out eczema data is unecessary. This adds the patients that had been removed who had received an eczema diagnosis
library(tidyverse)
library(dplyr)

#Load the metadata in and rename column sample.id
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- rename(metadata, sampleid = sample.id)
metadata

#Remove Edinburgh, Buenos Aires, and San Sebastian
str(metadata)
colnames(metadata)
unique(metadata$site_x)
table(metadata$site_x)

remove_site <- c("Edinburgh", "Buenos Aires", "San Sebastian")
metadata_usa <- metadata[ !metadata$site_x %in% remove_site, ]

unique(metadata_usa$site_x)

#Check the different types of MS, see the counts for each group

unique(metadata_usa$disease_course)

disease_type <- (metadata_usa %>% count(disease_course))
print(disease_type)

disease_sex <- table(metadata_usa$disease_course, metadata_usa$sex)
print(disease_sex)

#Export edited metadata tsv file
write_tsv(metadata_usa, "metadata_usa.tsv")

#Edit the manifest data to match
keep_samples <- metadata_usa$sample.id

manifest <- read.table("manifest.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
manifest <- rename(manifest, sampleid = sample.id)
manifest

manifest_usa <- manifest[manifest$`sample.id` %in% keep_samples, ]

write_tsv(manifest_usa, "manifest_usa.tsv")

#Check if the manifest and the metadata match :)
metadata_ids <- metadata_usa$sample.id
manifest_ids <- manifest_usa$sample.id

all(manifest_ids %in% metadata_ids)
