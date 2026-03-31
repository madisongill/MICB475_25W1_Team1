#### Load in all packages ####

library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(picante)
library(ggplot2)

#### Load in data ####
metaload <- "metadata_usa.tsv"
meta <- read_delim(metaload, delim="\t")

featureload <- "feature-table.txt"
feature <- read_delim(file = featureload, delim="\t", skip=1)

taxload <- "taxonomy.tsv"
tax <- read_delim(taxload, delim="\t")

phylotreeload <- "tree.nwk"
phylotree <- read.tree(phylotreeload)

#### Adjust feature table for phyloseq object ####
feature_matrix <- as.matrix(feature[,-1])

rownames(feature_matrix) <- feature$`#OTU ID`

FEATURE <- otu_table(feature_matrix, taxa_are_rows = TRUE) 
class(FEATURE)

#### Adjust sample metadata for phyloseq object ####

meta_new_df <- as.data.frame(meta[,-1])

rownames(meta_new_df)<- meta$'sampleid'

META <- sample_data(meta_new_df)
class(META)

#### Adjust taxonomy file for phyloseq object ####

tax_new <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 

tax_new <- tax_new[,-1]

rownames(tax_new) <- tax$`Feature ID`

TAX <- tax_table(tax_new)
class(TAX)

#### Making the phyloseq object ####
ms_phyloseq <- phyloseq(FEATURE, META, TAX, phylotree)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
ms_filt <- subset_taxa(ms_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
ms_filt_nolow <- filter_taxa(ms_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
ms_filt_nolow_samps <- prune_samples(sample_sums(ms_filt_nolow)>100, ms_filt_nolow)

# Remove samples where disease course (diagnosis) is na

rrms_disease_course <- subset_samples(ms_filt_nolow_samps, !is.na(disease_course) )

# Remove samples where sex is na
rrms_clean <- subset_samples(rrms_disease_course, !is.na(sex) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rrms_rare <- rarefy_even_depth(rrms_clean, rngseed = 1, sample.size = 10022)

##### Saving #####
save(rrms_clean, file="rrms_final.RData")
save(rrms_rare, file="rrms_rare.RData")

#### Alpha diversity ######
#filtered the rareified object to only include RRMS subtype

rrms_alone <-subset_samples(rrms_rare, disease_course =="RRMS") 

plot_richness(rrms_alone) 

plot_richness(rrms_alone, measures = c("Shannon")) 

estimate_richness(rrms_alone)

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(rrms_alone)), phy_tree(rrms_alone),
                 include.root=F) 

# add PD to metadata table
sample_data(rrms_alone)$PD <- phylo_dist$PD

# plot male vs female for Shannon
gg_Shannon <- plot_richness(rrms_alone, x = "sex", measures = c("Shannon")) +
  xlab("Sex") +
  geom_boxplot()
gg_Shannon

#### Save the Shannon figure as .png file #####
ggsave(filename = "rrms_shannon.png"
       , gg_Shannon
       , height=4, width=6)

# plot male vs female for Faith's phylogenetic diversity

plot.pd <- ggplot(sample_data(rrms_alone), aes(sex, PD)) + 
  geom_boxplot() +
  xlab("Sex") +
  ylab("Phylogenetic Diversity")
plot.pd

#### Save the Faith's phylogenetic diversity figure as .png file #####
ggsave(filename = "rrms_pd.png"
       , plot.pd
       , height=4, width=6)

#### Beta diversity #####
#### Produce an ordination figure #####
ord_rrms <- distance(rrms_alone, method="unifrac")

rrms_ord <- ordinate(rrms_alone, method="PCoA", ord_rrms)

rrms_ord_plot <- plot_ordination(rrms_alone, rrms_ord, color = "sex") +
  labs(color = "Sex") +
  scale_color_manual(values=c("F"="pink", "M"="blue"))
rrms_ord_plot

#### Save the ordination figure as .png file #####
ggsave(filename = "rrms_pcoa.png"
       , rrms_ord_plot
       , height=4, width=6)
