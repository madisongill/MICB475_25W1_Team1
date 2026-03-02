
#Feb 26, 2026 MS Subtype SPMS Aim 1 Workflow

#1. Loading the needed packages
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)

#2 Loading the data

meta_usa_filtered<- "phyloseq_export_files/metadata_usa.tsv"
meta<- read_delim(meta_usa_filtered, delim="\t")

otutable<-"phyloseq_export_files/feature-table.txt"
otu<-read_delim(file= otutable, delim="\t", skip=1)

taxonomy<-"phyloseq_export_files/taxonomy.tsv"
tax<-read_delim(taxonomy, delim="\t")

phylotreefp<-"phyloseq_export_files/tree.nwk"
phylotree<-read.tree(phylotreefp)


#3. Format OTU table 
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#4. Format metadata_usa.tsv

# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sampleid'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#5. Format taxonomy table

# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)


#6. Create phyloseq object 
# Merge all into a phyloseq object
spms <- phyloseq(OTU, SAMP, TAX, phylotree)


#7. Ensuring the data is ready for diversity analysis by cleaning/quality control, removing low reads, contaminated samples, etc
# Remove non-bacterial sequences, if any
spms_filt <- subset_taxa(spms,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
spms_filt_nolow <- filter_taxa(spms, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
spms_filt_nolow_samps <- prune_samples(sample_sums(spms_filt_nolow)>100,spms_filt_nolow)

# Remove samples where disease course (diagnosis) is na
spms_clean_disease_course <- subset_samples(spms_filt_nolow_samps, !is.na(disease_course) )

# Remove samples where sex is na
spms_clean <- subset_samples(spms_clean_disease_course, !is.na(sex) )

# 8. Rarefy samples to 10,022 based on proposal information
# rngseed sets a random number. Set rngseed to 1
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(spms_clean))), cex=0.1)
spms_rare <- rarefy_even_depth(spms_clean, rngseed = 1, sample.size = 10022)

#9. Save data created
save(spms_clean, file="spms_clean.RData")
save(spms_rare, file="spms_rare.RData")


#10. Alpha diversity Shannon Diversity on SPMS 

spms_only<-subset_samples(spms_rare, disease_course =="SPMS") #filtered the rareified object to only include SPMS subtype

plot_richness(spms_only) 

plot_richness(spms_only, measures = c("Shannon")) 

estimate_richness(spms_only)

gg_richness <- plot_richness(spms_only, x = "sex", measures = c("Shannon")) +
  xlab("Sex") +
  geom_boxplot()
gg_richness

ggsave(filename = "spms_shannon_plot_richness.png"
       , gg_richness
       , height=4, width=6)


#11. Feb 27,2026 Tasks
#Faith's phylogenetic diversity

phylo_dist <- pd(t(otu_table(spms_only)), phy_tree(spms_only),
                 include.root=F) 


# add PD to metadata table
sample_data(spms_only)$PD <- phylo_dist$PD

# plot sex against the PD
plot.pd <- ggplot(sample_data(spms_only), aes(sex, PD)) + 
  geom_boxplot() +
  xlab("Sex") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd

# Save the Faith's phylogenetic diversity figure as .png file 
ggsave(filename = "spms_pd.png"
       , plot.pd
       , height=4, width=6)

#12. Beta diversity making an ordination figure using PCoA method and unifrac distance metric

unifrac_dm <- distance(spms_only, method="unifrac")

pcoa_unifrac <- ordinate(spms_only, method="PCoA", distance=unifrac_dm)

plot_ordination(spms_only, pcoa_unifrac, color = "sex")


gg_pcoa <- plot_ordination(spms_only, pcoa_unifrac, color = "sex") +
  labs(col = "Sex")+
  scale_color_manual(values=c("F"="pink", "M"="blue"))
gg_pcoa

#Save the plot as .png file into RProjects folder

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)




