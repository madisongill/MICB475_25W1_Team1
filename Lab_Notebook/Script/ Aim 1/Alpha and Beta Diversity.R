#PPMS Alpha and Beta Diversity 

#Load necessary packages
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(picante)
library(ggplot2)

#Load metadata, feature table, taxonomy and rooted tree files
metaFP <- "metadata_usa.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

otuFP <- "feature-table.txt"
otu <- read.delim(file=otuFP, sep = "\t", skip=1)
colnames(otu)

taxaFP <- "taxonomy.tsv"
taxa <- read.delim(file=taxaFP,sep = "\t")

phylotreeFP <- "tree.nwk"
phylotree <- read.tree(phylotreeFP)
class(phylotree)

# Adjust files to be read into phyloseq object
otu_mat <- as.matrix(otu[,-1])
head(otu_mat)
rownames(otu_mat) <- otu$X.OTU.ID
otu_mat
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

meta$sampleid <- paste0("X", gsub("-", ".", meta$sampleid))
samp_df <- as.data.frame(meta)
rownames(samp_df) <- samp_df$sampleid
samp_df$sampleid <- NULL
SAMP <- sample_data(samp_df)
class(SAMP)

taxa_mat <- taxa |> 
  select(-Confidence) |>
  separate(col=Taxon, sep=";", 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) |>
  as.matrix()
taxa_mat <- taxa_mat[, -1]
rownames(taxa_mat) <- taxa$Feature.ID
taxa_mat
TAXA <- tax_table(taxa_mat)
class(TAXA)

# Create phyloseq object
ms <- phyloseq(OTU, SAMP, TAXA, phylotree)

# Rarefy phyloseq object
ms_rare <- rarefy_even_depth(ms, rngseed = 1, sample.size = 10022)

#Subset the rarefied phyloseq object by MS subtype (RRMS, SPMS, PPMS),
ppms <- subset_samples(ms_rare, disease_course =="PPMS")
ppms

#Calculate Shannon's Diversity index and plot
sd_ppms <- estimate_richness(ppms, measures = "Shannon")
sd_ppms$sex <- sample_data(ppms)$sex

sd_plot <- plot_richness(ppms, x = "sex", measures = "Shannon")

# Box plot
ggplot(sd_ppms, aes(x = sex, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Shannon Diversity by Group",
       x = "sex",
       y = "Shannon ") +
  theme_minimal()

## Calculate Faiths PD 
pd_ppms <- pd(t(otu_table(ppms)),
            phy_tree(ppms),
            include.root = FALSE)
pd_ppms$sex <- sample_data(ppms)$sex

#Boxplot
ggplot(pd_ppms, aes(x = sex, y = PD)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Faith's Phylogenetic Diversity by Group",
       x = "sex",
       y = "PD") +
  theme_minimal()

#Calculate weighted unifrac for beta diversity
dist_uf <- distance(ppms, method = "wunifrac")

pcoa_uf <- ordinate(ppms, method="PCoA", distance = dist_uf)

plot_ordination(ppms, pcoa_uf, color = "sex", shape = "sex")  +
  geom_point(size = 2)

adonis2(dist_uf ~ sex,

        data = data.frame(sample_data(ppms)))
