# Calculating the alpha and beta diversity for females vs males with MS

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

#Filter phyloseq object
ms_filt <- subset_taxa(ms,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
ms_filt_nolow <- filter_taxa(ms_filt, function(x) sum(x)>5, prune = TRUE)
ms_filt_nolow_samps <- prune_samples(sample_sums(ms_filt_nolow)>100, ms_filt_nolow)

# Rarefy phyloseq object
ms_rare <- rarefy_even_depth(ms_filt_nolow_samps, rngseed = 1, sample.size = 10022)

#Subset the rarefied phyloseq object to only include MS patients 
ms <- subset_samples(ms_rare, disease == "MS")

#Calculate Shannon's Diversity index and plot
sd_ms <- estimate_richness(ms, measures = "Shannon")
sd_ms$sex <- sample_data(ms)$sex

sd_plot <- plot_richness(ms, x = "sex", measures = "Shannon")

ggplot(sd_ms, aes(x = sex, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Shannon Diversity by Sex for MS patients",
       x = "sex",
       y = "Shannon ") +
  theme_minimal()

wilcox.test(Shannon ~ sex, data = sd_ms)

## Calculate Faiths PD 
pd_ms <- pd(t(otu_table(ms)),
              phy_tree(ms),
              include.root = FALSE)
pd_ms$sex <- sample_data(ms)$sex

ggplot(pd_ms, aes(x = sex, y = PD)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Faith's Phylogenetic Diversity by Sex for MS patients",
       x = "sex",
       y = "PD") +
  theme_minimal()

wilcox.test(PD ~ sex, data = pd_ms)

#Calculate weighted unifrac for beta diversity
dist_uf <- distance(ms, method = "wunifrac")

pcoa_uf <- ordinate(ms, method="PCoA", distance = dist_uf)

plot_ordination(ms, pcoa_uf, color = "sex", shape = "sex", title = "Weighted Unifrac")  + 
  geom_point(size = 2)

adonis2(dist_uf ~ sex,
        data = data.frame(sample_data(ms)))

plot_ordination(ms, pcoa_uf, color = "sex", shape = "allergies")  + 
  geom_point(size = 2)

# Unweighted Unifrac
unifrac_dist <- distance(ms, method = "unifrac", weighted = FALSE)
ordu_unifrac <- ordinate(ms, method = "PCoA", distance = unifrac_dist)
plot_ordination(ms, ordu_unifrac, color = "sex", title = "Unweighted Unifrac") +
  geom_point(size = 2)

adonis2(unifrac_dist ~ sex, data = data.frame(sample_data(ms)))

#Bray-Curtis Test
bray_dist <- distance(ms, method = "bray")

ordu <- ordinate(ms, method = "PCoA", distance = bray_dist)
plot_ordination(ms, ordu, color = "sex", title = "Bray Curtis") +
  geom_point(size = 2)

adonis2(bray_dist ~ sex, data = data.frame(sample_data(ms)))

#Jacard test
jac_dist <- distance(ms, method = "jaccard")
ordu_jac <- ordinate(ms, method = "PCoA", distance = jac_dist)
plot_ordination(ms, ordu_jac, color = "sex", title = "Jacard") +
  geom_point(size = 2)

adonis2(jac_dist ~ sex, data = data.frame(sample_data(ms)))
