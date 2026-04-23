#### Load libraries ####

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(dplyr)
library(ggrepel)

#### Load non-rarefied phyloseq object ####
load("ms_phyloseq_norarefy.RData")

# Subset to only include MS patients
ms_only <- subset_samples(ms_filt_nolow_samps, disease == "Control")

# glom to Genus
ms_genus <- tax_glom(ms_only, "Genus", NArm = FALSE)

# Rename the glommed phyloseq object using the genus names
genus_names <- as.character(tax_table(ms_genus)[,"Genus"])
genus_names[is.na(genus_names)] <- "Unknown"
taxa_names(ms_genus) <- make.unique(genus_names)

#### DESeq ####
## Got a zeros error, will add '1' count to all reads
ms_plus1 <- transform_sample_counts(ms_genus, function(x) x+1)
ms_deseq <- phyloseq_to_deseq2(ms_plus1, ~`sex`)
DESEQ_ms <- DESeq(ms_deseq)
res <- results(DESEQ_ms, tidy=TRUE) 
View(res)

res_padj <- res %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2)
view(res_padj)

## Volcano plot: effect size VS significance
ggplot(res_padj) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
volcano_plot <- ggplot(res_padj, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  geom_text_repel(
    data = subset(res_padj, significant),
    aes(label = row),
    size = 4.5
  ) +
  theme_classic() +
  theme(text = element_text (size = 16))
volcano_plot

ggsave(filename="DESeq_volcano_plot_healthy_new.png",volcano_plot)


# To get a table of results from DESEQ_ms
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
ms_phyloseq_norarefy_DESeq <- prune_taxa(sigASVs_vec,ms_genus)
sigASVs <- tax_table(ms_phyloseq_norarefy_DESeq ) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

#Make bar plot 

DESeq_barplot<-ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
DESeq_barplot

# Save bar plot
ggsave(filename="DESeq_barplot_healthy.png",DESeq_barplot)

