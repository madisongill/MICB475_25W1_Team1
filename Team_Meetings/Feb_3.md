# Feb 3, 2026 Meeting Agenda

## 1. Potential Research Topics

1) Multiple Sclerosis Dataset (MS)
Variables interested in: Gender, Method of birth, Allergies if a specific allergen type is more likely to develop MS?


Research Question: Compare gender differences of microbiome of: (1) non-MS patients and MS patients, (2) female MS patients and male MS patients


Previous literature:
There have been strong links in gender differences in MS patients:  
* Men diagnosed with MS seem to progress cognitively in symptoms more rapidly
* Women are more likely to be diagnosed with MS
* Observed sex bias in symptoms, diagnosis of MS, disease progression  


References  
[Sex-based differences in MS](https://pmc.ncbi.nlm.nih.gov/articles/PMC8537319/)  
[Men with MS experience accelerated cognitive decline](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1175874/full)  
[Differences in Gut Microbiome in MS vs Healthy Patients](https://pmc.ncbi.nlm.nih.gov/articles/PMC9965298/)  
[MS is more frequent in women: Looking at the oral and gut microbiome](https://www.mdpi.com/2076-3417/13/10/5881)


Notes from meeting:

When look at Male vs Female do within different MS types (there are 3 types, one is severe the other two not as severe)
Confounding variables impact 
Divide dataset into 4 because need controls as well 
Look at Metadata for other confounding variables such as region would need to control for 
Make table summarize: Region, type of MS, gender
Recommend focusing on one region as diet has big impact on microbiome
Filter out Manifest keep everyone in US, get rid of the other regions 

Proposal notes:
Step 1: Control for region (delete and filter out samples not from same region)
Filter out non-US cities in the metadata and manifest files
Only process whatever you filter out

Step 2: QIIME2 Pipeline
Just need to get to table.qzv file for proposal

Step 3:
RRMS biggest sample size least severe
compare female vs male
and then control female vs male
then PPMS female vs male
SPMS female vs male
Bessie suggests looking for n value of each group if not enough samples in severe combine the two

1) Diversity metrics
2) Taxonomic barplot
3) Core microbiome (venn diagram shared vs unique?)
4) Indicator taxa
5) Deseq (bar graph show if rich in one condition vs another)
6) Functional analysis (more involved) uses tool picrust 2, it predicts functional profiles about communities, which pathways are up and down regulated, can say which pathways are pathogenic

Proposal is due in 2 weeks
Bessie will not review draft, but can help build outline

Extract from server don't edit directly on Metadata file


