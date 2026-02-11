# Feb 10, 2026 Meeting Agenda

## To do items from last week

* Filter out the European and South American data from the metadata
* Possibly filter out the [eczema data](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/199481) (there were only 32 eczema cases in the USA cohort, the filtered data uploaded does not include the eczema patients)
  * Should we filter out N/A cases? 
* Ensure that the manifest and the metadata match to prevent errors in the QIIME pipeline
* Determine if separate analyses should be done based on disease progression (PPMS, SPMS, RRMS)
* Outline for the proposal

## Meeting Agenda

### Research Question: Are there differences in gut microbiome composition between female and male MS patients?

#### Disease distribution
* The samples for men and women small for PPMS and SPMS
  *   Particularly for men diagnosed PPMS, is this a significant sample to do analysis on?
  *   Keeping in mind we would have to run 6 separate analyses, not including control
* [Forms of MS](https://pmc.ncbi.nlm.nih.gov/articles/PMC4132635/):
  * PPMS: Steady decline of neurological function without recovery
  * SPMS: Progressive neurological decline
  * RRMS: Initial stage of MS, there are alternating periods of neurological disability and recovery
<img width="129" height="102" alt="image" src="https://github.com/user-attachments/assets/d05ae992-d01a-4c72-a100-0a5fd1bfdb51" />

## Including eczema patients
<img width="120" height="89" alt="image" src="https://github.com/user-attachments/assets/bb7f0fa0-57c3-4760-a796-afc19f7f65a5" />


#### Proposal Outline 
* Title
* Background info
  * Focus on literature describing strong sex-differences in MS patients
    * Sex-bias in diagnosis, progression of disease, and onset of symptoms
  * Draw from literature that focuses on gut microbiome differences in other major disease 
* Hypothesis
* Aims
  * Look at diversity differences
  * Look at composition differences
  * Look at functional differences   
* Approach
* Metrics
  * Diversity metrics (alpha and beta diversity)
  * Taxonomic barplots
  * Core microbiome
  * Indicator taxa
  * Deseq
  * Functional analysis (picrust2)
* Dataset overview
  * Explanation for removing eczema data
  * Explanation for removing patients outside USA cohort 
* Gantt chart
* Participation report


## Action Items 
* Meeting next week? Reading break



## To do items from last week

* Filter out the European and South American data from the metadata
* Possibly filter out the [eczema data](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/199481) (there were only 32 eczema cases in the USA cohort, the filtered data uploaded does not include the eczema patients)
  * Should we filter out N/A cases? 
* Ensure that the manifest and the metadata match to prevent errors in the QIIME pipeline
* Determine if separate analyses should be done based on disease progression (PPMS, SPMS, RRMS)
* Outline for the proposal


## Notes from meeting
Not necessary to filter out eczema patients
Separate the three MS subtypes if we think based on literature search that there are sex based differences
Could group the two male severity subtypes, 8 not too low to work with
RRMS and SPMS group together?
Include eczema group not worth it to remove
Want to balance number of samples per comparison group, if not feasible it is okay just need enough samples per group
A lot of comparison groups will not see significant differences
Manifest because subsetting data
Do not alter actual FASTQ file, keep files you want and not the ones you are not running in analysis own
# For proposal:
* Need QIIME2 pipeline done in proposal
* Do not filter out eczema patients
* Do not group the MS groups suggested
* Mention diversity, we are going to use alpha and beta diversity, picrust 2 functionality, 3 main parts
* Aims had 6 describe we propose to do these analyses and how are they going to answer our research question

* Proposed title
do not want to sound conclusive, because no results obtained
avoid being too broad

*Introduction/background
define key terms, or anything that is not common knowledge
how MS will affect the microbiome as well as sex
does not have to be human data could be mouse models
other autoimmune disorders could be used if cannot find anything on MS
if has been done before note that
T1D sex differences been studied so ensure to include in introduction
Why do you think different between male and female: ex: more inflammed environment, different immune system

* Go from broad to specific introduction

  Why do we think microbiome can answer sex difference in MS

* Hypothesis/Research question
  We hypothesize that there will be significance differences in MS between men and women
  rationales grounded by past research

* Experimental Aims
  6 Aims in total
  could also do 3 aims, diversity, taxonomy and functionality differences

  Beta select 2 of them Weighted unifrac is gold standard, people normally do Bray-curtis answer what it tells you

* Proposed approach
  detail in table format: have aims as columns, rows describe different steps
  Aim 2-5 dependent on QIIME 2 pipeline
  do not be redundant and repeat
 * citations are important cite all tools you use, QIIME2 and DADA2 as it is separate R packages cite it, picrust2
Need Gantt chart, leave time for writing the manuscript, making the slides, do not need to follow too closely

  

## Action Items 
* Meeting next week? Reading break
* QIIME 2 assignment 2: read length, seq depth, rarefaction, justify primers used, truncation
* Make sure to send to Bessie before proceeding
* References just need to be consistent style (in text numbers) reference section use Zotero
* Participation report: how everyone contributed
* Have 4 weeks to do analysis
* 1 person per aim
* For each aim have different folder, scripts go into each individual aim, subfolder for figures and tables

