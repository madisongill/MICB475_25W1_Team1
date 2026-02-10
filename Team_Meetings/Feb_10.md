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

#### Proposal Outline 
* Title
* Background info
  * Focus on literature describing strong sex-differences in MS patients
    * Sex-bias in diagnosis, progression of disease, and onset of symptoms
  * Draw from literature that focuses on gut microbiome differences in other major disease 
* Hypothesis
* Aims
  * Look at taxonomic differences
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
