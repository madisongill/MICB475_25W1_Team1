
# Team Meeting Agenda March 17, 2026


## Findings from Aim 2:


* Core Microbiome:
  * ASVs:
    <img width="4800" height="2400" alt="venn_ms_sex" src="https://github.com/user-attachments/assets/ef6b672d-cd63-4eba-a24a-5ff756865b25" />
    * 2 types of bacteroids and blautia are shared (different ASVs)
    * bacteroides and blautia genus also show up again for males --> could be different at the species level but we have not resolved that 
  * Unique Genera (removed duplicates):
    <img width="4800" height="2400" alt="venn_ms_sex_filtered" src="https://github.com/user-attachments/assets/ba16f143-89aa-47d0-b74d-52857a1e68f6" />

    <img width="4800" height="2400" alt="venn_genera_sex" src="https://github.com/user-attachments/assets/01ed8101-92d0-4386-8e7f-8d58cac52ac1" />

* Indicator Taxa:
 *  These results are when I lowered the indicator value to 0.5
 * There were no results when I had an indicator value of 0.7 
<img width="1288" height="59" alt="isa_revised_lower" src="https://github.com/user-attachments/assets/6f778db1-a669-4244-bb0d-5ce038f662bb" /> 

* DESeq2:
  
<img width="2903" height="1703" alt="DESeq_barplot" src="https://github.com/user-attachments/assets/8ac16810-45d2-44c6-9625-b04227d855e2" />

<img width="2903" height="1703" alt="DESeq_volcano_plot" src="https://github.com/user-attachments/assets/8bd32a84-bd38-4711-b41c-e5a3d3f02dfb" />

Found in the literature that women have higher abundance than men of Akkermansia, which is consistent with our results of men having a negative log2fold change (Kim et al., 2019). However, in our proposal we found lots of papers supporting an increased abundance of bacteroides spp. in men which does not match the obtained results as we see a log2fold decrease but then we also see an increase in Bacteroides.2 therefore not conclusive on this genus's distribution between men and women (Kim et al., 2019). For the bacterial genus Alistipes found that male mice had more than female mice when on a high fat diet instead of chow diet;however, it was also found that during puberty girls have more Alistipes species than boys (McGee & Huttenhower, 2021; Yuan et al., 2020). For the genus Tyzzerella I was unable to find sex based differences in abundance and distribution;however, it was found that Tyzzerella nexilis strains are more prevalent in patients with progressive MS (Takewaki et al., 2024). Faecalibacterium play a role in producing anti-inflammatory metabolites including butyrate (Ferreira-Halder et al., 2017). Ruminococcus (UCG-002) plays a role in pro-inflammatory response and we expected there to be a higher abundance in women compared to men which correlate with the data obtained. Found that Barnesiella species can enhance intestinal barrier integrity, and modulate inflammation (Liu et al., 2025). No notable sex differences in the literature for E. shigella. Anaerobutyircum species involved in producing butyrate (Kumari et al., 2021). Thomasclavelia species secrete immunoglobulin A proteases to avoid host immune response (Tran et al., 2025). Found that Lachnospiraceae species are responsible for increased incidence of MS (Yoon et al., 2025). Segatella degrades fiber in gut (Wang et al., 2026).

## Questions about the data
1. What is the best method to present the indicator taxa findings?
2. How to determine when reading the bar plot for DeSeq2 analysis if the ratio is women to men or men to women ?
3. How to interpret the log2fold change?

## References Used
1. Ferreira-Halder, C. V., Faria, A. V. de S., & Andrade, S. S. (2017). Action and function of Faecalibacterium prausnitzii in health and disease. Best Practice & Research Clinical Gastroenterology, 31(6), 643–648. https://doi.org/10.1016/j.bpg.2017.09.011
2. Kim, Y. S., Unno, T., Kim, B.-Y., & Park, M.-S. (2019). Sex Differences in Gut Microbiota. The World Journal of Men’s Health, 38(1). https://doi.org/10.5534/wjmh.190009
3. Kumari, M., Singh, P., Nataraj, B. H., Kokkiligadda, A., Naithani, H., Azmal Ali, S., Behare, Pradip. V., & Nagpal, R. (2021). Fostering next-generation probiotics in human gut by targeted dietary modulation: An emerging perspective. Food Research International, 150, 110716. https://doi.org/10.1016/j.foodres.2021.110716
4. Liu, X., Wang, L., Huang, B., Jiao, Y., Guan, Y., & Nuli, R. (2025). Barnesiella intestinihominis improves gut microbiota disruption and intestinal barrier integrity in mice with impaired glucose regulation. Frontiers in Pharmacology, 16. https://doi.org/10.3389/fphar.2025.1635579
5. McGee, J. S., & Huttenhower, C. (2021). Of Mice and Men and Women: Sexual Dimorphism of the Gut Microbiome. International Journal of Women’s Dermatology. https://doi.org/10.1016/j.ijwd.2021.10.007
6. Takewaki, D., Kiguchi, Y., Masuoka, H., Manu, M. S., Raveney, B. J. E., Narushima, S., Kurokawa, R., Ogata, Y., Kimura, Y., Sato, N., Ozawa, Y., Yagishita, S., Araki, T., Miyake, S., Sato, W., Suda, W., & Yamamura, T. (2024). Tyzzerella nexilis strains enriched in mobile genetic elements are involved in progressive multiple sclerosis. Cell Reports, 43(10), 114785. https://doi.org/10.1016/j.celrep.2024.114785
7. Tran, N., Frenette, A., & Holyoak, T. (2025). Structure of the Thomasclavelia ramosa immunoglobulin A protease reveals a modular and minimizable architecture distinct from other immunoglobulin A proteases. Proceedings of the National Academy of Sciences, 122(35). https://doi.org/10.1073/pnas.2503549122
8. Wang, S., Zhou, T., Wang, X., Zhao, J., & Wang, X. (2026). Bridging the gap:Prevotella/Segatella’s
impact on gut barrier function and advanced cultivation strategies to realize the uses in gut health. Gut Microbes, 18(1). https://doi.org/10.1080/19490976.2026.2638001
9. Yoon, H., Gerdes, L. A., Beigel, F., Sun, Y., Kövilein, J., Wang, J., Kuhlmann, T., Flierl-Hecht, A., Haller, D., Hohlfeld, R., Baranzini, S. E., Wekerle, H., & Peters, A. (2025). Multiple sclerosis and gut microbiota: Lachnospiraceae from the ileum of MS twins trigger MS-like disease in germfree transgenic mice—An unbiased functional study. Proceedings of the National Academy of Sciences, 122(18). https://doi.org/10.1073/pnas.2419689122
10. Yuan, X., Chen, R., Zhang, Y., Lin, X., & Yang, X. (2020). Sexual dimorphism of gut microbiota at different pubertal status. Microbial Cell Factories, 19(1). https://doi.org/10.1186/s12934-020-01412-2


# Notes from Team Meeting

- For core microbiome analysis do not need to resolve to species level
- Do not need to include indicator
- For DESeq2 glom to genus level for both plots
- Male ASVs decreased
- add horizontal line as p-adjusted value cut off
- okay to only focus on specific disease, can do healthy and compare male and female in context of MS
- Next step: run healthy analysis vs MS
- reach different conclusions from the different analyses if they don't match that is fine, diff stats to define if core taxa or not
- Core takes in other factors beyond abundance
- indicator species significant to pass p-value but threshold to low, do not need to report findings from ISA
- can say performed ISA but did not generate significant results

#Action Items
 * Fix DESeq plot
 * Glom DESeq
 * Run all Aim 2 analyses on Healthy subsets
 * Run PiCRUST analysis
 * Research on the results of Aim 2 and Aim 3 
