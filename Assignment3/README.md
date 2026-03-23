## **Assignment3**


## **Introduction**

The human gut microbiome is a complex community of microorganisms that is critical to host metabolism, immune regulation, and overall health. Increasing evidence suggests that diet is one of the most influential factors shaping microbial composition and diversity. In particular, plant-based diets have been associated with increased microbial diversity and enrichment of taxa involved in fiber fermentation, while omnivorous diets may favor different metabolic profiles (David et al., 2014). 

Several sequencing approaches have been used to study microbial communities, mostly 16S rRNA gene sequencing and shotgun metagenomic sequencing. 16S rRNA sequencing is limited in taxonomic resolution and cannot distinguish closely related species or provide functional insights (Callahan et al., 2017). In contrast, shotgun metagenomic sequencing enables comprehensive analysis of all genetic material in a sample, allowing for higher-resolution taxonomic classification and the potential to infer functional capabilities of microbial communities (Quince et al., 2017). 

Within shotgun metagenomic workflows, K-mer–based classifiers such as Kraken2 enable rapid assignment of sequencing reads by matching short DNA sequences to reference databases (Wood et al., 2019), while downstream Bracken refines abundance estimates by redistributing reads assigned to higher taxonomic levels (Lu et al., 2017). Compared to alignment methods, k-mer–based methods are computationally efficient but may be influenced by database completeness and parameter choices. Additionally, the use of reduced reference databases, such as mini Kraken2 databases, improves computational efficiency but may reduce classification accuracy and taxonomic resolution. By analyzing microbial composition at the genus level using shotgun metagenomic data, it is possible to investigate how dietary patterns influence microbiome structure while considering the methodological trade-offs inherent in different analytical approaches.

The objective of this study was to compare gut microbiome composition between omnivore and vegan individuals using shotgun metagenomic data. Specifically, this study aimed to examine taxonomic composition, assess alpha diversity, evaluate beta diversity, and identify differentially abundant taxa between dietary groups.



## **Methods**

**Data and Preprocessing**

Six shotgun metagenomic samples were analyzed in this study, consisting of three omnivore samples (SRR8146935, SRR8146936, SRR8146938) and three vegan samples (SRR8146937, SRR8146939, SRR8146940). The raw data were downloaded in SRA format and converted to FASTQ files using the SRA Toolkit (fasterq-dump : 3.4.0).

Quality control was performed using FastQC (version 0.12.1) for access to sequence quality scores, GC content, sequence length distribution, and the presence of adapter contamination. Minor warnings in GC and base content were observed, but the sequencing data were acceptable for downstream analysis, and no additional trimming or filtering was applied.



**Taxonomic Classification and Abundance Estimation**

Taxonomic classification was performed using Kraken2 (version 2.1.3), which assigns reads by matching sequence k-mers to a reference database. A mini Kraken2 database (~8 GB) was used. 

Paired-end reads from each sample were classified separately using Kraken2. For each sample, both a classification output file and a report file were generated. For example:
kraken2 \
  --db /mnt/e/kraken_db \
  --paired SRR8146935_1.fastq SRR8146935_2.fastq \
  --threads 8 \
  --report SRR8146935.report \
  --output SRR8146935.kraken

--db to specify the database, --paired for paired-end reads, --threads 8 to improve runtime, --report to generate a taxonomic summary, and --output to save read-level classifications. Default classification settings were used.


To improve abundance estimates, Bracken was applied to the Kraken2 report files at the genus level. Bracken redistributes reads assigned to higher taxonomic levels. Bracken was run with -l G, so abundance was estimated at the genus level.
bracken \
  -d /mnt/e/kraken_db \
  -i SRR8146935.report \
  -o SRR8146935.bracken.genus \
  -r 150 \
  -l G

 -d for the database path, -i for the Kraken2 report input, -o for the output file, -r 150 for the 150 bp read length, and -l G for genus-level estimation. The Bracken genus output files were then used for downstream abundance analysis in R.




**Construction of Abundance Table**

Genus abundance data from all six samples were merged into a single abundance table using R. Missing values were replaced with zeros to ensure consistent representation of taxa across samples. The resulting matrix contained samples as rows and genera as columns, and the values represent relative abundance.




**Alpha Diversity Analysis**

Alpha diversity was calculated using the vegan R package (Oksanen et al., 2022). The following diversity indices were computed: Shannon diversity index, Simpson diversity index and Observed richness. These metrics were calculated based on the genus abundance data and compared between omnivore and vegan groups.




**Beta Diversity Analysis**

Beta diversity was assessed using Bray–Curtis dissimilarity. A Bray–Curtis distance matrix was calculated using the vegdist function in the vegan package (Oksanen et al., 2022).

Principal Coordinates Analysis (PCoA) was performed using classical multidimensional scaling to visualize patterns of similarity and dissimilarity among samples. The proportion of variance explained by each principal coordinate was calculated from the eigenvalues and used to label axes in the ordination plot. Group differences in microbial community composition were statistically evaluated using PERMANOVA by adonis2 function with 999 permutations.



**Heatmap Visualization**

To further examine patterns in microbial community structure, a heatmap of the top 20 most abundant genera was generated using the pheatmap package in R (Kolde, 2019). Genera were ranked based on mean relative abundance across all samples, and the top 20 taxa were selected for visualization.




**Differential Abundance Analysis**

Beside pheatmap, differential abundance analysis was performed using ANCOMBC2. Genus count data from Bracken estimated read counts were used as input. phyloseq was constructed to integrate the abundance table with sample metadata. Differential abundance was tested using diet group (omnivore vs vegan) as the main variable.

The ANCOMBC2 model will be using the following parameters:

Fixed effect: diet group
Multiple testing correction: Holm method
Structural zero detection: enabled
Pseudo-count sensitivity analysis: enabled

Plus, genera with adjusted p-values (q-values) below 0.05 were considered statistically significant.


These analyses collectively provided inspections on microbial community structure, including within-sample diversity (alpha diversity), between-sample variation (beta diversity), and taxon-specific differences (differential abundance).






## **Result**
**Taxonomic Composition**

Figure 1. 
![Model](../image/top10_genus_barplot.png)


Figure-1 is the relative abundance of the top 10 bacterial genera across omnivore and vegan samples. Each bar represents an individual sample, and colors indicate different genera. The result revealed variation across samples, with some dominant taxa contributing to the microbial community. Alistipes showed relatively high abundance across all samples, with values from 0.035 in SRR8146936 to 0.249 in SRR8146935.

Differences between dietary groups were observed. Bacteroides showed higher abundance in some omnivore samples (0.209 in SRR8146938) compared to vegan samples (0.0487 in SRR8146940). Similarly, Phocaeicola showed higher mean abundance in omnivore samples (0.0633) compared to vegan samples (0.0254). In contrast, some genera such as Parabacteroides appeared consistent across groups, suggesting both shared core microbiota and diet-associated variation.

Overall, omnivore samples showed stronger dominance by specific taxa, while vegan samples displayed more evenly distributed genus-level abundances.



**Alpha Diversity**

Figure-2
![Model](../image/shannon_boxplot.png)


Figure 2. Shannon diversity index comparison between omnivore and vegan groups. Points represent individual samples, and boxplots summarize group distributions.

Alpha diversity analysis demonstrated that vegan samples exhibited higher Shannon diversity compared to omnivore samples. Specifically, Shannon diversity values for vegan samples were 2.88, 2.91, and approximately 2.90, whereas omnivore samples showed greater variability, with values of 1.45, 2.28, and 2.80.

The lowest diversity was observed in omnivore sample SRR8146936 (Shannon = 1.45), indicating a less diverse and more uneven microbial community. In contrast, vegan samples consistently showed high diversity (all > 2.8), suggesting a more complex and evenly distributed microbiome.

Simpson diversity also showed that vegan samples had high values (0.90–0.92), compared to more variable values in omnivores (0.56–0.87). Observed richness also varied widely, with one omnivore sample showing particularly high richness (838 genera), indicating intra-group variability.

Overall, these results suggest that vegan diets are associated with more consistent and higher microbial diversity, although variability within omnivore samples highlights individual differences.



**Beta Diversity**

Figure-3
![Model](../image/pcoa.png)


Figure 3 is the Principal Coordinates Analysis based on Bray–Curtis dissimilarity, showing differences in microbial community composition between omnivore and vegan samples. Each point represents one sample, and distances between points reflect differences. 

Beta diversity analysis showed partial separation between omnivore and vegan samples. Vegan samples tended to cluster more closely in ordination space, indicating greater similarity in microbial composition within this group.

In contrast, omnivore samples were more dispersed, showing higher variability in microbial communities among individuals. This pattern is consistent with the alpha diversity results. 

The slight overlap between the two groups was observed, indicating that diet is not the only determinant of microbiome composition. In addiiton, PERMANOVA analysis itself did not detect statistically significant differences between groups, and the limitation could be the small sample size (n = 3 per group) and within-group variability.



**Differential Abundance**

Figure-4
![Model](../image/Pheat.png)


Figure 4 shows the Differential abundance analysis of bacterial genera between vegan and omnivore groups. Each point represents a genus, with log fold change indicating differences in relative abundance. The analysis did not identify significant difference wiht tesing correction (FDR ≈ 0.87). For this sample size, no taxa showed differential abundance within dietary groups.

On the other hand, Phocaeicola showed higher mean abundance in omnivore samples (0.0633) compared to vegan samples (0.0254), indicating enrichment in omnivores. Prevotella and Barnesiella also showed reduced abundance in vegan samples.




**Differential Abundance continue**

Figure-5
![Model](../image/ANCOMBC2_volcano.png)


The Volcano plot showed abundance results from ANCOMBC2 analysis. Each point represents a genus, with log fold change indicating differences between vegan and omnivore groups under p value. Similarly, differential abundance analysis using ANCOMBC2 did not see genera significance (q > 0.05). 

ANCOMBC2 identified trends in genus abundance between omnivore and vegan diets. Genera on the positive side of the plot showed higher abundance in vegan samples, whereas those on the negative side were enriched in omnivore samples. However, most taxa were clustered near the center of the plot with low p value, showing weak statistical prove for differential abundance. 


## **Discussion**
The study investigated the dietary patterns on gut microbiome composition using shotgun metagenomic data. Overall, the results suggest that diet influences microbial diversity and community structure, although the magnitude of these effects was actually modest. 

To start with, parameter and database choices may influence the results. The mini Kraken2 database reduced runtime, but also lowered taxonomic sensitivity to some unclassified reads. Overall, Kraken2 classified approximately 36–44% of reads, with the remainder unclassified. These rates were relatively low but still sufficient for downstream genus-level analysis because dominant microbial taxa were recovered across samples, including major groups of interest such as Bacteroides and Alistipes, supported by high read counts. However, the unclassified reads could reduce taxonomic resolution and may have limited detection of low abundance organisms.


In addition, alpha diversity analysis indicated that the vegan samples exhibited higher Shannon diversity compared to omnivore samples. This observation aligns with previous studies demonstrating that the plant diets and dietary fibers had increased microbial diversity (David et al., 2014). Higher diversity is associated with greater ecosystem stability and metabolic flexibility, and therefore vegan diets may support a better gut microbiome. In contrast, the lower and more variable diversity observed in omnivore samples may reflect uneven microbial dominance patterns associated with dietary heterogeneity.

Beta diversity analysis further supported the diet differences since vegan samples were clustered more closely than omnivore samples. This suggests that individuals consuming plant diets may harbor more similar microbial communities. This might be due to some shared dietary components such as fiber and complex carbohydrates. However, the overlap observed between dietary groups indicates that diet alone does not fully determine microbiome composition. This is consistent with previous findings that inter-individual variability, host genetics, and environmental factors also play significant roles in shaping the gut microbiome (De Filippis et al., 2019).

At the taxonomic level, several genera exhibited differences in relative abundance between dietary groups. For example, Bacteroides appeared more abundant in some omnivore samples. Bacteroides is related more to the animal-based diets and protein metabolism.  The finding is consistent with prior studies linking Bacteroides to rich fat and animal protein (Wu et al., 2011). Conversely, genera(Faecalibacterium) with plant polysaccharide degradation and short chain fatty acid production were present across samples and may contribute to beneficial metabolic functions, including butyrate production. (Louis & Flint, 2017).

The pheatmap of the top 20 genera showed visible patterns in relative abundance and suggested partial clustering of samples by diet. These patterns indicate potential diet differences in microbial composition. Despite these patterns, differential abundance analysis using ANCOMBC2 did not identify any statistically significant taxa after testing correction. This unexpected result might be the small sample size (n = 3 per group), which reduces statistical power and limits the ability to detect small differences. Further analysis needs larger cohorts and more stable abundance patterns. The volcano plot of ANCOMBC2 showed a wide distribution of log fold changes, indicating that some taxa may differ between dietary groups, although not significantly. The result is consistent with microbiome studies that the diet microbial shifts are often modest and require large sample sizes to achieve statistical significance (De Filippis et al., 2019). On top of that, pheatmap highlights relative patterns across samples, whereas ANCOMBC2 applies rigorous statistical modeling but should be dealing with a more refined dataset. 

Several method designs may also influence the results. The use of a mini Kraken2 database may reduce taxonomic resolution and sensitivity, affecting downstream analyses. Additionally, the use of relative abundance data rather than absolute counts may obscure differences in total microbial load. Besides, the variability during sample processing and sequencing may contribute to the observed heterogeneity.

In conclusion, this study demonstrates that dietary patterns are associated with differences in gut microbiome diversity and composition, with vegan samples showing higher diversity and more consistent community structure. While differential abundance analysis identifies trends but not statistically significant taxa. The observed clustering patterns in genera suggest some relevant differences, but that requires further investigation. Future studies with larger sample sizes, improved taxonomic resolution, and integration of functional data will be more helpful to understand the relationship between diet and the gut microbiome.





## **Reference**

David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., Ling, A. V., Devlin, A. S., Varma, Y., Fischbach, M. A., Biddinger, S. B., Dutton, R. J., & Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. Nature, 505(7484), 559–563. https://doi.org/10.1038/nature12820

De Filippis, F., Pellegrini, N., Vannini, L., Jeffery, I. B., La Storia, A., Laghi, L., Serrazanetti, D. I., Di Cagno, R., Ferrocino, I., Lazzi, C., Turroni, S., Cocolin, L., Brigidi, P., O’Toole, P. W., & Ercolini, D. (2019). High-level adherence to a Mediterranean diet beneficially impacts the gut microbiota. Cell Host & Microbe, 25(4), 596–606. https://doi.org/10.1016/j.chom.2019.03.007

Louis, P., & Flint, H. J. (2017). Formation of propionate and butyrate by the human colonic microbiota. Environmental Microbiology, 19(1), 29–41. https://doi.org/10.1111/1462-2920.13589

Oksanen, J., Blanchet, F. G., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., Minchin, P. R., O’Hara, R. B., Simpson, G. L., Solymos, P., Stevens, M. H. H., Szoecs, E., & Wagner, H. (2022). vegan: Community Ecology Package (R package version 2.6-4).

Wu, G. D., Chen, J., Hoffmann, C., Bittinger, K., Chen, Y.-Y., Keilbaugh, S. A., Bewtra, M., Knights, D., Walters, W. A., Knight, R., Sinha, R., Gilroy, E., Gupta, K., Baldassano, R., Nessel, L., Li, H., Bushman, F. D., & Lewis, J. D. (2011). Linking long-term dietary patterns with gut microbial enterotypes. Science, 334(6052), 105–108. https://doi.org/10.1126/science.1208344

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104





