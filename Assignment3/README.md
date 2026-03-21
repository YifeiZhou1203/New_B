## **Introduction

The human gut microbiome is a complex community of microorganisms that is critical to host metabolism, immune regulation, and overall health. Increasing evidence suggests that diet is one of the most influential factors shaping microbial composition and diversity. In particular, plant-based diets have been associated with increased microbial diversity and enrichment of taxa involved in fiber fermentation, while omnivorous diets may favor different metabolic profiles (David et al., 2014). 

Several sequencing approaches have been used to study microbial communities, mostly 16S rRNA gene sequencing and shotgun metagenomic sequencing. 16S rRNA sequencing is limited in taxonomic resolution and cannot distinguish closely related species or provide functional insights (Callahan et al., 2017). In contrast, shotgun metagenomic sequencing enables comprehensive analysis of all genetic material in a sample, allowing for higher-resolution taxonomic classification and the potential to infer functional capabilities of microbial communities (Quince et al., 2017). 

Within shotgun metagenomic workflows, K-mer–based classifiers such as Kraken2 enable rapid assignment of sequencing reads by matching short DNA sequences to reference databases (Wood et al., 2019), while downstream Bracken refine abundance estimates by redistributing reads assigned to higher taxonomic levels (Lu et al., 2017). Compared to alignment-based approaches, k-mer–based methods are computationally efficient but may be influenced by database completeness and parameter choices. Additionally, the use of reduced reference databases, such as mini Kraken2 databases, improves computational efficiency but may reduce classification accuracy and taxonomic resolution.

By analyzing microbial composition at the genus level using shotgun metagenomic data, it is possible to investigate how dietary patterns influence microbiome structure while considering the methodological trade-offs inherent in different analytical approaches.

The objective of this study was to compare gut microbiome composition between omnivore and vegan individuals using shotgun metagenomic data. Specifically, this study aimed to examine taxonomic composition, assess alpha diversity, evaluate beta diversity, and identify differentially abundant taxa between dietary groups.


## **Methods

**Data and Preprocessing

Six shotgun metagenomic samples were analyzed in this study, consisting of three omnivore samples (SRR8146935, SRR8146936, SRR8146938) and three vegan samples (SRR8146937, SRR8146939, SRR8146940). The raw data were downloaded in SRA format and converted to FASTQ files using the SRA Toolkit.

Quality control was performed using FastQC (version 0.11.x). FastQC was used to assess sequence quality scores, GC content, sequence length distribution, and the presence of adapter contamination. The QC results indicated that the sequencing data were of acceptable quality for downstream analysis, and no additional trimming or filtering steps were applied.



##Taxonomic Classification and Abundance Estimation

Taxonomic classification was performed using Kraken2 (version 2.x), a k-mer–based classifier that assigns sequencing reads to taxa by exact matching of k-mers to a reference database (Wood et al., 2019). In this analysis, a mini Kraken2 database(8gb) was used rather than the full standard database(50gb-100gb). This reduced database decreases memory usage but may also reduce classification sensitivity and taxonomic resolution.

Default Kraken2 parameters were used for classification, including the default confidence threshold. Reads were classified at multiple taxonomic levels, and genus-level assignments were extracted for downstream analysis.

Following classification, abundance estimation was refined using Bracken (version 2.x), which applies a Bayesian approach to redistribute reads assigned to higher taxonomic levels (Lu et al., 2017). Bracken was run at the genus level, and the output files (.bracken.genus) provided both estimated read counts and relative abundance values.

The parameter settings for Bracken included:
Taxonomic level: genus
Read length: matched to sequencing data (default setting)
Database: corresponding mini Kraken2 database



**Construction of Abundance Table

Genus abundance data from all six samples were merged into a single abundance table using R. Missing values were replaced with zeros to ensure consistent representation of taxa across samples. The resulting matrix contained samples as rows and genera as columns, with values representing relative abundance.



**Alpha Diversity Analysis

Alpha diversity was calculated using the vegan R package (Oksanen et al., 2022). The following diversity indices were computed: Shannon diversity index, Simpson diversity index and Observed richness. These metrics were calculated based on the genus abundance data and compared between omnivore and vegan groups.
