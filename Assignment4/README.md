## **Assignemnt4**



## **Introduction**

Influenza A virus (IAV) infection remains a major cause of respiratory disease and triggers complex immune responses in host tissues. Following infection, both epithelial and immune cell populations undergo dynamic changes in gene expression, which are critical for viral clearance and tissue recovery (Iwasaki & Pillai, 2014; Krammer et al., 2018). In particular, innate immune cells (macrophages) play a key role in recognizing viral components and initiating inflammatory signaling. On the other hand, epithelial cells contribute to both barrier function and antiviral responses. Understanding how these cell populations change over time is important for interpreting the progression of infection and host defense mechanisms.

Traditional bulk RNA sequencing approaches average gene expression across many cells, making it difficult to resolve cell type–specific responses. In contrast, single cell RNA sequencing (scRNA-seq) allows gene expression to be measured at the individual cell level, enabling identification of heterogeneous cell populations and their transcriptional states (Tang et al., 2009; Luecken & Theis, 2019). This approach is particularly useful in infection studies when different cell types can respond differently to the same stimulus.

Analysis of scRNA-seq data typically involves a series of computational steps, including quality control, normalization, dimensionality reduction, clustering, and marker gene identification. The Seurat framework provides an integrated pipeline for performing these steps, including graph-based clustering and visualization using dimensionality reduction techniques such as UMAP (Butler et al., 2018; Stuart et al., 2019). These methods enable identification of distinct cell populations and comparison of their transcriptional profiles across experimental conditions.

In this study, scRNA-seq data from mouse nasal tissues following IAV infection were analyzed to investigate changes in cell populations and gene expression over time. Clustering and marker gene analysis were used to identify major cell types, followed by differential expression analysis within a selected cluster to compare early and late time points. Functional enrichment analysis was then applied to interpret the biological processes associated with observed transcriptional changes.




## *Mathod

**Data source and preprocessing** 

A pre-processed Seurat object containing single-cell RNA sequencing (scRNA-seq) data and associated metadata was provided for analysis. The metadata included variables such as tissue type and time after infection, and those classifications were used in downstream comparisons. The object was imported into R and inspected to confirm the structure of the expression data and metadata before analysis began.

**Quality control and cell filtering**

Quality control was performed using three metrics: the number of detected genes per cell (nFeature_RNA), the total number of RNA counts per cell (nCount_RNA), and the proportion of mitochondrial transcripts (percent.mt). Mitochondrial percentage was calculated from genes matching the pattern "^mt-". Cells were retained only if they met all of the following criteria: more than 500 but fewer than 4000 detected genes, fewer than 10,000 total RNA counts, and less than 10% mitochondrial expression. Cells outside these thresholds were removed prior to downstream analysis.

Before filtering, several large internal slots in the Seurat object were cleared and garbage collection was run to reduce memory usage during processing. This step was used only for computational efficiency and did not change the analytical workflow.

**Normalization, feature selection, and dimensionality reduction**

The RNA assay was used for downstream analysis. Expression data were normalized using Seurat’s NormalizeData function with default settings. Highly variable genes were then identified using the variance-stabilizing transformation (selection.method = "vst"), with the top 2000 variable features retained for downstream analysis. The data were scaled, and principal component analysis (PCA) was performed using 30 principal components. An elbow plot was generated to visualize the variance structure, and the first 15 principal components were used for neighbor finding, clustering, and UMAP visualization.

**Clustering and visualization**

Cells were clustered using Seurat’s graph-based clustering workflow. A shared nearest-neighbor graph was constructed using the first 15 principal components, and clusters were identified at a resolution of 0.5. Uniform Manifold Approximation and Projection (UMAP) was then performed using the same 15 principal components to visualize the data in two dimensions. UMAP plots were generated to show cluster structure overall, as well as sample grouping by time point and tissue type.

**Marker identification and manual annotation**

To identify cluster-specific marker genes, cluster identities were first set using the Seurat cluster assignments. Because marker testing on the full dataset was more computationally demanding, cells were randomly downsampled so that each cluster contributed at most 100 cells to the marker analysis. Marker genes were identified using FindAllMarkers with the Wilcoxon rank-sum test, restricting results to positive markers with min.pct = 0.25 and logfc.threshold = 0.25. Marker tables were exported, and the top markers for each cluster were selected based on average log2 fold change.

Manual annotation was then performed by examining known marker genes across clusters. Feature plots were generated for selected neuronal, epithelial, immune, stromal, and endothelial markers, including genes such as Cnga2, Omp, Epcam, Ptprc, Nkg7, Lyz2, Col1a1, and Pecam1. Based on these marker patterns, clusters were assigned to cell type labels such as olfactory neurons, macrophages, basal epithelial cells, B cells, dendritic cells, endothelial cells, NK cells, and neutrophils. The annotated cell types were visualized on the UMAP embedding.

**Differential expression analysis**

Differential expression analysis was performed for cluster 2, which was annotated as macrophages in the manual cell type assignment. Cells from this cluster were subsetted, and only cells from the D02 and D14 time points were retained for comparison. To reduce computational load and keep the comparison manageable, each time point was randomly downsampled to a maximum of 200 cells before testing. Differential expression between D14 and D02 macrophages was carried out using Seurat’s FindMarkers function with the Wilcoxon rank-sum test, using min.pct = 0.1 and logfc.threshold = 0.25. The resulting genes were ordered by adjusted p-value and exported for downstream interpretation.

**Functional enrichment analysis**

To examine the biological functions associated with transcriptional changes in macrophages, over-representation analysis (ORA) was performed on genes upregulated in D14 relative to D02. Significant genes were selected from the differential expression results using an adjusted p-value cutoff of 0.1 and an absolute log2 fold change threshold of 0.25. Genes with positive fold change were submitted to Gene Ontology enrichment analysis using the enrichGO function from the clusterProfiler package, with mouse annotation from org.Mm.eg.db, gene symbols as the input key type, Biological Process as the ontology, and Benjamini-Hochberg adjustment for multiple testing. Enrichment results were exported and visualized with a dot plot.



## **Results**

Figure-1

![Model](../Assignment4/scRNA.png)

Figure-1 shows violin plots for the distribution of quality control metrics across all cells. The plot includes the number of detected genes (nFeature_RNA), total RNA counts (nCount_RNA), and percentage of mitochondrial gene expression (percent.mt). These metrics were used to filter low-quality cells prior to downstream analysis.

The QC metrics show a wide distribution of gene counts and sequencing depth across cells. Most cells fall within a moderate range of gene detection and RNA counts, but a small subset displays higher values. Mitochondrial gene expression is generally low across cells. 



Figure-2

![Model](../Assignment4/UMAP_CELL_TYPE.png)


Figure-2 shows the UMAP projection of all cells colored by annotated cell types. Clusters were identified using graph-based clustering and annotated based on known marker gene expression.

The UMAP embedding reveals clear separation of major cell populations, including epithelial cells, immune cells, and neuronal cell types. Distinct clusters corresponding to macrophages, neutrophils, NK cells, and B cells are observed, along with epithelial and olfactory neuron populations. The separation suggests that clustering successfully captured biologically meaningful differences in transcriptional profiles.



FIgure-3

![Model](../Assignment4/UMAP_TIME.png)

Figure-3 is the UMAP projection colored by time point (Naive, D02, D05, D08, D14), showing the distribution of cells across different stages following infection.

Cells from different time points are broadly distributed across the same clusters. Cell identity is the primary driver of clustering. However, subtle shifts in density within certain clusters suggest temporal changes in cell abundance or transcriptional states, further comparison should be applied between early and late time points.



Figure-4

![Model](../Assignment4/expression.png)

Figure-4 shows the expression levels of selected marker genes (Cnga2, Omp, Ptprc, and Epcam) across clusters, used for manual cell type annotation.

Marker gene expression patterns are consistent with expected cell type identities. Neuronal markers such as Cnga2 and Omp are enriched in clusters corresponding to olfactory neurons. Epcam is broadly expressed in epithelial populations. The immune marker Ptprc is highly expressed in immune cell clusters, supporting the annotation of macrophages, lymphocytes, and other immune cell types.


Figure-5

![Model](../Assignment4/ORA_GO.png)

Dot plot showing Gene Ontology (GO) Biological Process enrichment for genes upregulated in macrophages at D14 compared to D02. Dot size represents gene count, and color indicates adjusted p-value.

Enrichment analysis highlights processes related to translation and ribosome biogenesis, including cytoplasmic translation and rRNA processing. These results suggest increased protein synthesis activity in macrophages at the later stage of infection. The results also reflect enhanced cellular activation or functional response.


Figure-6
![Model](../Assignment4/celltype_composition.png)

Figure-6 shows the proportion of each annotated cell type across different time points (Naive, D02, D05, D08, D14).

The relative abundance of cell types varies across time points. Immune cell populations such as macrophages, neutrophils, and NK cells show moderate changes in proportion at later stages (D08–D14) compared to earlier time points. In contrast, epithelial and neuronal populations remain stable across conditions. 


Figure-7
![Model](../Assignment4/Volcano.png)

Volcano plot shows differentially expressed genes between D14 and D02 macrophages. Genes with positive log2 fold change are upregulated at D14, while negative values are upregulated at D02.

Several genes upregulated at D14 are ribosomal proteins (Rps21, Rps28, Rpl37a dominating), indicating increased expression of translation-related genes. In contrast, fewer genes are upregulated at D02, suggesting a stonger transcriptional shift at the later stage. 


## **Discussion**

This study used single-cell RNA sequencing to characterize cell populations and transcriptional changes in mouse nasal tissues following influenza A virus (IAV) infection. Clustering and marker-based annotation identified many epithelial, neuronal, and immune cell types, including macrophages, neutrophils, NK cells, and B cells. The separation of these populations in the UMAP embedding indicates that the clustering approach captured true biological variation across the dataset.

The overall UMAP structure remained stable, suggesting that cell identity is the primary source of variation. However, differences in the distribution of cells across clusters indicate that infection might alter the abundance of specific populations rather than generating entirely new cell types. This is consistent with previous studies that showed the immune responses to IAV involve dynamic transcriptional changes within existing cell populations (Iwasaki & Pillai, 2014; Krammer et al., 2018).

Differential expression analysis of macrophages revealed clear differences between early (D02) and later (D14) stages of infection. Genes upregulated at D14 were enriched for processes related to translation and ribosome biogenesis. This pattern suggests an increase in protein synthesis activity at later stages of infection. One explanation is that macrophages transition to a more active functional state, producing cytokines, signaling molecules, and some immune-related proteins required for sustained antiviral responses. Increased translational activity has been associated with immune cell activation in viral infection contexts, supporting this interpretation (Iwasaki & Pillai, 2014).

However enrichment of translation-related pathways is commonly observed in scRNA-seq datasets and may partly reflect high expression of ribosomal genes rather than a specific biological process (Luecken & Theis, 2019). As a result, the observed enrichment may not exclusively indicate functional activation but could also arise from technical or baseline transcriptional effects. Additional validation or comparison with other cell types would be required to confirm whether this signal is specific to macrophages in this context.

Several limitations of this analysis should also be pointed out. First, differential expression was performed on downsampled cells, which reduces computational burden but that mmight decrease statistical power and limit detection of subtle differences. Second, the analysis focused on a single cluster, meaning that changes in other cell types were not explored. In addition, functional enrichment analysis provides indirect inference of biological processes and thus it does not establish causal relationships. Besides, from the volcano plot, many of the significant genes are ribosomal proteins, and they are highly expressed in scRNA-seq data commonly. This may partially explain the strong enrichment of translational pathways.

Future work could extend this analysis by examining additional clusters, comparing responses across tissue types, or applying trajectory-based approaches to better capture temporal changes in cellular states. Incorporating cell to cell interaction analysis can also provide information into how different populations coordinate during infection.

Overall, this study demonstrates how scRNA-seq can be used to resolve cell type–specific responses to viral infection. The observed increase in translational related processes in macrophages at later stages indicates functional changes, although further analysis is needed to determine the biological significance of this pattern.








