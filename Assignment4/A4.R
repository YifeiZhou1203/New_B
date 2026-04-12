#Assignment 4 scRNA-seq analysis

####Outline#####
#Load data 
#Read Seurat object, check dimensions and metadata 

#Quality control (QC)
#Calculate mitochondrial gene percentage
#Visualize QC metrics

#Filter out low-quality cells
#Preprocessing
#Normalize data
#Select variable genes
#Scale data
#Run PCA

#Clustering & visualization
#Find neighbors and clusters

#Run UMAP
#Plot clusters by different groupings
#Marker detection & annotation
#Identify cluster marker genes 

#Select top markers
#Manually assign cell types
#Focused analysis 
#Subset cluster 2
#Compare D14 vs D02
#Perform differential expression

#Functional analysis
#Filter significant genes
#Run GO enrichment (ORA)

#Cell composition
#Volcano plot of DE genes

#Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(org.Mm.eg.db)

#Set working directory
setwd("~/Desktop/A4")

######Load Seurat object######
seurat_obj <- readRDS("seurat_ass4.rds")

#Qc basic checks
dim(seurat_obj)
head(seurat_obj@meta.data)
table(seurat_obj$organ_custom)
table(seurat_obj$time)
table(seurat_obj$disease__ontology_label, useNA = "ifany")


#######Quality control#########
#selects genes with names start with "mt-"
#stores as a new metadata column
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Remove heavy matrices before filtering to reduce memory use
seurat_obj[["RNA"]]@scale.data <- matrix(0, 0, 0)
seurat_obj@reductions <- list()
seurat_obj@graphs <- list()
seurat_obj@neighbors <- list()
gc()

#Filter cells conditions 
md <- seurat_obj@meta.data
keep_cells <- rownames(md)[
  md$nFeature_RNA > 500 &
    md$nFeature_RNA < 4000 &
    md$nCount_RNA < 10000 &
    md$percent.mt < 10
]

seurat_obj <- seurat_obj[, keep_cells]
gc()

dim(seurat_obj)

######Normalization and clustering######
DefaultAssay(seurat_obj) <- "RNA"

#Normalization
#Selects top 2000 most variable genes
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

ElbowPlot(seurat_obj)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
#graph-based clustering algorithm
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

table(Idents(seurat_obj))

######UMAP plots#######
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "time")
DimPlot(seurat_obj, reduction = "umap", group.by = "organ_custom")


#######Marker discovery for annotation; Downsample to reduce memory load#######

Idents(seurat_obj) <- "seurat_clusters"

set.seed(1234)
#loop and get all cell IDs belonging to cluster cl
#Randomly select cells till 100; Combine all sampled cells from all clusters into one vector

cells_use <- unlist(lapply(levels(Idents(seurat_obj)), function(cl) {
  cs <- WhichCells(seurat_obj, idents = cl)
  sample(cs, min(length(cs), 100))
}))

seu_small <- subset(seurat_obj, cells = cells_use)

#Gene must be expressed in ≥25% of cells in that cluster
cluster_markers_small <- FindAllMarkers(
  seu_small,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(cluster_markers_small, "all_cluster_markers_small.csv", row.names = FALSE)

top10 <- cluster_markers_small %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10, "top10_markers_small.csv", row.names = FALSE)
print(top10, n = 200)

#Feature plots for manual annotation
#Which cluster shows strong expression?
FeaturePlot(seurat_obj, features = c(
  "Cnga2", "Gng8", "Omp",
  "Krt5", "Krt14", "Trp63",
  "Epcam", "Krt8", "Krt18",
  "Ptprc",
  "Cd3d", "Nkg7", "Ms4a1",
  "Lyz2", "Csf1r", "Adgre1",
  "Col1a1", "Dcn",
  "Pecam1", "Kdr", "Emcn"
), ncol = 3)


#####Manual annotation for clusters 0-9########
cluster_annotation <- c(
  "0" = "Olfactory neurons",
  "1" = "Olfactory neurons (subtype)",
  "2" = "Macrophages",
  "3" = "Basal epithelial cells",
  "4" = "B cells",
  "5" = "Immature neurons",
  "6" = "Dendritic cells",
  "7" = "Endothelial cells",
  "8" = "NK cells",
  "9" = "Neutrophils"
)

#convert cluster numbers to cell type labels.
seurat_obj$celltype <- as.character(seurat_obj$seurat_clusters)
seurat_obj$celltype[seurat_obj$celltype %in% names(cluster_annotation)] <-
  cluster_annotation[seurat_obj$celltype[seurat_obj$celltype %in% names(cluster_annotation)]]

DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)


######Differential expression: cluster 2 macrophages, D14 vs D02########
Idents(seurat_obj) <- "seurat_clusters"


#contains only macrophages
cells_cl2 <- WhichCells(seurat_obj, idents = "2")
sub2 <- subset(seurat_obj, cells = cells_cl2)


#time point (D02, D14)
Idents(sub2) <- "time"
cells_use <- WhichCells(sub2, idents = c("D02", "D14"))
sub2 <- subset(sub2, cells = cells_use)

table(sub2$time)

#Cutdown sample for memory safety
set.seed(1234)
cells_d02 <- WhichCells(sub2, idents = "D02")
cells_d14 <- WhichCells(sub2, idents = "D14")

keep_cells <- c(
  sample(cells_d02, min(length(cells_d02), 200)),
  sample(cells_d14, min(length(cells_d14), 200))
)

sub2_small <- subset(sub2, cells = keep_cells)

de_cluster2 <- FindMarkers(
  sub2_small,
  ident.1 = "D14",
  ident.2 = "D02",
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

de_cluster2$gene <- rownames(de_cluster2)
de_cluster2 <- de_cluster2 %>% arrange(p_val_adj)

write.csv(de_cluster2, "de_cluster2_D14_vs_D02.csv", row.names = FALSE)
head(de_cluster2, 20)

#Significant DE genes for enrichment
de_sig2 <- de_cluster2 %>%
  filter(p_val_adj < 0.1 & abs(avg_log2FC) > 0.25)

head(de_sig2, 20)


######ORA: genes upregulated in D14#######
genes_up_D14_2 <- de_sig2 %>%
  filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  unique()

ego_D14_2 <- enrichGO(
  gene = genes_up_D14_2,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ora_results <- as.data.frame(ego_D14_2)
write.csv(ora_results, "ORA_cluster2_D14_up.csv", row.names = FALSE)

head(ora_results, 10)

dotplot(ego_D14_2, showCategory = 10) +
  ggtitle("ORA of genes upregulated in macrophages at D14 vs D02")

saveRDS(seurat_obj, "seurat_obj_processed.rds")



#####Cell type composition across time######
comp_df <- as.data.frame(table(seurat_obj$time, seurat_obj$celltype))
colnames(comp_df) <- c("time", "celltype", "Freq")

#Convert counts to proportions within each time point
comp_df <- comp_df %>%
  group_by(time) %>%
  mutate(Proportion = Freq / sum(Freq))

#Stacked bar plot
ggplot(comp_df, aes(x = time, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Cell type composition across time points",
    x = "Time point",
    y = "Proportion of cells",
    fill = "Cell type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


######Volcano plot for macrophage DE genes######
volcano_df <- de_cluster2

#Add significance category
volcano_df$change <- "Not significant"
volcano_df$change[volcano_df$p_val_adj < 0.1 & volcano_df$avg_log2FC > 0.25] <- "Up in D14"
volcano_df$change[volcano_df$p_val_adj < 0.1 & volcano_df$avg_log2FC < -0.25] <- "Up in D02"

#Add -log10 p-value
volcano_df$neg_log10_padj <- -log10(volcano_df$p_val_adj + 1e-300)

#Select top genes for labeling
top_labels <- volcano_df %>%
  filter(change != "Not significant") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

#plot
ggplot(volcano_df, aes(x = avg_log2FC, y = neg_log10_padj, color = change)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_text(
    data = top_labels,
    aes(label = gene),
    size = 3,
    vjust = -0.5,
    show.legend = FALSE
  ) +
  labs(
    title = "Differentially expressed genes in macrophages (D14 vs D02)",
    x = "Average log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Group"
  ) +
  theme_minimal()
