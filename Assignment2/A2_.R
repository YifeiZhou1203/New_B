library(tximport)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(dplyr)


samples <- data.frame(
  sample = c("SRR10551665","SRR10551664","SRR10551663",
             "SRR10551662","SRR10551661","SRR10551660",
             "SRR10551659","SRR10551658","SRR10551657"),
  stage = factor(c("early","early","early",
                   "thin","thin","thin",
                   "mature","mature","mature"),
                 levels = c("early","thin","mature"))
)



files <- file.path("quants", samples$sample, "quant.sf")
names(files) <- samples$sample

txi <- tximport(files, type="salmon", txOut=TRUE)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ stage
)

dds <- DESeq(dds)



####PCA and Pheatmap
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "stage")



mat <- assay(vsd)
sampleDists <- dist(t(mat))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(
  sampleDistMatrix,
  annotation_col = as.data.frame(colData(dds)[, "stage", drop = FALSE])
)




res_thin_early <- results(dds, contrast=c("stage","thin","early"))
res_mature_early <- results(dds, contrast=c("stage","mature","early"))
res_mature_thin <- results(dds, contrast=c("stage","mature","thin"))



write.csv(as.data.frame(res_thin_early),
          "results/DE_thin_vs_early.csv")

write.csv(as.data.frame(res_mature_early),
          "results/DE_mature_vs_early.csv")

write.csv(as.data.frame(res_mature_thin),
          "results/DE_mature_vs_thin.csv")



#####ORA
res_mature_early <- read.csv("results/DE_mature_vs_early.csv",
                             row.names=1)


sig_genes <- res_mature_early %>%
  filter(padj < 0.05) %>%
  rownames()


ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)


head(ego)

write.csv(as.data.frame(ego),
          "results/GO_mature_vs_early.csv")

dotplot(ego, showCategory=20)



sum(res_mature_early$padj < 0.05, na.rm=TRUE)
head(rownames(res_mature_early))

rownames(res_mature_early) <- sub("_mRNA$", "", rownames(res_mature_early))
sig_genes <- rownames(res_mature_early)[res_mature_early$padj < 0.05]
ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
head(ego)

dotplot(ego, showCategory=20)
write.csv(as.data.frame(ego),
          "results/GO_mature_vs_early.csv")



######GSEA

#Remove NA values and get a ranked vector 
res_clean <- res_mature_early %>%
  filter(!is.na(stat))


ranked_list <- res_clean$stat
names(ranked_list) <- rownames(res_clean)

#Sort the rank list 
ranked_list <- sort(ranked_list, decreasing = TRUE)

head(ranked_list)

gseaRes <- gseGO(
  geneList     = ranked_list,
  OrgDb        = org.Sc.sgd.db,
  keyType      = "ORF",
  ont          = "BP",
  pAdjustMethod= "BH",
  verbose      = FALSE
)
head(gseaRes@result)

write.csv(gseaRes@result,
          "results/GSEA_mature_vs_early.csv")

dotplot(gseaRes, showCategory=20)
gseaplot2(gseaRes, geneSetID = 1)



summary(res_mature_early)
summary(res_thin_early)
summary(res_mature_thin)
sum(res_mature_early$padj < 0.05, na.rm=TRUE)
sum(res_mature_early$padj < 0.05 & res_mature_early$log2FoldChange > 0, na.rm=TRUE)
sum(res_mature_early$padj < 0.05 & res_mature_early$log2FoldChange < 0, na.rm=TRUE)



rownames(res_thin_early)   <- sub("_mRNA$", "", rownames(res_thin_early))
rownames(res_mature_early) <- sub("_mRNA$", "", rownames(res_mature_early))


####pair-wise analysis 
make_volcano <- function(res, title_text) {
  
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  df <- df %>%
    filter(!is.na(padj))
  
  df$neglog10FDR <- -log10(df$padj)
  
  df$significant <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1,
                           "2-fold change",
                           "Not significant")
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10FDR)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    scale_color_manual(values = c("2-fold change" = "red",
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = title_text,
         x = "log2 Fold Change",
         y = "-log10(FDR)") +
    theme_minimal()
}


volcano_A <- make_volcano(res_thin_early,
                          "Volcano Plot: Thin vs Early")
volcano_A


volcano_B <- make_volcano(res_mature_early,
                          "Volcano Plot: Mature vs Early")
volcano_B
