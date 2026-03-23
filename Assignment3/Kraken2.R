####Assignment3

######install the packages#####
install.packages("tidyverse")
install.packages("vegan")
install.packages("pheatmap")
BiocManager::install("phyloseq")
install.packages("BiocManager")
BiocManager::install("ANCOMBC")
BiocManager::install("microbiome")
library(tidyverse)
library(vegan)
library(ggplot2)
library(pheatmap)
library(phyloseq)
library(ANCOMBC)
library(microbiome)

#read .genus files from Braken
# 35, 36, 38 are omnivore; 37, 39, 40 are vegan
files <- c(
  "SRR8146935.bracken.genus",
  "SRR8146936.bracken.genus",
  "SRR8146937.bracken.genus",
  "SRR8146938.bracken.genus",
  "SRR8146939.bracken.genus",
  "SRR8146940.bracken.genus"
)

sample_names <- c(
  "SRR8146935",
  "SRR8146936",
  "SRR8146937",
  "SRR8146938",
  "SRR8146939",
  "SRR8146940"
)

group <- c(
  "omnivore",
  "omnivore",
  "vegan",
  "omnivore",
  "vegan",
  "vegan"
)

metadata <- data.frame(
  Sample = sample_names,
  Group = group,
  stringsAsFactors = FALSE
)


#generate metadata table with the 2 groups
metadata$Group <- factor(metadata$Group, levels = c("omnivore", "vegan"))
print(metadata)


######Function to read one Bracken file:name, fraction_total_reads######
read_bracken <- function(file, sample_name) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  print(paste("Reading:", file))
  print(colnames(df))
  
# Keep only genus name and relative abundance
  df2 <- df[, c("name", "fraction_total_reads")]
  colnames(df2) <- c("Genus", sample_name)
  
  return(df2)
}


####Read all files using the function#####
bracken_list <- list()

for (i in 1:length(files)) {
  bracken_list[[i]] <- read_bracken(files[i], sample_names[i])
}


####Merge all samples by Genus####

abundance_table <- bracken_list[[1]]

for (i in 2:length(bracken_list)) {
  abundance_table <- full_join(abundance_table, bracken_list[[i]], by = "Genus")
}

# Replace NA with 0 and generate csv
abundance_table[is.na(abundance_table)] <- 0

print(head(abundance_table))

write.csv(abundance_table, "genus_abundance_table.csv", row.names = FALSE)


#####Make matrix for AAlpha and Beta density analysis: rows = samples, columns = genera####

abundance_matrix <- abundance_table %>%
  column_to_rownames("Genus") %>%
  t()

abundance_matrix <- as.data.frame(abundance_matrix)
abundance_matrix$Sample <- rownames(abundance_matrix)

abundance_matrix <- left_join(metadata, abundance_matrix, by = "Sample")


print(dim(abundance_matrix))
print(head(abundance_matrix[, 1:6]))

write.csv(abundance_matrix, "genus_abundance_matrix_with_metadata.csv", row.names = FALSE)


#####Top 10 genera bar plot#####
genus_only <- abundance_matrix[, !(colnames(abundance_matrix) %in% c("Sample", "Group"))]

mean_abundance <- colMeans(genus_only)
top10 <- names(sort(mean_abundance, decreasing = TRUE))[1:10]

plot_data <- abundance_matrix %>%
  select(Sample, Group, all_of(top10)) %>%
  pivot_longer(
    cols = all_of(top10),
    names_to = "Genus",
    values_to = "RelativeAbundance"
  )

ggplot(plot_data, aes(x = Sample, y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Group, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Top 10 genus relative abundance",
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("top10_genus_barplot.png", width = 10, height = 6)



####Alpha diversity######
alpha_input <- genus_only
alpha_input <- as.data.frame(alpha_input)

shannon <- diversity(alpha_input, index = "shannon")
simpson <- diversity(alpha_input, index = "simpson")
observed <- specnumber(alpha_input)

alpha_div <- data.frame(
  Sample = abundance_matrix$Sample,
  Group = abundance_matrix$Group,
  Shannon = shannon,
  Simpson = simpson,
  Observed = observed
)

print("Alpha diversity:")
print(alpha_div)

write.csv(alpha_div, "alpha_diversity.csv", row.names = FALSE)

# Shannon boxplot
ggplot(alpha_div, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 3) +
  theme_bw() +
  labs(title = "Shannon diversity by diet group")

ggsave("shannon_boxplot.png", width = 6, height = 5)



######Beta diversity (PCoA plot)#####

#Apply Bray-Curtis distance
bray <- vegdist(genus_only, method = "bray")

#PCoA
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

#percent variance
eig_vals <- pcoa$eig
percent_var <- round(100 * eig_vals[1:2] / sum(eig_vals[eig_vals > 0]), 1)

pcoa_df <- data.frame(
  Sample = abundance_matrix$Sample,
  Group = abundance_matrix$Group,
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

print(pcoa_df)

#plot
p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(aes(label = Sample), vjust = -1, size = 3.5, show.legend = FALSE) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.68, color = NA, show.legend = FALSE) +
  stat_ellipse(level = 0.68, linewidth = 0.8, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity",
    x = paste0("PCoA1 (", percent_var[1], "%)"),
    y = paste0("PCoA2 (", percent_var[2], "%)")
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

print(p_beta)
ggsave("pcoa_braycurtis.png", p_beta, width = 7.5, height = 5.5, dpi = 300)

#PERMANOVA
permanova <- adonis2(bray ~ Group, data = metadata, permutations = 999)
print(permanova)
capture.output(permanova, file = "permanova_results.txt")


######Heatmap of top 20 genera#####
top20 <- names(sort(mean_abundance, decreasing = TRUE))[1:20]

heatmap_data <- genus_only[, top20, drop = FALSE]

#set sample names as rownames
rownames(heatmap_data) <- abundance_matrix$Sample

#annotation table
annotation_df <- data.frame(Group = abundance_matrix$Group)
rownames(annotation_df) <- abundance_matrix$Sample

# check if they match
print(rownames(heatmap_data))
print(rownames(annotation_df))

pheatmap(
  as.matrix(heatmap_data),
  scale = "row",
  annotation_row = annotation_df,
  main = "Top 20 genera heatmap"
)



#####Differential abundance with ANCOMBC2#####
count_list <- list()

for (i in 1:length(files)) {
  df <- read.delim(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
# Use estimated read counts, not fraction_total_reads
  df2 <- df[, c("name", "new_est_reads")]
  colnames(df2) <- c("Genus", sample_names[i])
  
  count_list[[i]] <- df2
}

count_table_counts <- count_list[[1]]

for (i in 2:length(count_list)) {
  count_table_counts <- full_join(count_table_counts, count_list[[i]], by = "Genus")
}

count_table_counts[is.na(count_table_counts)] <- 0

#matrix for phyloseq and metadata table 
otu_mat <- count_table_counts %>%
  column_to_rownames("Genus") %>%
  as.matrix()

otu_mat <- round(otu_mat)

meta_df <- metadata
rownames(meta_df) <- meta_df$Sample

#make sure sample order matches
otu_mat <- otu_mat[, rownames(meta_df)]

#create simple taxonomy table
tax_mat <- matrix(rownames(otu_mat), ncol = 1)
rownames(tax_mat) <- rownames(otu_mat)
colnames(tax_mat) <- "Genus"

#build phyloseq object
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(meta_df),
  tax_table(tax_mat)
)

print(ps)



####Run ANCOMBC2#####
ancombc.out <- ancombc2(
  data = ps,
  tax_level = "Genus",
  fix_formula = "Group",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0,
  lib_cut = 0,
  s0_perc = 0.05,
  group = "Group",
  struc_zero = TRUE,
  neg_lb = TRUE
)


res <- ancombc.out$res
print(colnames(res))

write.csv(res, "ANCOMBC2_full_results.csv", row.names = TRUE)

#show significant taxa
sig_res <- subset(res, q_Groupvegan < 0.05)
write.csv(sig_res, "ANCOMBC2_significant_results.csv", row.names = TRUE)

print(sig_res)

#get rid of structural zeros 
if (!is.null(ancombc.out$zero_ind)) {
  write.csv(ancombc.out$zero_ind, "ANCOMBC2_structural_zeroes.csv", row.names = TRUE)
  print("Structural zero table saved.")
}


######Volcano-like plot for ANCOMBC2####
plot_df <- res
plot_df$Taxon <- rownames(plot_df)
plot_df$Significant <- plot_df$q_Groupvegan < 0.05

ggplot(plot_df, aes(x = lfc_Groupvegan, y = -log10(p_Groupvegan))) +
  geom_point(aes(color = Significant), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "ANCOMBC2 differential abundance",
    x = "Log fold change (vegan vs omnivore)",
    y = "-log10(p-value)"
  )

ggsave("ANCOMBC2_volcano.png", width = 8, height = 6, dpi = 300)