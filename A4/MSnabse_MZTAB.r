# BINF 6210 – Assignment 4
# Pipeline for MSnbase aggregation and normalization comparison 
# Author: Yifei Zhou
# Date: 2025-12-05
##############################################################


# Load required packages
library("MSnbase")       
library("reshape2")      
library("ggplot2")       
library("RColorBrewer")
library("rpx")

######1 -Overview----
#1. load ProteomeXchange dataset
#2. Read peptide-level quantification data and Assign TMT reporter channels
#3. Inspect peptide expression data
#4. Filter and preprocess peptide data
#5. Aggregate peptides to protein level, and plot protein-level intensities
#combineFeatures() using sum, median, or mean
#Compare different aggregation methods
#6. Normalize peptide intensities (DIV_median versus VSN)
#7. Select subset of proteins of interest, and prepare data for ggplot2
#8. Melt data to long format; Visualize peptide-level trends per protein



dir.create("figs",    showWarnings = FALSE, recursive = TRUE)


# ######2. -Load a ProteomeXchange dataset PXD000001 using PXDataset-----
px1 <- PXDataset("PXD000001") 

px1  

# List all files available in the PXDataset
pxfiles(px1)  

# Download and access a specific mzTab file from the dataset
mztab <- pxget(px1, "F063721.dat-mztab.txt")  

# "F063721.dat-mztab.txt" is a processed file containing peptide and protein quantification data

mztab  # Print the path to the downloaded mzTab file

# Read peptide-level quantification data from the mzTab file
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")  
# "PEP" is the  peptide-level function 


#######3 -Assign sample names based on TMT6 reporter channels-----
sampleNames(qnt) <- reporterNames(TMT6)  
# Reporter channels for a 6-plex TMT experiment are assigned to columns of qnt

# Inspect the first few rows of the expression matrix
head(exprs(qnt))  



#######4 -Filter and process peptide-level data -----
qnt <- filterNA(qnt)      
processingData(qnt)       


#######5 -Apply further processing: transform, normalization----
# Protein aggregation comparison: sum vs median vs mean 
# Combine peptide-level data into protein-level data using different aggregation methods
# Sum of all peptides for each protein
prot_sum <- combineFeatures(qnt,
                            groupBy = fData(qnt)$accession,  # Group by protein accession
                            method = "sum")                  # Sum intensities of peptides

# Median of peptides for each protein
prot_med <- combineFeatures(qnt,
                            groupBy = fData(qnt)$accession,
                            method = "median")

# Mean of peptides for each protein
prot_mean <- combineFeatures(qnt,
                             groupBy = fData(qnt)$accession,
                             method = "mean")

#######6 -Line plot for last 5 proteins----- 
cls <- brewer.pal(5, "Set1")  # Choose 5 distinct colors from Set1 palette
#par(mfrow = c(1, 3)) 

# Plot protein intensity using sum aggregation
sum_plot <-matplot(t(tail(exprs(prot_sum), 5)), type = "b",   # Transpose matrix and plot last 5 proteins
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(prot_sum), 5),  # Add legend for the last 5 proteins
       lty = 1, bty = "n", cex = .8, col = cls)
ggsave("figs/sum_plot.png", sum_plot)

# Plot protein intensity using median aggregation
med_plot <-matplot(t(tail(exprs(prot_med), 5)), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (median peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(prot_med), 5),
       lty = 1, bty = "n", cex = .8, col = cls)
ggsave("figs/med_plot.png", med_plot)

# Plot protein intensity using mean aggregation
mean_plot <-matplot(t(tail(exprs(prot_mean), 5)), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (mean peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(prot_med), 5),
       lty = 1, bty = "n", cex = .8, col = cls)
ggsave("figs/mean_plot.png", mean_plot)


# ######7 -Normalize using median ---each peptide intensity is divided by the median across samples----
qntM <- normalise(qnt, "div.median")  # Median normalization
qntV <- normalise(qnt, "vsn")  

#Select small subset 
acc <- c("P00489", "P00924", "P02769", "P62894", "ECA")  # List of protein accessions of interest

# Extract indices of peptides belonging to each protein
idx <- sapply(acc, grep, fData(qntM)$accession)
idx2 <- sapply(idx, head, 3)           # Select first 3 peptides for each protein
smallM <- qntM[unlist(idx2), ]         # Subset normalized data

# Repeat for second normalized dataset
idx_vsn <- sapply(acc, grep, fData(qntV)$accession)
idx2_vsn <- sapply(idx_vsn, head, 3)
small_vsn <- qntV[unlist(idx2_vsn), ]

#Prepare data frame for ggplot
#Convert to data frame including Protein and Feature identifiers
dfrM <- data.frame(
  exprs(smallM),                      
  Protein = as.character(fData(smallM)$accession),
  Feature = featureNames(smallM),
  stringsAsFactors = FALSE
)

dfr_vsn <- data.frame(
  exprs(small_vsn),
  Protein = as.character(fData(small_vsn)$accession),
  Feature = featureNames(small_vsn),
  stringsAsFactors = FALSE
)

# Rename TMT reporter columns for clarity
colnames(dfrM) <- c("126", "127", "128", "129", "130", "131", "Protein", "Feature")
colnames(dfr_vsn) <- c("126", "127", "128", "129", "130", "131", "Protein", "Feature")

# Simplify protein names for plotting readability
dfrM$Protein[dfrM$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
dfrM$Protein[dfrM$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
dfrM$Protein[dfrM$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
dfrM$Protein[dfrM$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
dfrM$Protein[grep("ECA", dfrM$Protein)] <- "Background"

dfr_vsn$Protein[dfr_vsn$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
dfr_vsn$Protein[dfr_vsn$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
dfr_vsn$Protein[dfr_vsn$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
dfr_vsn$Protein[dfr_vsn$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
dfr_vsn$Protein[grep("ECA", dfr_vsn$Protein)] <- "Background"

#######8 -Melt for ggplot ---
# Convert wide to long format for ggplot compatibility
dfrM_melt <- melt(dfrM, id.vars = c("Protein", "Feature"))
dfr_vsn_melt <- melt(dfr_vsn, id.vars = c("Protein", "Feature"))

# Plot using ggplot ---
# Median-normalized data
#par(mfrow = c(1, 2)) 
Median_ggplot <- ggplot(dfrM_melt, aes(x = variable, y = value, colour = Protein)) +
  geom_point() +
  geom_line(aes(group = as.factor(Feature)), alpha = 0.5) +  # Connect points per peptide
  facet_grid(. ~ Protein) +                                   # Facet by protein
  theme(legend.position = "none") +
  labs(x = "Reporters", y = "Median-normalised intensity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    # Rotate x-axis labels
ggsave("figs/Median_ggplot.png", Median_ggplot)



# VSN/second normalized data
VSN_ggplot <- ggplot(dfr_vsn_melt, aes(x = variable, y = value, colour = Protein)) +
  geom_point() +
  geom_line(aes(group = as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) +
  theme(legend.position = "none") +
  labs(x = "Reporters", y = "Median-normalised intensity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figs/VSN_ggplot.png", VSN_ggplot)




