# BINF 6210 – Assignment 4
# Pipeline for missed-cleavage analysis using OrgMassSpecR 
# Author: Yifei Zhou
# Date: 2025-12-05




BiocManager::install("BRAIN")
BiocManager::install("Rdisop")
BiocManager::install("OrgMassSpecR")
BiocManager::install("cleaver")


library(BRAIN)
library(Rdisop)
library(OrgMassSpecR)
library(OrgMassSpecR)
library(Biobase)



###############################################
######1 -Overview----
##MISSED CLEAVAGES
## Test MC = 0, 1, 2 for ENO1 digestion

#1. Load required packages and Retrieve experimental proteomics data
#2. Load target protein sequence from FASTA
#3. Digest protein with different missed-cleavage settings (0, 1, 2)
#4. Extract experimental peptides and compare theoretical vs experimental peptides:
    #Count matches between filtered theoretical peptides and experimental peptides
    #calculate sequence coverage: matched residues / total protein residues
#5. Plot results
 #Barplot: number of matched peptides vs missed cleavages
 #Histogram: peptide length distributions
 #Density plot: visualize length distribution differences across missed-cleavage settings
 #Table of missed cleavages, theoretical peptides, filtered peptides, experimental matches

###############################################


#######2 -Load a ProteomeXchange dataset PXD000001 using PXDataset---- 
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



########3 --Load ENO1 sequence (download if needed)----
if (!file.exists("P00924.fasta")) {
  download.file(
    "http://www.uniprot.org/uniprot/P00924.fasta",
    destfile = "P00924.fasta"
  )
}

eno_seq <- paste(readLines("P00924.fasta")[-1], collapse = "")

########4 --Digest ENO with different missed-cleavage settings----
eno_mc0 <- Digest(eno_seq, missed = 0)
eno_mc1 <- Digest(eno_seq, missed = 1)
eno_mc2 <- Digest(eno_seq, missed = 2)

#Extract experimental peptides from qnt
exppep <- as.character(fData(qnt[grep("ENO", fData(qnt)[, 2]), ])[, 1])

#Minimum length of experimental peptides
minlength <- min(nchar(exppep))


########5 --processing----
#count length-filtered peptides
filter_len <- function(df, thr) df[nchar(df$peptide) >= thr, ]

#count matches
count_matches <- function(theoretical, experimental) {
  sum(theoretical$peptide %in% experimental)
}

## Apply length filter
mc0_f <- filter_len(eno_mc0, minlength)
mc1_f <- filter_len(eno_mc1, minlength)
mc2_f <- filter_len(eno_mc2, minlength)

#######6 --Build summary table----
results <- data.frame(
  MissedCleavages = c(0, 1, 2),
  
  Total_Theoretical = c(
    nrow(eno_mc0),
    nrow(eno_mc1),
    nrow(eno_mc2)
  ),
  
  LengthFiltered = c(
    nrow(mc0_f),
    nrow(mc1_f),
    nrow(mc2_f)
  ),
  
  ExperimentalMatches = c(
    count_matches(mc0_f, exppep),
    count_matches(mc1_f, exppep),
    count_matches(mc2_f, exppep)
  )
)

##Print results
print(results)

########7 -visualize ----
# barplot of matched peptides
Peptide_bar_plot <-barplot(
  results$ExperimentalMatches,
  names.arg = results$MissedCleavages,
  xlab = "Missed Cleavages",
  ylab = "# Experimental Peptides Matched",
  main = "Effect of Missed-Cleavages on ENO1 Matches"
)

dir.create("figs",    showWarnings = FALSE, recursive = TRUE)
ggsave("figs/Peptide_bar_plot.png", Peptide_bar_plot)



#####8 - Build a combined data frame of peptide lengths----

len_df <- rbind(
  data.frame(
    length = nchar(eno_mc0$peptide),
    MC = "MC = 0"
  ),
  data.frame(
    length = nchar(eno_mc1$peptide),
    MC = "MC = 1"
  ),
  data.frame(
    length = nchar(eno_mc2$peptide),
    MC = "MC = 2"
  )
)


## Histogram of peptide lengths

Peptide_his <- ggplot(len_df, aes(x = length, fill = MC)) +
  geom_histogram(alpha = 0.8, bins = 25, position = "dodge") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peptide Length Distribution Across Missed-Cleavage Settings",
    x = "Peptide Length (amino acids)",
    y = "Count"
  )
ggsave("figs/Peptide_his.png", Peptide_his)


## Density plot 
Peptide_density <-ggplot(len_df, aes(x = length, color = MC, linetype = MC)) +
  geom_density(size = 1.1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Peptide Length Density Plot",
    x = "Peptide Length (amino acids)",
    y = "Density"
  )
ggsave("figs/Peptide_density.png", Peptide_density)
