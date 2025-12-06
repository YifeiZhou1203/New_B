# BINF 6210 – Assignment 4
# MSnbase 
# Author: Yifei Zhou
# Date: 2025-12-05

library("rpx")          
library("Biobase")      
library("MSnbase")      
library("RColorBrewer")    
library("reshape2")      
library("ggplot2")      

############################################################


dir.create("figs",    showWarnings = FALSE, recursive = TRUE)



#######1 -load peptide level data----
#Download PXDataset PXD000001 TMT benchmark dataset(same procedures as previous analysis)
px1 <- PXDataset("PXD000001")
px1


pxfiles(px1)

#Load mzTab peptide quantification file



######2 -peptide-level quantification---- 
mztab <- pxget(px1, "F063721.dat-mztab.txt")
mztab

#Read peptide (PEP) section of mzTab as an MSnSet
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")

#Assign TMT6 reporter channel names
sampleNames(qnt) <- reporterNames(TMT6)

#Inspect 
head(exprs(qnt))

#Remove missing data
qnt <- filterNA(qnt)

#processing metadata
processingData(qnt)


############################################################
#######3 -load MS2 data----
#Load raw mzXML MS/MS data for direct reporter quantification (MS2 spectra)

mzxml <- pxget(px1, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML")

#Read MS data
rawms <- readMSData(mzxml, centroided = TRUE, verbose = FALSE)
rawms


length(rawms)

#inspect
table(msLevel(rawms))


#####4 -Quantify reporter ions using MSnbase------
#extract TMT reporter intensities from MS2 spectra

qntms <- quantify(rawms, reporters = TMT7, method = "max")
qntms



######5 -Prepare data for incomplete dissociation analysis-----
# Signal = total intensity of TMT6 reporter channels; Incomplete = channel for incomplete dissociation (reporter 7)
d <- data.frame(
  Signal = rowSums(exprs(qntms)[, 1:6]),
  Incomplete = exprs(qntms)[, 7]
)

# Log-transform intensities for linearity
d <- log(d)

# Default plotting aesthetics
cls <- rep("#00000050", nrow(qnt))
pch <- rep(1, nrow(qnt))



#6- QC Plot-----
# examine relationship between TMT signal and incomplete dissociation/plot 

Reference_QC <- plot(Signal ~ Incomplete, data = d,
     xlab = expression(Incomplete~dissociation),
     ylab = expression(Sum~of~reporters~intensities),
     pch = 19,
     col = "#4582B380")
    
grid()

abline(0, 1, lty = "dotted")
title(main = "Reference incomplete dissociation plot", 
      cex.main = 1.5, font.main = 2, col.main = "black") 

ggsave("figs/Reference_QC.png", Reference_QC)

# Regression line to evaluate linear relationship
abline(lm(Signal ~ Incomplete, data = d), col = "darkblue")


############################################################
#6- Annotation-----
# Annotate the proteins in this dataset: for evaluating quantification accuracy


# set background transparent 
d$color <- "#4582B300"


# BSA, ENO, CYT, and PHO are spike-in controls used for QC (same as the previous Mat-plot colors)
d$color[grep("P02769", fData(qnt)$accession)] <- "gold4"        # BSA
d$color[grep("P00924", fData(qnt)$accession)] <- "dodgerblue"   # ENO
d$color[grep("P62894", fData(qnt)$accession)] <- "springgreen4" # CYT
d$color[grep("P00489", fData(qnt)$accession)] <- "darkorchid2"  # PHO



#Fit linear model to quantify influence of incomplete dissociation
fit <- lm(Signal ~ Incomplete, data = d)
r2 <- summary(fit)$r.squared


#7- Final QC plot with protein highlighted-----
#QC Plot with Standards + Regression + R²
Final_QC <- plot(Signal ~ Incomplete, data = d,
     pch = 19,
     col = d$color,
     xlab = expression(Incomplete~dissociation),
     ylab = expression(Sum~of~reporters~intensities))
     main = "QC Plot: Highlighting Standard Proteins"
grid()

#add reference line 
abline(0, 1, lty = "dotted")

#add regression line 
abline(fit, col = "darkblue", lwd = 2)

#print R² 
text(
  x = max(d$Incomplete, na.rm=TRUE)*0.6,
  y = max(d$Signal, na.rm=TRUE)*0.9,
  labels = paste0("R² = ", round(r2, 3)),
  col = "darkblue",
  cex = 1.2
)

#add legend for standard proteins
legend("topleft", legend=c("BSA","ENO","CYT","PHO","Other"),
       col=c("gold4","dodgerblue","springgreen4","darkorchid2","#4582B3"),
       pch=19,
       bty="n")

title(main = "QC Plot: Highlighting Standard Proteins", 
      cex.main = 1.5, font.main = 2, col.main = "black")  

ggsave("figs/Final_QC <-.png", Final_QC)
