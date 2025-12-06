# BINF 6210 – Assignment 4
# Pipeline for MALDIquant and MyclassifR package
# Author: Yifei Zhou
# Date: 2025-12-05


# Load required packages
library(MALDIquant)        
library(MALDIquantForeign) 
library(MSclassifR)       


######1 -Overview----
#Load Libraries & Dataset
#Preprocessing of Spectra
#For each spectrum: Variance Stabilization; Smoothing;Baseline Correction
#Peak Detection
#Alignment & Matrix Generation
#Visualization of Preprocessing Effects
#Exploratory Data Analysis (PCA)
#Visualize PC1 vs PC2 to observe clustering patterns.
#Assess spectral similarity and identify outliers.




dir.create("figs",    showWarnings = FALSE, recursive = TRUE)


######2 -load data----
data("CitrobacterRKIspectra", "CitrobacterRKImetadata", package="MSclassifR")
spectra <- CitrobacterRKIspectra

#view number of spectra
length(spectra)


######3 -pre-processing----
#Preprocessing specta:reducing noise, stabilizing variance, smoothing, and removing baseline
#Variance stabilization: sqrt reduces heteroscedasticity
#Smoothing: Moving Average 
#Baseline correction: SNIP method 
spectra_proc <- lapply(spectra[1:5], function(s) {
  
  
  s_var <- transformIntensity(s, method="sqrt")
  
  
  s_smooth <- smoothIntensity(s_var, method="MovingAverage", halfWindowSize=2)
  
 
  removeBaseline(s_smooth, method="SNIP")
})

#Peak detection 

peaks_list <- lapply(spectra_proc, detectPeaks, SNR=2, halfWindowSize=2)

# Create peak intensity matrix: rows=peaks, columns=spectra
peakMatrix <- intensityMatrix(peaks_list, spectra_proc)



######3 -Visualize raw vs processed spectrum----
par(mfrow=c(2,2))  

# Raw spectrum
plot(spectra[[1]], main="1: Raw Spectrum", xlab="m/z", ylab="Intensity")

# Variance stabilized
plot(transformIntensity(spectra[[1]], method="sqrt"),
     main="2: Variance Stabilization", xlab="m/z", ylab="Intensity")

# Smoothed spectrum
plot(smoothIntensity(transformIntensity(spectra[[1]], method="sqrt"),
                     method="MovingAverage", halfWindowSize=2),
     main="3: Smoothing", xlab="m/z", ylab="Intensity")

# Baseline corrected
plot(spectra_proc[[1]], main="4: Baseline Corrected", xlab="m/z", ylab="Intensity")

# Overlay detected peaks
points(peaks_list[[1]], col="red", pch=19)
legend("topright", legend="Detected Peaks", col="red", pch=19)



######4 -Overlay preprocessed spectra----
colors <- rainbow(3)

Overlay <- plot(spectra_proc[[1]], type="l", col=colors[1],
     main="Overlay of First 3 Preprocessed Spectra",
     xlab="m/z", ylab="Intensity")

for(i in 2:3){
  lines(spectra_proc[[i]], col=colors[i])
}

legend("topright", legend=paste0("Spectrum ", 1:3), col=colors, lty=1)

ggsave("figs/Overlay.png", Overlay)

# Principal Component Analysis: PCA Reduces dimensionality and explore spectral differences
pca_res <- prcomp(t(peakMatrix), scale.=TRUE)  # spectra as rows

# Summary of variance explained
summary(pca_res)



######5 -PCA plot----
# Visualize clustering of first 3 spectra in PCA space
shapes <- c(19, 17, 15)

PCA_plot <- plot(pca_res$x[1:3,1], pca_res$x[1:3,2],
     col=colors, pch=shapes,
     xlab="PC1", ylab="PC2",
     main="PCA of First 3 Spectra")

# Add labels to points
text(pca_res$x[1:5,1], pca_res$x[1:3,2], labels=paste0("Spec", 1:3), pos=3, col=colors)

# Add legend
legend("topright",
       legend=paste0("Spectrum ", 1:3),
       col=colors,
       pch=shapes,
       title="Spectra")

ggsave("figs/PCA_plot.png", PCA_plot)





 