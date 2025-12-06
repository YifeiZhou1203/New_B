# BINF 6210 – Assignment 4
# Pipeline for MSnbase DDA analysis  
# Author: Yifei Zhou
# Date: 2025-12-05
##############################################################
#install and Load packages


BiocManager::install("msdata")

library(MSnbase)
library(msdata)
library(rgl)
library(plotly)

######1 -Overview----
#MS1 Pre-processing
#2. MS2 Precursor Extraction
#3. MS2 → MS1 Matching (Precursor Validation)
#Compute m/z tolerance. 
#Compare MS2 precursor m/z to MS1 peaks within m/z window.
#Apply retention time tolerance 
#4. Spectrum Visualization
#Plot the first MS2 spectrum to inspect fragmentation quality.
#Generate Total Ion Chromatogram (TIC) from MS1 scans:
#5. 3D MS1 Visualization
#Select first 50 MS1 spectra.
#Define a common m/z grid and interpolate each spectrum.
#Generate a 3D waterfall/surface plot (m/z × RT × intensity).

dir.create("figs",    showWarnings = FALSE, recursive = TRUE)


######1 -load data----
#Load example DDA mzML file: "TMT_Erwinia_DDA.mzML"
f <- proteomics(full.names = TRUE)
data_file <- f[4]

raw <- readMSData(data_file, mode="onDisk", centroided = FALSE)


# Check how many scans are MS1 vs MS2
table(msLevel(raw))

######2 -Pre-processing----
#Separate MS1 and MS2 scans using filterMsLevel(): extracts scans MS levels at 1 or 2
ms1 <- filterMsLevel(raw, 1)
ms2 <- filterMsLevel(raw, 2)

cat("MS1 scans:", length(ms1), "\n")
cat("MS2 scans:", length(ms2), "\n")


# chromatograms, precursor matching require MS1; while MS2 is for precursor analysis



# Smooth MS1 spectra + peak picking

ms1_smooth <- smooth(ms1, method="SavitzkyGolay", halfWindowSize=4)
ms1_peaks <- pickPeaks(ms1_smooth, SNR=5, method="MAD")

# Pick peaks using median absolute deviation estimator; SNR = 5 for default; increases sensitivity
# That creates the MS1 peak lists for matching MS2 precursor m/z



#Inspect the first MS1 spectrum:
#Select the first MS1 spectrum
sp <- ms1[[1]]

# Plot raw intensities as profile data
MS1 <- plot(mz(sp), intensity(sp), type="l")
title("First MS1 Spectrum")

ggsave("figs/MS1.png", MS1)

######3 -Prepare MS1 data vectors: extract all MS1 spectra-------
sp_list <- spectra(ms1_peaks)

# Create flat vectors of all m/z and their retention times
#matching between MS2 precursor m/z and MS1 peaks across time.
all_mz <- unlist(lapply(sp_list, mz))
all_rt <- unlist(lapply(sp_list, rtime))




############4 -MS2 → MS1 precursor matching (DDA linking)-----
ppm <- 20
rt_window <- 15   # RT tolerance

# Extract precursor m/z and MS2 retention times
prec <- precursorMz(ms2)
ms2_rt <- rtime(ms2)

# Remove MS2 scans that have missing precursor m/z
valid <- !is.na(prec)

# Loop through each valid MS2 scan to compare its precursor m/z against MS1 peaks within ppm + RT window selected earlier; Determine whether the precursor selected in MS2 is visible in the MS1 peak map

matches <- lapply(which(valid), function(i) {
  
  pmz <- prec[i]          # precursor m/z
  rt  <- ms2_rt[i]        # MS2 retention time
  tol <- pmz * ppm / 1e6  # ppm → absolute m/z tolerance
  
  # Find all MS1 peaks that fall within tolerance windows
  idx <- which(
    abs(all_mz - pmz) < tol &
      abs(all_rt - rt) < rt_window
  )
  
  # Return summary for that MS2 scan
  data.frame(
    MS2_scan = acquisitionNum(ms2)[i],
    Precursor_mz = pmz,
    MS2_RT = rt,
    Matched_MS1_peaks = length(idx)
  )
})

matches <- do.call(rbind, matches)
head(matches)


#Extract MS1 chromatogram for a precursor
#presence, intensity, and retention behavior for selected MS1 ion

i <- 1  # first MS2 scan
pmz <- precursorMz(ms2)[i]
tol <- pmz * ppm / 1e6


chr <- chromatogram(ms1_peaks, mz = c(pmz - tol, pmz + tol))

if (length(chr) > 0 && !isEmpty(chr[[1]])) {
  plot(chr[[1]])
  title(paste("Chromatogram for precursor m/z", round(pmz, 4)))
} else {
  cat("No MS1 peaks found in this m/z window\n")
}



#######5 -Plot first MS2 spectrum-------
MS2 <- plot(ms2[[1]])
title(paste("MS2 spectrum scan", acquisitionNum(ms2)[1]))

ggsave("figs/MS2.png", MS2)

# Total Ion Chromatogram: Detect instrument drift, spray instability, or column issues

tic_values <- tic(ms1)        # total intensity per MS1 scan
rt_ms1 <- rtime(ms1)          # retention time for each MS1 scan

plot(rt_ms1, tic_values, type = "l",
     xlab = "Retention Time (sec)", ylab = "TIC",
     main = "Total Ion Chromatogram (TIC)")




####
#######6 -3D Waterfall MS1 Plot------
####
# Convert first 50 MS1 spectra for a clean waterfall
n <- 50
sp_list <- spectra(ms1)[1:n]

mz_grid <- seq(400, 1600, by = 1)    # define m/z grid

# interpolate each spectrum onto common m/z grid
int_matrix <- sapply(sp_list, function(sp) {
  approx(mz(sp), intensity(sp), xout = mz_grid, rule = 2)$y
})



plot_ly() %>%
  add_surface(
    x = mz_grid,
    y = seq_len(n),           # spectrum index / RT
    z = int_matrix,
    colorscale = "Viridis"
  ) %>%
  layout(
    title = "3D MS1 Waterfall Plot",
    scene = list(
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Spectrum index (RT order)"),
      zaxis = list(title = "Intensity")
    )
  )

