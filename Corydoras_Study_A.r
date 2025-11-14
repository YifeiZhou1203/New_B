library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)
library(here)
library(forcats)

#commenting the old code to read the df
#df_fish <- read_tsv(file = "../jade/data/result.tsv")

# 1st edit (adding new reproducible code to read the df)
# Use a project-root-aware path (via here::here) and check that the data file exists.
# This prevents path errors and ensures the script runs reproducibly on any machine.
file_path <- here("result.tsv")

if (!file.exists(file_path)) {
  stop("Error: 'result.tsv' not found in the project root. Please place the file in the correct location.")
}

df_fish <- read_tsv(file_path)

head(df_fish)

head(df_fish$coord)


##get and separate coord into lat and lon
df_fish$coord <- gsub("\\[|\\]", "", df_fish$coord)


de_fishA <- df_fish %>%
  separate(coord, into = c("lat", "lon"), sep = ",", remove = FALSE)

de_fishA$lat <- as.numeric(trimws(de_fishA$lat))
de_fishA$lon <- as.numeric(trimws(de_fishA$lon))

head(de_fishA[, c("coord", "lat", "lon")])



df_fishC <- de_fishA %>%
  filter(!is.na(bin_uri) & !is.na(lat) & !is.na(lon))

bin_summary <- df_fishC %>%
  group_by(species) %>%
  summarise(n_bins = n_distinct(bin_uri), n_samples = n())
print(bin_summary)


#bin_matrix <- table(df_fishC$bin_uri, df_fishC$species)
#bin_matrix_fix <- (bin_matrix > 0)*1

bin_plot_south <- df_fishC %>%
  filter(lon >=-90 & lon<=-30,
         lat >=-40 & lat<=15)

# commenting the old code due to title spelling issues, not proportional, and more.
# ggplot(bin_plot_south, aes(x=lon,y=lat, color = bin_uri))+
#   borders(xlim=c(-80, -30), ylim=c(-40, 15), fill = "gray90")+
#   geom_point(size=2, alpha=0.8)+
# 
#   labs(title="Distribution of BINs acorss south America", x= "longitude", y="latitude")+
#   theme(legend.position= "right")

# edit #2
# Improved map visualization: cleaner theme, fixed title, proportional axes, clearer legend
ggplot(bin_plot_south, aes(x = lon, y = lat, color = bin_uri)) +
  borders(xlim = c(-80, -30), ylim = c(-40, 15), fill = "gray90") +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_viridis_d(name = "BIN ID") +
  coord_fixed() +
  labs(title = "Distribution of BINs Across South America",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "right")


bin_matrix_south <- table(bin_plot_south$bin_uri, bin_plot_south$species)
bin_matrix_fix_south <- (bin_matrix_south > 0)*1




bin_region <- bin_plot_south %>%
  group_by(country_iso, bin_uri) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  pivot_wider(names_from = bin_uri, values_from = n_records, values_fill = 0)

comm_matrix <-as.matrix(bin_region[,-1])
rownames(comm_matrix) <- bin_region$country_iso
comm_matrix_pa <-(comm_matrix >0)*1
any(is.na(comm_matrix_pa))


set.seed(42)
nmds <- metaMDS(comm_matrix_pa, distance = "bray", k =2)

nmds_points <- as.data.frame(nmds$points)
nmds_points$country_iso <-rownames(nmds_points)
ggplot(nmds_points, aes(x=MDS1, y=MDS2, shape = country_iso))+
  geom_point(size = 4, alpha = 0.8)+
  labs(title = "NMDS OF BIN Composition Across South America",
       x= "NMDS1", y="NMDS2")

# commenting the old code here (It hides the biodiversity pattern and could use improvements to better match the scientific question)
## bargram
# bin_his <- bin_plot_south %>%
#   group_by(country_iso) %>%
#   summarise(unique_bins = n_distinct(bin_uri))
# 
# ggplot(bin_his, aes(x=country_iso, y=unique_bins, pattern = country_iso))+
#   geom_bar (stat="identity")+
#   labs(title="Distribution of Bin Across countries", x="Countries", y="Bin numbers")+
#   theme(legend.position = "none")


# Given the Reorder countries by BIN richness instead of alphabetical order.
# The original plot used the default A–Z ordering, which hid biodiversity patterns.
# Using fct_reorder() arranges countries by unique BIN count, making geographic trends clearer.
bin_his <- bin_plot_south %>%
  group_by(country_iso) %>%
  summarise(unique_bins = n_distinct(bin_uri)) %>%
  mutate(country_iso = fct_reorder(country_iso, unique_bins))

ggplot(bin_his, aes(x = country_iso, y = unique_bins)) +
  geom_col(fill = "steelblue") +
  labs(title = "BIN Richness by Country",
       x = "Country", y = "Number of Unique BINs") +
  theme_minimal()

# Overall: File handling and reproducibility were improved by replacing hard-coded paths
# with here::here() and adding basic checks for the input file.

# Overall: Visualization clarity was enhanced by correcting titles, using a cleaner theme,
# enforcing proportional axes, and reordering countries by BIN richness.

# Overall: Visualization was improved to better support the scientific question.

# Note: The script would benefit from more inline comments explaining each processing step,
# especially during data cleaning and matrix construction, to improve readability for collaborators.


