library(tidyverse)
library(ggplot2)
library(viridis)
library(vegan)

df_fish <- read_tsv(file = "../jade/data/result.tsv")
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

ggplot(bin_plot_south, aes(x=lon,y=lat, color = bin_uri))+
  borders(xlim=c(-80, -30), ylim=c(-40, 15), fill = "gray90")+
  geom_point(size=2, alpha=0.8)+
  
  labs(title="Distribution of BINs acorss south America", x= "longitude", y="latitude")+
  theme(legend.position= "right")

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



## bargram
bin_his <- bin_plot_south %>%
  group_by(country_iso) %>%
  summarise(unique_bins = n_distinct(bin_uri))

ggplot(bin_his, aes(x=country_iso, y=unique_bins, pattern = country_iso))+
  geom_bar (stat="identity")+
  labs(title="Distribution of Bin Across countries", x="Countries", y="Bin numbers")+
  theme(legend.position = "none")



