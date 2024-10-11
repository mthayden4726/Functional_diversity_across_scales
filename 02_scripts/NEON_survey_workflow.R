## This script uses Adam Mahood's package for summarizing plant species richness 
## for NEON subplots

# Updated Oct. 11 2024

devtools::install_github("admahood/neonPlantEcology")
library(neonPlantEcology)

# Specify sites of interest
sites <- c("SRER", 
           "TEAK",
           "SERC",
           "KONZ",
           "BART",
           "OSBS",
           "ONAQ",
           "HEAL",
           "PUUM",
           "UNDE",
           "TALL",
           "WOOD",
           "CLBJ",
           "YELL",
           "NIWO",
           "WREF",
           "TOOL")
year <- 2019

# Specify sites to download
sites <- npe_download(sites = sites,
                      product = "plant_diversity")

# Convert to matrix and filter to year
species_occurrence_matrix <- npe_community_matrix(sites, binary=TRUE)

# Get coordinates
data("plot_centroids")
centroids <- npe_plot_centroids(species_occurrence_matrix, spatial_only = FALSE) %>%
  filter(subtype == "basePlot") %>%
  select(plotID, plotSize, latitude, longitude, 
         elevation, minElev, maxElev, slope, aspect,
         nlcdClass) %>%
  unique()

# Get summary for site and select variables of interest
species_summary <- npe_summary(sites) %>%
  filter(eventID == year) %>%
  select(site, plotID, shannon_total, evenness_total, nspp_total, invaded)

# Bind files
plot_data <- merge(species_summary, centroids, by = "plotID")

# Write to file
write.csv(plot_data, "NEON_species_summaries.csv")
