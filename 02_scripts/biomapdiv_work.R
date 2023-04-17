# BiodivmapR script for assessing scale
# Meghan Hayden - 3/23/23

# From git
install.packages("remotes")
remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR')
install.packages("neonUtilities")
remotes::install_github('earthlab/neonhs')

# clean environment
rm(list=ls(all=TRUE));gc()
# load biodivMapR and useful libraries
library(biodivMapR)
library(labdsv)
library(tools)
library(ggplot2)
library(gridExtra)
library(neonUtilities)
library(jsonlite)
library(httr)
library(raster)
library(neonhs)
# library(utils)
# library(stars)
# library(zip)
# library(rgdal)

# ===============================================================================
# create a directory where to store data
Datadir <- '01_DATA'
dir.create(path = Datadir,recursive = T,showWarnings = F)
# name your binary raster with the same name as the online file
NameRaster <- 'NEON_D14_SRER_DP3_520000_3528000_reflectance.h5'
easting = 520000
northing = 3528000
byTileAOP(dpID = "DP3.30006.001", site = "SRER", 
                  year = "2021", 
                  easting = easting,
                  northing = northing, 
                  check.size = T,
                  savepath = Datadir)
Datadir = "~/DP3.30006.001/neon-aop-products/2021/FullSite/D14/2021_SRER_4/L3/Spectrometer/Reflectance/"
Im_path = "DP3.30006.001/neon-aop-products/2021/FullSite/D14/2021_SRER_4/L3/Spectrometer/Reflectance/NEON_D14_SRER_DP3_520000_3528000_reflectance.h5"
image <- raster("DP3.30006.001/neon-aop-products/2021/FullSite/D14/2021_SRER_4/L3/Spectrometer/Reflectance/NEON_D14_SRER_DP3_520000_3528000_reflectance.h5")
hs_paths <- list.files(path='.', pattern = 'reflectance.h5', 
                       recursive = TRUE, full.names = TRUE)

file <- hs_read(Im_path, bands = c(1,50,100,400))
header <- raster::hdr(file, format = "ENVI", filename = NameRaster)
file
plot(file)
################################################################################
##                      Set parameters for biodivMapR                         ##
## https://jbferet.github.io/biodivMapR/articles/biodivMapR_2.html            ##
################################################################################
# Define path for image file to be processed
Input_Image_File <- file
# Define path for corresponding mask file
# Set to FALSE if no mask available
Input_Mask_File <- FALSE
# Define path for master output directory where files produced during the process are saved
Output_Dir <- '02_RESULTS'
dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# Define levels for radiometric filtering
NDVI_Thresh <- 0.8
Blue_Thresh <- 500
NIR_Thresh <- 1500
# Apply normalization with continuum removal?
Continuum_Removal <- TRUE
# Type of dimensionality reduction
TypePCA <- 'SPCA'
# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- FALSE
# window size forcomputation of spectral diversity
window_size <- 10
# computational parameters
nbCPU <- 4
MaxRAM <- 0.2
# number of clusters (spectral species)
nbclusters <- 50
file
BandName <- file@data@names
SpectralBands <- c(382, 627, 878, 2381)
WLunits <- 'Nanometers'
create_hdr(ImPath = Im_path, Sensor = 'MyOwnSensor', 
           SpectralBands = SpectralBands,BandName = BandName, WLunits = WLunits)

print("PERFORM RADIOMETRIC FILTERING")
Input_Mask_File <- perform_radiometric_filtering(Image_Path = Im_path, Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                 NIR_Thresh = NIR_Thresh)

print("PERFORM DIMENSIONALITY REDUCTION")
PCA_Output <- perform_PCA(Input_Image_File = Im_path,
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM,
                          Continuum_Removal = Continuum_Removal)
