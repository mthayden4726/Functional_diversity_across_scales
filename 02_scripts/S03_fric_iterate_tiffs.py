#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for iterating through each file per site to calculate
spectral diversity metrics.
- Alpha diversity
- Functional richness

Uses functions written in 'S01_specdiv_functions.py'


Author: M. Hayden
Last Updated: 6/8/2023
"""
#################################  SET-UP  ###################################

## Import packages ##
import hytools as ht
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import requests
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import kneed
from kneed import KneeLocator
from scipy.spatial import ConvexHull
import subprocess
from urllib.request import urlretrieve
import multiprocessing as mp
import os, glob
import csv
import rasterio
from osgeo import gdal
import rioxarray as rxr
import xarray as xr
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
import copy
import re
# Import functions defined in S01_specdiv_functions.py
from S01_specdiv_functions import * # add scripts folder to python path manager

## Set global parameters ##
window_sizes = [10, 50, 100, 150, 200, 300, 400]   # list of window sizes to test
ndvi_threshold = 0.4 # ndvi threshold for radiometric filtering
# shade_threshold = 500 # should be low threshold across NIR region (15%)
# cloud_threshold = 1500 # should be high threshold in blue band
bad_bands = [[300,400],[1300,1450],[1780,2000],[2450,2600]] # bands to be masked out
sample_size = 0.1 # proportion of pixels to subsample for fitting PCA
comps = 4 # default component numbers for PCA 
nclusters = 15 # default component numbers for K-means clustering
nbCPU = 4 # computational parameters
MaxRAM = 0.6 # computational parameters

## Set parameters for desired imagery & local output ##
# Define path for master output directory where files produced during the process are saved
Data_Dir = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/02_processed/TEAK_flightlines'
Output_Dir = '01_data/02_processed'

#files = os.listdir(Data_Dir)
#file = glob.glob(Data_Dir + "/TEAK*")
#file = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/02_processed/TEAK_flightlines/TEAK_20190615_171836.tif'
#naip_csf = rxr.open_rasterio(file, masked=True)
#naip_csf


all_bands = glob.glob(os.path.join(Data_Dir,"20190615_171836_output_*.tif"))
all_bands.sort()
post_all_bands = []
for i, aband in enumerate(all_bands):
    cleaned = rxr.open_rasterio(aband,
                                     masked=True).squeeze()
    # This line below is only needed if you wish to stack and plot your data
    cleaned["band"] = list(map(int, re.findall(r'\d+', str(all_bands[i]))))[6]
    post_all_bands.append(cleaned)

ndvi = es.normalized_diff(post_all_bands[400], post_all_bands[360])
ep.plot_bands(ndvi, cmap="RdYlGn", cols=1, vmin=-1, vmax=1)

# Calculate NDVI
#ndvi = es.normalized_diff(post_stack[400], post_stack[360])
# Plot NDVI

# Create classes and apply to NDVI results
ndvi_class_bins = [-np.inf, 0.4, np.inf]
ndvi_class = np.digitize(ndvi, ndvi_class_bins)

# Apply the nodata mask to the newly classified NDVI data
ndvi_class = np.ma.masked_where(
    np.ma.getmask(ndvi), ndvi_class
)
np.unique(ndvi_class)

# Plot classes
# Define color map
nbr_colors = ["gray", "yellowgreen"]
nbr_cmap = clr.ListedColormap(nbr_colors)

# Define class names
ndvi_cat_names = [
    "No Vegetation",
    "Vegetation",
]

# Get list of classes
classes = np.unique(ndvi_class)
classes = classes.tolist()
# The mask returns a value of none in the classes. remove that
classes = classes[0:2]

# Plot your data
fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(ndvi_class, cmap=nbr_cmap)

ep.draw_legend(im_ax=im, classes=classes, titles=ndvi_cat_names)
ax.set_title(
    "Normalized Difference Vegetation Index (NDVI) Classes",
    fontsize=14,
)
ax.set_axis_off()

# Auto adjust subplot to fit figure size
plt.tight_layout()

# Mask based on classes
masked_bands = []
for i, aband in enumerate(post_all_bands):
    veg = xr.where(
     ndvi_class == 2, aband, "nan"
    )
    masked_bands.append(veg)
    
# Adapting functions for stack of bands 
# needs a numpy array

post_stack = xr.concat(masked_bands, dim="band")
post_stack

 
#####################  CALCULATE FRIC FOR SELECTED SITE  ##########################

## Create dicts for storing scale-FRic and scale-CV output for each file ##
scale_fric = {} 
#scale_cv = {}


# Loop through all files for the site
index = 0
for index, i in enumerate(files):
    
    # Read image
    neon_image = file
    #neon_image = files[index]
    neon = ht.HyTools()
    neon.read_file(neon_image)
    
    # Mask out bad bands and non-vegetation (eventually it's own function)
    neon.create_bad_bands(bad_bands)
    ndvi = neon.ndi()
    neon.mask['sample'] = (neon.mask['no_data']) & (ndvi > ndvi_threshold)
    # Calculate alpha diversity based on CV
    #adiv = calc_cv(neon, window_sizes)
    #scale_cv[i] = adiv
    
    # Scale, center and PCA transform
    X = subsample(neon,sample_size,bad_bands)
    x_mean, x_std, pca = scale_transform(X, comps)
    
    # Calculate functional richness based on PCA components
    volumes = calc_fun_rich_no_iter(neon, window_sizes,
                            x_mean, x_std, pca, comps)
    scale_fric[index] = volumes
    index += 1
    
# Save to .txt
# open file for writing
f = open("01_data/02_processed/TEAK_FRic.txt","w")
# write file
f.write(str(scale_fric))
# close file
f.close()

## Export file ##
# Save as numpy file
#np.save('CPER_FRic_scales.npy', scale_fric) 
# Save to .csv
# w = csv.writer(open("output.csv", "w"))
#dict = scale_fric
# loop over dictionary keys and values
#for key, val in dict.items():
#    # write every key and value to file
#    w.writerow([key, val])
#w.close()
