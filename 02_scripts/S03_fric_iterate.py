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
import os
import csv
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
Data_Dir = 'Desktop/BioSCape_across_scales/01_data/01_rawdata'
Output_Dir = '/BioSCape_across_scales/01_data/02_processed'

## Set parameters for data import ##
SITECODE = 'TEAK' # NEON site of interest
PRODUCTCODE = 'DP3.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2019-06' # Timeframe of desired imagery

# If files are already downloaded:
# files = os.listdir(Data_Dir)

## Otherwise...
## Create list of data files to process in remainder of script ##
file_paths = find_neon_files(SITECODE,
                             PRODUCTCODE,
                            YEAR)
# Download files to data directory on OS - WARNING: LARGE STORAGE REQUIREMENT
files = retrieve_neon_files(file_paths, Data_Dir)
 
#####################  CALCULATE FRIC FOR SELECTED SITE  ##########################

## Create dicts for storing scale-FRic and scale-CV output for each file ##
scale_fric = {} 
#scale_cv = {}

# Loop through all files for the site
index = 0
for index, i in enumerate(files):
    
    # Read image
    neon_image = Data_Dir + '/' + files[index]
    #neon_image = files[index]
    neon = ht.HyTools()
    neon.read_file(neon_image,'neon')
    
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
f = open("TEAK_FRic.txt","w")
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
