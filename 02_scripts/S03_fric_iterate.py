#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:19:30 2023

@author: meha3816
"""

#################################  SET-UP  #################################

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
import os
# Import functions defined in S01_specdiv_functions.py
from S01_specdiv_functions import * # add scripts folder to python path manager

## Set global parameters ##
window_sizes = [10, 15, 20, 30, 40, 50, 75, 100, 150]  # list of window sizes to test
ndvi_threshold = 0.4 # ndvi threshold for radiometric filtering
# shade_threshold = 500 # should be low threshold across NIR region (15%)
# cloud_threshold = 1500 # should be high threshold in blue band
bad_bands = [[300,400],[1300,1450],[1780,2000],[2450,2600]] # bands to be masked out
comps = 4 # default component numbers for PCA 
nclusters = 15 # default component numbers for K-means clustering
nbCPU = 4 # computational parameters
MaxRAM = 0.2 # computational parameters

## Set parameters for desired imagery & local output ##
# Define path for master output directory where files produced during the process are saved
Output_Dir = 'Users/meha3816/Desktop/BioSCape_across_scales/01_data/02_processed'
Temp_Dir = 'Users/meha3816/Desktop/BioSCape_across_scales/01_data/temp_neon' # for files downloaded from NEON API

# Set parameters for data import
SITECODE = 'CPER' # NEON site of interest
PRODUCTCODE = 'DP3.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2021-06' # Timeframe of desired imagery

## Create list of data files to process in remainder of script ##
file_paths = find_neon_files(SITECODE,
                             PRODUCTCODE,
                             YEAR)
# Download files to temp OS location - WARNING: LARGE STORAGE REQUIREMENT
files = retrieve_neon_files(file_paths)

#####################  CALCULATE FRIC FOR THE SITE  #####################
scale_fric = {} # dict for storing output
for i in range(len(files)):
    # Read image
    neon_image = files[i]
    neon = ht.HyTools()
    neon.read_file(neon_image,'neon')
    # Mask out bad bands and non-vegetation (eventually it's own function)
    neon.create_bad_bands(bad_bands)
    ndvi = neon.ndi()
    neon.mask['sample'] = (neon.mask['no_data']) & (ndvi > ndvi_threshold)
    # Scale, center and PCA transform
    X = subsample(neon, sample_size = 0.1)
    x_mean, x_std = scale_transform(X)
    # Calculate functional richness based on PCA components
    scale_fric[i] = calc_fun_rich(neon, window_sizes, x_mean)
    
 
