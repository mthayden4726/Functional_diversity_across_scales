#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:22:47 2023

@author: meha3816
"""

## Load necessary packages ##
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
import parmap
import os
import tqdm 
from progress.bar import Bar
from tqdm.contrib.concurrent import process_map
from multiprocessing import Pool, cpu_count
from window_calcs import * # add scripts folder to python path manager
import glob
from matplotlib.pyplot import subplots, show

Data_Dir = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata'

SITECODE = 'SERC' # NEON site of interest
PRODUCTCODE = 'DP3.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2019-05'

files1 = os.listdir(Data_Dir)
files_TEAK = glob.glob(Data_Dir + "/NEON_D17*")
files_SRER = glob.glob(Data_Dir + "/NEON_D14*")
files_TALL = glob.glob(Data_Dir + "/NEON_D08*")
files_KONZ = glob.glob(Data_Dir + "/NEON_D06*") # never downloaded, too slow
files_SERC = glob.glob(Data_Dir + "/NEON_D02*") # never downloaded, too slow



file_paths = find_neon_files(SITECODE,
                             PRODUCTCODE,
                          YEAR)
files = retrieve_neon_files(file_paths, Data_Dir)

# count proportion of cells in scene with ndvi == 0
list = []
for file in files_SERC:
    neon = ht.HyTools() 
    neon.read_file(file,'neon')
    # show_rgb(neon, r=660,g=550,b=440, correct = [])
    ndvi = neon.ndi()
    zero_els = np.count_nonzero(ndvi==0)
    list.append(file)
    list.append((zero_els/(1000*1000))*100)
    #plt.imshow(ndvi)
    #plt.show()

# plot scenes 
fig, axs = plt.subplots(nrows=5, ncols=4)
for file, ax in zip(files_TALL, axs.ravel()):
    neon = ht.HyTools() 
    neon.read_file(file,'neon')
    # show_rgb(neon, r=660,g=550,b=440, correct = [])
    ndvi = neon.ndi()
    # filter df for ticker and plot on specified axes
    # show image
    shw = ax.imshow(ndvi)
      
    # make bar
    bar = plt.colorbar(shw)
      
    # show plot with labels
    plt.xlabel('Meters')
    plt.ylabel('Meters')
    
plt.show()


# For a single plot
# make plot
#fig, ax = plt.subplots()
  
# show image
#shw = ax.imshow(ndvi)
  
# make bar
#bar = plt.colorbar(shw)
  
# show plot with labels
#plt.xlabel('Meters')
#plt.ylabel('Meters')
#bar.set_label('NDVI')
#plt.show()
