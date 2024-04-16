#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:28:57 2023

@author: meha3816
"""

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
import os, glob
import tqdm 
from progress.bar import Bar
from tqdm.contrib.concurrent import process_map
from multiprocessing import Pool, cpu_count
from S01_Moving_Window_FRIC import * # add scripts folder to python path manager
from S01_Functions_KONZ import * # add scripts folder to python path manager
from osgeo import gdal, osr
from matplotlib.pyplot import subplots, show
import h5py
import sys
import rasterio
import boto3
import re
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError

# Set specifications 
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

nir_band = 90
red_band = 58
ndvi_threshold = 0.25
epsg = 32619

SITECODE = 'ONAQ'
YEAR = '2019-05'

SERVER = 'http://data.neonscience.org/api/v0/'

## Find files for:
# Elevation (DP3.30024.001)

PRODUCTCODE = 'DP3.30024.001'

# Build so that we can loop through reading files in
url = SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/'+YEAR
# Request the url
data_request = requests.get(url)
# Convert the request to Python JSON object
data_json = data_request.json()
print(data_json)
# Create list of file paths of interest for given site, product, year
file_paths = []
for file in data_json['data']['files'][:]:
  if 'DTM.tif' in file['name']:
    file_paths.insert(1, file['url'])
    print(file['url'])
print(file_paths)

# Canopy Height (DP3.30015.001)

# Slope/Aspect (DP3.30025.001) 



# Loop through all CLBJ files
#for i,file in enumerate(file_names):
    #print(file)
    #flight = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D01/2019_HARV_6/L1/Spectrometer/ReflectanceH5/2019082013/NEON_D01_HARV_DP1_' + file +'_reflectance.h5'
    #files = []
    #files.append(flight)
    #try:
    #    retrieve_neon_files(files, Data_Dir)
    #except Exception as e:
    #    continue 
    #img = Data_Dir + "/NEON_D01_HARV_DP1_" + file + '_reflectance.h5'
    #neon = ht.HyTools() 
    #neon.read_file(img,'neon')
    # Store map info for raster
    #refl_md, header_dict = store_metadata(neon, epsg)

    #print("uploading array")
    #upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    #os.remove(local_file_path)
    #os.remove(img)
    #print("flightline complete")
