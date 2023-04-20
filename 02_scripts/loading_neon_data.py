#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:02:37 2023

@author: meha3816
"""

import numpy as np
import hytools as ht
#from hytools_lite.io.envi import WriteENVI
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from skimage.util import view_as_blocks
from numba import jit
import matplotlib.pyplot as plt
#from maap.maap import maap
import IPython
import os
import requests
import json
import itertools
import kneed
from kneed import KneeLocator
import requests
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
base = importr('base')
utils = importr('utils')
stats = importr('stats')
# for hytools need msgpack, google-cloud, google-cloud-vision, mne
# brainpipe (not available for pip), kneed, numba, skimage

# Create file path
# 3. Create path to desired NEON files

# Server's URL will remain the same
SERVER = 'http://data.neonscience.org/api/v0/'
# Options to specify sitecode, productcode, and timeframe
SITECODE = 'SRER'
PRODUCTCODE = 'DP1.30006.001'
YEAR = '2019-09'
# Build so that we can loop through reading files in
url = SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/'+YEAR
# print(url)

# Request the url
data_request = requests.get(url)
print(data_request.status_code) #200 means success

# Convert the request to Python JSON object
data_json = data_request.json()

file_paths = []
for file in data_json['data']['files'][:20]:
  if 'reflectance.h5' in file['name']:
        file_paths.insert(1, file['url'])
        print(file['url'])

url = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090414/NEON_D14_SRER_DP1_20190904_165832_reflectance.h5'    
        
img = Image.open(url)
print(file_paths)

