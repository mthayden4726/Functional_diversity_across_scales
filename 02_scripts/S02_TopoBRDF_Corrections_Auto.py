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
#from osgeo import gdal, osr
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

# Select site
site = "HARV"

# Set parameters for NDVI threshold
nir_band = 90
red_band = 58
ndvi_threshold = 0.4

# Find correction coefficients (define search terms) - listing first dates for each for now
if site == "HARV":
  search_criteria = 'NEON_D01_HARV_DP1_20190811'
  date_id = '2019081112'
  site_id = '2019_HARV_6'
  domain = 'D01'
elif site == "OSBS":
  search_criteria = 'NEON_D03_OSBS_DP1_20190415'
  date_id = '2019041512'
  site_id = '2019_OSBS_5'
  domain = 'D03'
elif site == "UNDE":
  search_criteria = 'NEON_D05_UNDE_DP1_20190606'
  date_id = '2019060617'
  site_id = '2019_UNDE_3'	
  domain = 'D05'
elif site == "TALL":
  search_criteria = 'NEON_D08_TALL_DP1_20190427'
  date_id = '2019042713'
  site_id = '2019_TALL_5'
  domain = 'D08'
elif site == "WOOD":
  search_criteria = 'NEON_D09_WOOD_DP1_20190722'
  date_id = '2019072214'
  site_id = '2019_WOOD_3'	
  domain = 'D09'
elif site == "CLBJ":
  search_criteria = 'NEON_D11_CLBJ_DP1_20190419'
  date_id = '2019041914'
  site_id = '2019_CLBJ_4'
  domain = 'D11'
elif site == "YELL":
  search_criteria = 'NEON_D12_YELL_DP1_20190720'
  date_id = '2019072014'
  site_id = '2019_YELL_2'
  domain = 'D12'
elif site == "NIWO":
  search_criteria = 'NEON_D13_NIWO_DP1_20190814'
  date_id = '2019081416'
  site_id = '2019_NIWO_3'
  domain = 'D13'
elif site == "WREF":
  search_criteria = 'NEON_D16_WREF_DP1_20190712'
  date_id = '2019071216'
  site_id = '2019_WREF_3'
  domain = 'D16'
elif site == "TOOL":
  search_criteria = 'NEON_D18_TOOL_DP1_20190707'
  date_id = '2019070717'
  site_id = '2019_TOOL_3'
  domain = 'D18'
elif site == "PUUM":
  search_criteria = 'NEON_D20_PUUM_DP1_20190104'
  date_id = '2019010419'
  site_id = '2019_PUUM_1'
  domain = 'D20'

# List objects in the S3 bucket in the matching directory
dirpath = "NEON BRDF-TOPO Corrections/2019_" + site + "/"
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.json') and (search_criteria1 in obj['Key'])]
file_names = set()
for i,file in enumerate(files):
    match = re.search(r'DP1_(.*?)_reflectance', file)
    if match:
        file_name = match.group(1)
        print(file_name)
        file_names.add(file_name)
    else:
        print("Pattern not found in the URL.")
file_names = list(file_names)  # Convert set back to a list if needed
print(file_names)

# Loop through all site files
for i,file in enumerate(file_names):

    # Set to none to reduce memory use
    img = None
    neon = None
    topo_coeffs = None
    brdf_coeffs = None
    refl_md = None
    header_dict = None
    wavelength = None
    good_wl = None
    good_wl_list = None
    arrays = None
    fullarraystack = None
    ndvi = None
    mask = None
    
    print(file)
    flight = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D15/' + site_id + '/L1/Spectrometer/ReflectanceH5/' + date_id + '/NEON_' + domain + '_' + site + '_DP1_' + file +'_reflectance.h5'
    files = []
    files.append(flight)
    try:
        retrieve_neon_files(files, Data_Dir)
    except Exception as e:
        continue 
    img = Data_Dir + '/' + file + '_reflectance.h5'
    neon = ht.HyTools() 
    neon.read_file(img,'neon')
    print("file loaded")
    topo_file = "NEON BRDF-TOPO Corrections/2019_" + site + "/NEON_" + domain + "_" + site + "_DP1_" + file + "_reflectance_topo_coeffs_topo.json"
    print(topo_file)
    brdf_file = "NEON BRDF-TOPO Corrections/2019_" + site + "/NEON_" + domain + "_" + site + "_DP1_" + file + "_reflectance_brdf_coeffs_topo_brdf.json"
    s3.download_file(bucket_name, topo_file, Data_Dir + '/topo.json')
    s3.download_file(bucket_name, brdf_file, Data_Dir + '/brdf.json')
    print("Files downloaded successfully.")
    topo_coeffs = Data_Dir + "/topo.json"
    brdf_coeffs = Data_Dir + "/brdf.json"
    neon.load_coeffs(topo_coeffs,'topo')
    neon.load_coeffs(brdf_coeffs, 'brdf')
    print("corrections loaded")
    # Store map info for raster
    refl_md, header_dict = store_metadata(neon)
    # Export with corrections
    wavelength = header_dict['wavelength']
    good_wl = np.where((wavelength < 1340) | (wavelength > 1955), wavelength, np.nan)
    good_wl_list = good_wl[~np.isnan(good_wl)]
    print("creating arrays")
    arrays = [neon.get_wave(wave, corrections= ['topo','brdf'], mask = None) for wave in good_wl_list]
    print("stacking arrays")
    fullarraystack = np.dstack(arrays)
    print("Shape of fullarraystack:", fullarraystack.shape)
    print("calculating ndvi")
    ndvi = np.divide((fullarraystack[:, :, nir_band] - fullarraystack[:, :, red_band]), (fullarraystack[:, :, nir_band] + fullarraystack[:, :, red_band]), 
                     where=(fullarraystack[:, :, nir_band] + fullarraystack[:, :, red_band]) != 0)
    print("Shape of ndvi array:", ndvi.shape)
    # Apply NDVI threshold mask
    mask = ndvi < ndvi_threshold
    print("Shape of mask array:", mask.shape)
    print("masking by ndvi")
    fullarraystack[mask, :] = np.nan
    destination_s3_key = site + '_flightlines/'+ str(file)+'_output_' + '.tif'
    local_file_path = Out_Dir + '/output_fullarray_' + file + '.tif'
    print(local_file_path)
    print("rasterizing array")
    array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
    print("uploading array")
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)
    print("flightline complete")
