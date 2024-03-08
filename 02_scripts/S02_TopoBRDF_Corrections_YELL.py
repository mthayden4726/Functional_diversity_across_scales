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
epsg = 32612

# Find correction coefficients (define search terms)
search_criteria = "NEON_D12_YELL_DP1_20190720"
dirpath = "NEON BRDF-TOPO Corrections/2019_YELL/"

# List objects in the S3 bucket in the matching directory
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.json') and (search_criteria in obj['Key'])]
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

file_names = ['20190720_165734', '20190720_192029', '20190720_221856', '20190720_175406', '20190720_225212', '20190720_174630', '20190720_223348', '20190720_220845', '20190720_171401', '20190720_163308', '20190720_184146', 
              '20190720_214100', '20190720_225949', '20190720_220340', '20190720_164916', '20190720_225602', '20190720_184916', '20190720_182608', '20190720_221346', 
              '20190720_224256', '20190720_222358', '20190720_173850', '20190720_223826', '20190720_215314', '20190720_214709', '20190720_170543', '20190720_185644', 
              '20190720_191342', '20190720_180157', '20190720_181808', '20190720_181005', '20190720_172212', '20190720_224726', '20190720_162450', '20190720_222854', 
              '20190720_215834', '20190720_183346', '20190720_173026', '20190720_190404', '20190720_161642']


# Loop through all UNDE files
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
    flight = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D12/2019_YELL_2/L1/Spectrometer/ReflectanceH5/2019072014/NEON_D12_YELL_DP1_' + file +'_reflectance.h5'
    files = []
    files.append(flight)
    try:
        retrieve_neon_files(files, Data_Dir)
    except Exception as e:
        continue 
    img = Data_Dir + "/NEON_D12_YELL_DP1_" + file + '_reflectance.h5'
    neon = ht.HyTools() 
    neon.read_file(img,'neon')
    print("file loaded")
    topo_file = "NEON BRDF-TOPO Corrections/2019_YELL/NEON_D12_YELL_DP1_" + file + "_reflectance_topo_coeffs_topo.json"
    print(topo_file)
    brdf_file = "NEON BRDF-TOPO Corrections/2019_YELL/NEON_D12_YELL_DP1_" + file + "_reflectance_brdf_coeffs_topo_brdf.json"
    try:
        s3.download_file(bucket_name, topo_file, Data_Dir + '/topo.json')
        s3.download_file(bucket_name, brdf_file, Data_Dir + '/brdf.json')
    except Exception as e:
        continue
    print("Files downloaded successfully.")
    topo_coeffs = Data_Dir + "/topo.json"
    brdf_coeffs = Data_Dir + "/brdf.json"
    neon.load_coeffs(topo_coeffs,'topo')
    neon.load_coeffs(brdf_coeffs, 'brdf')
    print("corrections loaded")
    # Store map info for raster
    refl_md, header_dict = store_metadata(neon, epsg)
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
    destination_s3_key = 'YELL_flightlines/'+ str(file)+'_output_' + '.tif'
    local_file_path = Out_Dir + '/output_fullarray_' + file + '.tif'
    print(local_file_path)
    print("rasterizing array")
    array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
    print("uploading array")
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)
    os.remove(img)
    print("flightline complete")
