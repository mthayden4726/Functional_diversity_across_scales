#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to perform BRDF and Topographical corrections on NEON flightlines using correction coefficients from Kovach et al. 

Author: M. Hayden
Updated: May 20, 2024

User input:
1. Name of the NEON site (e.g., BART)
2. Domain of the NEON site (e.g., D01)
3. ID of desired flights (e.g., 5)
4. Date of desired flights (e.g., 20190825)
5. Date and ID of desired flights (e.g., 2019082513)
6. EPSG of NEON site (e.g., 32619)
7. Desired NDVI threshold (e.g., 0.25)

"""

# Load required libraries
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
from S01_Functions import * # add scripts folder to python path manager
from osgeo import gdal, osr
from matplotlib.pyplot import subplots, show
import h5py
import sys
import rasterio
import boto3
import re
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
import argparse

# Set specifications 
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# Set global parameters
nir_band = 90
red_band = 58

# Use arg parse for local variables
# Create the parser
parser = argparse.ArgumentParser(description="Input script for BRDF/Topo corrections.")

# Add the arguments
parser.add_argument('--SITECODE', type=str, required=True, help='SITECODE (All caps)')
parser.add_argument('--DOMAIN', type=str, required=True, help='DOMAIN (D##)')
parser.add_argument('--ID_NO', type=str, required=True, help='ID (#)')
parser.add_argument('--DATE', type=str, required=True, help='YEAR (YYYYMMDD)')
parser.add_argument('--DATE_ID', type=str, required=True, help='YEAR (YYYYMMDD##)')
parser.add_argument('--EPSG', type=int, required=True, help='EPSG (#####)')
parser.add_argument('--NDVI', type=float, required=True, help='NDVI threshold(0-1)')
parser.add_argument('--NIR', type=float, required=True, help='NIR threshold(0-1)')

# Parse the arguments
args = parser.parse_args()

# Assign the arguments to variables
SITECODE = args.SITECODE
DOMAIN = args.DOMAIN
ID_NO = args.ID_NO
DATE = args.DATE
DATE_ID = args.DATE_ID
epsg = args.EPSG
ndvi_threshold = args.NDVI
nir_shade_threshold = args.NIR

SITE_STR = DOMAIN + '/' + '2019_' + SITECODE + '_' + ID_NO
SITE_STR_SHORT = DOMAIN + '_' + SITECODE

# Find correction coefficients (define search terms)
search_criteria = "NEON_" + DOMAIN + "_" + SITECODE + "_DP1_" + DATE
dirpath = "NEON BRDF-TOPO Corrections/2019_" + SITECODE + "/"

# List objects in the S3 bucket in the matching directory
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']

# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.json') and (search_criteria in obj['Key'])]
# Create empty set for file names
file_names = set()
# For each reflectance file, add ID to list
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

# Loop through all files in list created above
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

    # Retrieve and load reflectance as hytools object
    print(file)
    flight = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/'+ SITE_STR + '/L1/Spectrometer/ReflectanceH5/' + DATE_ID + '/NEON_' + SITE_STR_SHORT + '_DP1_' + file +'_reflectance.h5'
    files = []
    files.append(flight)
    # Download file
    try:
        retrieve_neon_files(files, Data_Dir)
    except Exception as e:
        continue 
    img = Data_Dir + '/NEON_' + SITE_STR_SHORT + '_DP1_' + file + '_reflectance.h5'
    neon = ht.HyTools() 
    neon.read_file(img,'neon')
    print("Flightline loaded")
    
    # Load correction coefficients
    topo_file = 'NEON BRDF-TOPO Corrections/2019_' + SITECODE + '/NEON_' + SITE_STR_SHORT + '_DP1_' + file + '_reflectance_topo_coeffs_topo.json'
    print(topo_file)
    brdf_file = 'NEON BRDF-TOPO Corrections/2019_' + SITECODE + '/NEON_' + SITE_STR_SHORT + '_DP1_' + file + '_reflectance_brdf_coeffs_topo_brdf.json'
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
    print("Corrections loaded")

    # Correct and export reflectance as a raster
    # Store map info for raster
    refl_md, header_dict = store_metadata(neon, epsg)
    # Remove bad bands
    wavelength = header_dict['wavelength']
    good_wl = np.where(((wavelength < 1340) | (wavelength > 1445)) & ((wavelength < 1790) | (wavelength > 1955)) & (wavelength >= 400) & (wavelength <= 2400), wavelength, np.nan)
    good_wl_list = good_wl[~np.isnan(good_wl)]
    # Perform corrections on each 'good' band
    print("Creating arrays")
    arrays = [neon.get_wave(wave, corrections= ['topo','brdf'], mask = None) for wave in good_wl_list]
    # Stack all corrected bands
    print("Stacking arrays")
    fullarraystack = np.dstack(arrays)
    print("Shape of corrected array:", fullarraystack.shape)
    
    # Perform radiometric corrections using NDVI threshold
    print("Calculating ndvi")
    ndvi = np.divide((fullarraystack[:, :, nir_band] - fullarraystack[:, :, red_band]), (fullarraystack[:, :, nir_band] + fullarraystack[:, :, red_band]), 
                     where=(fullarraystack[:, :, nir_band] + fullarraystack[:, :, red_band]) != 0)
    print("Shape of ndvi array:", ndvi.shape)
    # Apply NDVI threshold mask
    ndvi_mask = ndvi < ndvi_threshold
    print("Shape of NDVI mask array:", ndvi_mask.shape)

    # Add NIR-based shade mask (shade areas where NIR values are low)
    print("Creating NIR shade mask")
    nir_shade_mask = fullarraystack[:, :, nir_band] < nir_shade_threshold  # Define your shade threshold here
    print("Shape of NIR shade mask array:", nir_shade_mask.shape)

    # Combine NDVI mask and NIR shade mask
    combined_mask = ndvi_mask | nir_shade_mask  # Combine with OR condition
    print("Shape of combined mask array:", combined_mask.shape)

    # Apply the combined mask to the full array stack
    print("Masking by NDVI and NIR shade")
    fullarraystack[combined_mask, :] = np.nan  # Apply mask, set masked pixels to NaN
    
    print("masking by ndvi")
    fullarraystack[mask, :] = np.nan
    
    # Rasterize and export array to S3
    print("Rasterizing array")
    destination_s3_key = SITECODE + '_flightlines/'+ str(file)+'_output_' + '.tif'
    local_file_path = Out_Dir + '/output_fullarray_' + file + '.tif'
    print(local_file_path)
    array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
    print("Uploading array")
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)
    os.remove(img)
    print("Flightline complete. Next flightline!")
