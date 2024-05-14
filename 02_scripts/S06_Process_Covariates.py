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
import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import fiona
import glob
import os
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
import geopandas as gpd
import pandas as pd
from fiona.crs import from_epsg
import argparse

#########################

# Set specifications 
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# global variables
SERVER = 'http://data.neonscience.org/api/v0/'

# experimenting with arg parse
# Create the parser
parser = argparse.ArgumentParser(description="Input script for environmental analysis.")

# Add the arguments
parser.add_argument('--SITECODE', type=str, required=True, help='SITECODE (All caps)')
parser.add_argument('--DOMAIN', type=str, required=True, help='DOMAIN (D##)')
parser.add_argument('--ID_NO', type=int, required=True, help='ID (#)')
parser.add_argument('--YEAR', type=str, required=True, help='YEAR (YYYY-MM)')
parser.add_argument('--ENV', type=str, required=True, choices=['DTM', 'CHM'], help='Environmental Covariate of Interest (DTM or CHM)')

# Parse the arguments
args = parser.parse_args()

# Assign the arguments to variables
SITECODE = args.SITECODE
DOMAIN = args.DOMAIN
ID_NO = args.ID_NO
YEAR = args.YEAR
ENV = args.ENV

# prompt script user to provide input instead
#SITECODE = input("SITECODE (All caps)")
#DOMAIN = input("DOMAIN (D##)")
#ID_NO = input("ID (#)")
#YEAR = input("YYYY-MM")
#ENV = input("Environmental Covariate of Interest (DTM, CHM, or slope)")
#shapefiles = input("List of shapefiles (['00#',...]")

# assign local variables

SITE_STR = DOMAIN + '/' + '2019_' + SITECODE + '_' + ID_NO
SITE_STR_SHORT = DOMAIN + '_' + SITECODE

if ENV == "DTM":
  PRODUCTCODE = 'DP3.30024.001'
  ENV_lab = "DTM"
elif ENV == "CHM":
  PRODUCTCODE = 'DP3.30015.001'
  ENV_lab = "CanopyHeightModel"
elif ENV == "slope":
  PRODUCTCODE = 'DP3.30025.001'
  ENV_lab = "Slope"

#################################

# Find and load files
url = SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/'+YEAR
# Request the url
data_request = requests.get(url)
# Convert the request to Python JSON object
data_json = data_request.json()
print(data_json)
# Create list of file paths of interest for given site, product, year
file_paths = []
file_string = ENV + ".tif"
for file in data_json['data']['files'][:]:
  if file_string in file['name']:
    file_paths.insert(1, file['url'])
    #print(file['url'])
print(file_paths)

file_names = set()
for i,file in enumerate(file_paths):
    match = re.search(r'DP3_(.*?).tif', file)
    if match:
        file_name = match.group(1)
        print(file_name)
        file_names.add(file_name)
    else:
        print("Pattern not found in the URL.")
file_names = list(file_names)  # Convert set back to a list if needed
print(file_names)

for i, file in enumerate(file_names):
  print(file)
  flight = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/' + SITE_STR + '/L3/DiscreteLidar/' + ENV_lab + 'Gtif/NEON_' + SITE_STR_SHORT + '_DP3_' + file +'.tif'
  files = []
  files.append(flight)
  retrieve_neon_files(files, Data_Dir)
  local_file_path = Data_Dir + "/NEON_" + SITE_STR_SHORT + "_DP3_" + file + '.tif'
  destination_s3_key = 'Environmental_Covariates/' + SITECODE + '/' + ENV + '_' + str(file) + '.tif'
  
  upload_to_s3(bucket_name, local_file_path, destination_s3_key)
  os.remove(local_file_path)
  print("flightline complete")

#################### 
# Clip files to shapefiles
search_criteria = ENV
dirpath = "Environmental_Covariates/" + SITECODE + "/"

# List objects in the S3 bucket in the matching directory
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria in obj['Key'])]
print(files)

# Clip files
for j,shape in enumerate(shapefiles):
    print(shape)
    ID = "Site_boundaries/" + SITECODE + "/" + SITECODE + "_" + str(shape)
    downloaded_files = download_shapefile(bucket_name, ID, Out_Dir)
    shapefile_path = next(file for file in downloaded_files if file.endswith('.shp'))
    
    # Open shapefile and access geometry
    with fiona.open(shapefile_path, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    print(shapes)
    
    for i, file in enumerate(files):
        print('Loading file from S3')
        s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
        flight = Out_Dir + '/file_' + str(i) + '.tif'
        print(flight)
        with rasterio.open(flight) as src:
            try:
                out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                out_meta = src.meta
                print("File Clipped")
                print(out_meta)
                out_meta.update({"driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform})
                local_file_path = Out_Dir + "/clip.tif"
                with rasterio.open(local_file_path, "w", **out_meta) as dest:
                    dest.write(out_image)
                    print("File Written")
                destination_s3_key = 'Environmental_Covariates/' + SITECODE + '/' + ENV + '_' + str(shape) + '_Clipped_file_' + str(i) + '.tif'
                upload_to_s3(bucket_name, local_file_path, destination_s3_key)
                print("File uploaded to S3")
                os.remove(local_file_path)
            except ValueError as e:
                # Handle the case where there is no overlap between the raster and the shapefiles
                print(f"Skipping file {i} as it does not overlap with the shapefile.")
            os.remove(flight)
        print("Cleared data files")

#################### 
# Mosaic files and produce output

src_files_to_mosaic = []
summary_data = pd.DataFrame({'Site': [],
                             'Plot': [],
                             'Env': [],
                             'Mean': [],
                             'Median': [],
                             'Max': [],
                             'Min': [],
                             'Std': [],
                             'Var': [],
                            })

for i,ID in enumerate(shapefiles):
    src_files_to_mosaic = []
    # List files associated with a single buffer shape
    search_criteria1 = str(ID)
    search_criteria2 = ENV
    search_criteria3 = "Clipped"
    dirpath = "Environmental_Covariates/" + SITECODE

    # List objects in the S3 bucket in the matching directory
    objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
    # Filter objects based on the search criteria
    files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria1 in obj['Key'])
            and (search_criteria2 in obj['Key']) and (search_criteria3 in obj['Key'])]
    print(files)
    for j,file in enumerate(files):
        flight  = Out_Dir + '/file_' + str(j) + '.tif'
        try:
            s3.download_file(bucket_name, file, flight)
            print(f"The file '{file}' exists.")
        except Exception as e:
            print(f"Error: {e}")
        src = rasterio.open(flight)
        src_files_to_mosaic.append(src)

    # Mosaic files
    print(src_files_to_mosaic)
    mosaic, out_trans = merge(src_files_to_mosaic, method = 'max')
    print('Merge complete')
    # Update metadata
    out_meta = src.meta.copy()
    print(out_meta)
    print(mosaic.shape)
    out_meta.update({"driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans})
    print(out_meta)

    # Write to computer, send to S3
    local_file_path = Out_Dir + "/mosaic_" + SITECODE + ".tif"
    with rasterio.open(local_file_path, "w", **out_meta) as dest:
        dest.write(mosaic)
    print("File written")

    # Push to S3 bucket
    destination_s3_key = 'Environmental_Covariates/' + SITECODE + '/' + ENV + '_Mosaic_' +str(ID)+ '.tif'
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    print("File uploaded to S3")

    # Get summary metrics
    with rasterio.open(local_file_path) as data_src:
      env_data = data_src.read(1, masked=True)

    # initialize summaries
    data = [{'Site': SITECODE,
             'Plot': str(ID), 
             'Env': ENV,
             'Mean': env_data.mean(),
             'Median': np.median(env_data),
             'Max': env_data.max(),
             'Min': env_data.min(),
             'Std': env_data.std(),
             'Var': env_data.var()
                            }]
    # append to dataframe
    summary_data = summary_data.append(data, ignore_index=True)
    print(summary_data)
    
    # Remove unneeded files (mosaic and shapefile)
    os.remove(local_file_path)
  
    mosaic = None
    src_files_to_mosaic = None 

# Create the pandas DataFrame
print(summary_data)
summary_data.to_csv('summary.csv')
local_csv_path = 'summary.csv'
    
# Push to S3 bucket
destination_s3_key = 'Environmental_Covariates/Summary_files/' + SITECODE + '_' + ENV + '_Summary.csv'
upload_to_s3(bucket_name, local_csv_path, destination_s3_key)
print("File uploaded to S3")

os.remove(local_csv_path)
