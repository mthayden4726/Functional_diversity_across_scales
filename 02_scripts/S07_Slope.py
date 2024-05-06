#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 6 2024

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
import richdem as rd

#########################

# Set specifications 
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# global variables
SERVER = 'http://data.neonscience.org/api/v0/'

SITECODES = ['CLBJ']

search_criteria1 = "DTM"
search_criteria2 = "Mosaic"

summary_data = pd.DataFrame({'Site': [],
                             'Plot': [],
                             'Mean_slope': [],
                             'Median_slope': [],
                             'Max_slope': [],
                             'Min_slope': [],
                             'Std_slope': [],
                             'Var_slope': [],
                            })

#################################
# Pull elevation mosaics for calculating slope for each plot
for site in SITECODES:
    dirpath = "Environmental_Covariates/" + site + "/"
    file_stem = dirpath + "DTM_Mosaic_"
    objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
    # Filter objects based on the search criteria
    files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria1 in obj['Key']) and (search_criteria2 in obj['Key'])]
    file_names = set()
    for i,file in enumerate(files):
        match = re.search(r'DTM_Mosaic_(.*?).tif', file)
        if match:
            file_name = match.group(1)
            print(file_name)
            file_names.add(file_name)
        else:
            print("Pattern not found in the URL.")
    file_names = list(file_names)  # Convert set back to a list if needed
    print(file_names)
    for j in file_names:
        file_name = file_stem + str(j) + '.tif'
        print(file_name)
        # Download plot mosaic
        s3.download_file(bucket_name, file_name, Data_Dir + '/mosaic.tif')
        local_file_path = Data_Dir + '/mosaic.tif'
        with rasterio.open(local_file_path) as data_src:
          dem_data = data_src.read(1,masked=True)
        dem = rd.rdarray(dem_data, no_data=-9999)
        mask = (dem != -9999)
        slope = rd.TerrainAttribute(dem, attrib='slope_riserun')
        masked_slope = np.ma.masked_array(slope, mask=~mask)
        print(masked_slope.min())

        # initialize summaries
        data = [{'Site': site,
            'Plot': str(j), 
            'Mean_slope': slope.mean(),
            'Median_slope': np.median(slope),
            'Max_slope': slope.max(),
            'Min_slope': slope.min(),
            'Std_slope': slope.std(),
            'Var_slope': slope.var()
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
        destination_s3_key = 'Environmental_Covariates/Summary_files/' + site + '_Slope_Summary.csv'
        upload_to_s3(bucket_name, local_csv_path, destination_s3_key)
        print("File uploaded to S3")

        os.remove(local_csv_path)
