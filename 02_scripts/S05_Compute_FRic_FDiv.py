## Import packages ##
import hytools as ht
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import requests
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import kneed
from kneed import KneeLocator
import scipy.spatial
from scipy.spatial import ConvexHull
import subprocess
from urllib.request import urlretrieve
import multiprocessing as mp
import os, glob
import csv
import rasterio
from osgeo import gdal
import rioxarray as rxr
import xarray as xr
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep
import copy
import re
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
import pandas as pd
from fiona.crs import from_epsg
import pycrs
import csv
from csv import writer
# Import functions defined in S01_specdiv_functions.py
from S01_Functions import * # add scripts folder to python path manager
from S01_Moving_Window_FRIC import *
from S01_Moving_Window_FDiv import *

# Set directories
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/03_output'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

## Set global parameters ##
# window_sizes = [10, 30, 60, 120]   # list of window sizes to test
window_sizes = [60, 120, 240, 480, 960, 1200, 1500, 2000, 2200]
comps = 3 # default component numbers for PCA

# Loop through clipped files
# Choose site and plots
file_stem = 'HEAL_flightlines/Mosaic_HEAL_HEAL_'
plots = [
  '002',
  '004',
  '005',
  '013',
  '015',
  '018',
  '024',
  '026'
]

# Loop through plots
for i in plots:
    clip_file = file_stem + str(i) + '_v2.tif'
    print(clip_file)
    # Download plot mosaic
    s3.download_file(bucket_name, clip_file, Data_Dir + '/mosaic.tif')
    file = Data_Dir + '/mosaic.tif'
    # Open as raster
    raster = rxr.open_rasterio(file, masked=True)
    print(raster)
    # Convert data array to numpy array
    veg_np = raster.to_numpy()
    shape = veg_np.shape
    print(shape)
    # Flatten features into one dimesnion
    dim1 = shape[1]
    dim2 = shape[2]
    bands = shape[0]
    X = veg_np.reshape(bands,dim1*dim2).T
    print(X.shape)
    X = X.astype('float32')
    X[np.isnan(X)] = np.nan
    X[X < 0] = np.nan
    X = X/10000 # rescale data
    # Impute values for NAs
    imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    X_transformed = imputer.fit_transform(X)
    # Change nan to 0 
    #X_no_nan = np.nan_to_num(X, nan=0)
    # Take mean 
    x_mean = X_transformed.mean(axis=0)[np.newaxis, :]
    # Scale & standardize array
    X_transformed -=x_mean
    x_std = np.nanstd(X_transformed,axis=0)[np.newaxis, :]
    X_transformed /=x_std
    # Perform initial PCA fit
    pca = PCA(n_components=comps) # set max number of components
    pca.fit(X_transformed)
    pca_x =  pca.transform(X_transformed)
    print(pca_x)
    pca_x = pca_x.reshape((dim1, dim2,comps))
    print(pca_x.shape)
    # Calculate FRic on PCA across window sizes
    results_FR = {}
    local_file_path_fric = Out_Dir + "/HEAL_fric_" + str(i) + "_v2.csv"
    window_batches = [(a, pca_x, results_FR, local_file_path_fric) for a in np.array_split(window_sizes, cpu_count() - 1) if a.any()]
    volumes = process_map(
        window_calcs,
        window_batches,
        max_workers=cpu_count() - 1
    )
    destination_s3_key_fric = "/HEAL_fric_veg_" + str(i) + "_v2.csv"
    #f = open(local_file_path,"w")
    # write file
    #f.write(str(volumes))
    # close file
    #f.close()
    upload_to_s3(bucket_name, local_file_path_fric, destination_s3_key_fric)
    print("FRic file uploaded to S3")
    # Calculate FDiv on PCA across window sizes
    results_FD = {}
    local_file_path_fdiv = Out_Dir + "/HEAL_fdiv_veg_" + str(i) + "_v2.csv"
    window_batches = [(a, pca_x, results_FD, local_file_path_fdiv) for a in np.array_split(window_sizes, cpu_count() - 1) if a.any()]
    volumes = process_map(
        window_calcs_fdiv,
        window_batches,
        max_workers=cpu_count() - 1
    )
    # open file for writing
    destination_s3_key_fdiv = "/HEAL_fdiv_veg_" + str(i) + "_v2.csv"
    upload_to_s3(bucket_name, local_file_path_fdiv, destination_s3_key_fdiv)
    print("FDiv file uploaded to S3")
    os.remove(file)
    X = None
    X_no_nan = None
    pca_x = None
    veg_np = None
    print("Mosaic Complete - Next...")
