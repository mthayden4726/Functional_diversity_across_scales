## Import packages ##
import hytools as ht
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import requests
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
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
# Import functions defined in other scripts
from S01_Moving_Window_FEve.py import * # add scripts folder to python path manager
from S01_Functions.py import * # add scripts folder to python path manager

# Set directories
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/03_output'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

## Set global parameters ##
window_sizes = [60, 120, 240, 480, 700, 960, 1200, 1500, 2000, 2200]
ndvi_threshold = 0.4 # ndvi threshold for radiometric filtering
comps = 3 # default component numbers for PCA
# Other potential options (not currently included in this script):
# shade_threshold = 500 # should be low threshold across NIR region (15%)
# cloud_threshold = 1500 # should be high threshold in blue band
# bad_bands = [[300,400],[1300,1450],[1780,2000],[2450,2600]] # bands to be masked out
# sample_size = 0.1 # proportion of pixels to subsample for fitting PCA
# nclusters = 15 # default component numbers for K-means clustering

# Loop through clipped files
file_stem = 'SRER_flightlines/Mosaic_SRER_shape_'
sites = [3,4]
for i in sites:
    clip_file = file_stem + str(i) + '.tif'
    print(clip_file)
    s3.download_file(bucket_name, clip_file, Data_Dir + '/mosaic.tif')
    file = Data_Dir + '/mosaic.tif'
    raster = rxr.open_rasterio(file, masked=True)
    print(raster)
    # Array-ize, transform with PCA
    pca_x = pca_steps(raster, comps)
    print(pca_x)
    print(pca_x.shape[0])
    # Paralellize calcs for different window sizes
    results_FE = {}
    local_file_path = Out_Dir + "/SRER_feve_" + str(i) + ".csv"
    window_batches = [(a, pca_x, results_FE, local_file_path) for a in np.array_split(window_sizes, cpu_count() - 1) if a.any()]
    volumes = process_map(
        window_calcs_feve,
        window_batches,
        max_workers=cpu_count() - 1
    )
    #print(volumes)
    # open file for writing
    # local_file_path = Out_Dir + "/TEAK_fric_" + str(i) + ".csv"
    destination_s3_key = "/SRER_feve_" + str(i) + ".csv"
    #f = open(local_file_path,"w")
    # write file
    #f.write(str(volumes))
    # close file
    #f.close()
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    print("File uploaded to S3")
    os.remove(file)
    print("Mosaic Complete - Next...")
