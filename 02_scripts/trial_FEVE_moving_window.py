import numpy as np
from multiprocessing import Pool
from S01_Moving_Window_FEve import *
from S01_Functions import *
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

def calculate_metrics(starting_point, raster, window_size):
    # Calculate metrics for a given starting point and window size
    # This function should include the logic for calculating Primm's evenness,
    # partial weighted evenness, and functional evenness for the specified window
    # around the starting point
  
    i, j = starting_point
    half_window = window_size // 2
    window = raster[max(0, i - half_window):min(raster.shape[0], i + half_window + 1),
                    max(0, j - half_window):min(raster.shape[1], j + half_window + 1)]
    
    # Calculate metrics for the window
    PEW, S = primms_mst(window)
    FEve = calculate_FEve_villager(PEW, S)

    return starting_point, FEve

def process_maps(raster, window_size):
    results = []
    pool = Pool(processes= cpu_count() - 1)
    
    for starting_point in generate_starting_points(raster, window_size):
        # Submit a job to the pool for each starting point
        results.append(pool.apply_async(calculate_metrics, (starting_point,)))
    
    pool.close()
    pool.join()

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Starting_Point', 'FEve']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        
        # Get results from worker processes
        for result in results:
            starting_point, metrics = result.get()
            writer.writerow({'Starting_Point': starting_point, 'FEve': metrics})

def generate_starting_points(raster, window_size):
    # Generate starting points for the moving window analysis
    # This function should iterate over the raster and yield starting points
    # based on the desired window size
    
    # Placeholder for demonstration purposes
    half_window = window_size // 2
    for i in range(0, pca_x.shape[0], window_size):
        for j in range(0, pca_x.shape[1], window_size):
            yield (i, j)

window_size = 60
s3.download_file(bucket_name, 'TEAK_flightlines/Mosaic_clip_site_0.tif', Data_Dir + '/mosaic.tif')
file = Data_Dir + '/mosaic.tif'
raster = rxr.open_rasterio(file, masked=True)
print(raster)
pca_x = pca_steps(raster, comps)
print(pca_x)
output_file = 'FEve_results.csv'
process_maps(pca_x, window_size, output_file)
