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
from window_calcs import * # add scripts folder to python path manager
from S01_specdiv_functions import * # add scripts folder to python path manager
from osgeo import gdal, osr
from matplotlib.pyplot import subplots, show
import h5py
import sys
import rasterio
import boto3

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'

# Array to raster
def array2raster(newRaster, reflBandArray, reflArray_metadata, Out_Dir, epsg):
    NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,
    }
    pwd = os.getcwd()
    os.chdir(Out_Dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = 1
    pixelWidth = float(refl_md['res']['pixelWidth'])
    pixelHeight = -float(refl_md['res']['pixelHeight'])
    originX = refl_md['ext_dict']['xMin']
    originY = refl_md['ext_dict']['yMax']
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    #for band in range(bands):
    #    outRaster.GetRasterBand(band + 1).WriteArray(array[:, :, band])
    outRaster.WriteArray(reflBandArray[:, :])
    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)
    
# Select files to pull from NEON
SITECODE = 'TEAK' # NEON site of interest
PRODUCTCODE = 'DP1.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2019-06'
#file_paths = find_neon_files(SITECODE,
#                             PRODUCTCODE,
#                          YEAR)
#file_paths_dates = glob.glob("https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/20190615*" + "/NEON_D17_TEAK_DP1_20190615*")
# Pull from NEON
#flights_TEAK = retrieve_neon_files(file_paths, Data_Dir)

# Read in image
neon_image = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata/NEON_D17_TEAK_DP1_20190615_171251_reflectance.h5'

neon = ht.HyTools()
neon.read_file(neon_image,'neon')

# Load correction coefficients
topo_coeffs = '/home/ec2-user/BioSCape_across_scales/home/BioSCape_across_scales/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_171251_reflectance_topo_coeffs_topo.json'
brdf_coeffs = '/home/ec2-user/BioSCape_across_scales/home/BioSCape_across_scales/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_171251_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')
neon.load_coeffs(topo_coeffs,'topo')

# Create dictionary containing metadata for creation of raster
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands] 
#Extract no data value & scale factor
#refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['noDataVal'] = -9999
refl_md['scaleFactor'] = 10000
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32611 # for wgs 84, UTM 11N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax

# Array to raster)
def array2raster(newRaster, reflBandArray, reflArray_metadata, ras_dir, epsg):
    NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,
    }
#
    pwd = os.getcwd()
    os.chdir(Out_Dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = 1
    pixelWidth = float(refl_md['res']['pixelWidth'])
    pixelHeight = -float(refl_md['res']['pixelHeight'])
    originX = refl_md['ext_dict']['xMin']
    originY = refl_md['ext_dict']['yMax']
#
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    # for band in range(bands):
    #    outRaster.GetRasterBand(band + 1).WriteArray(array[:, :, band])
    outRaster.WriteArray(reflBandArray[:, :])
    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)

def upload_to_s3(bucket_name, file_path, s3_key):
    """
    Upload a file from an EC2 instance to Amazon S3.

    :param bucket_name: Name of the S3 bucket
    :param file_path: Local path to the file on the EC2 instance
    :param s3_key: Destination key in the S3 bucket (e.g., folder/file_name.ext)
    """
    # Initialize the S3 client
    s3 = boto3.client('s3')
                                        
    try:
    # Upload the file
        s3.upload_file(file_path, bucket_name, s3_key)
        print(f'Successfully uploaded {file_path} to {bucket_name}/{s3_key}')
    except Exception as e:
        print(f"Error uploading file: {e}")


# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    print(i)
    refl_md['wavelength'] = i
    local_file_path = Out_Dir + '/TEAK_flightlines/20190615_171251_output_' + str(i) + '.tif'
    destination_s3_key = 'TEAK_flightlines/20190615_171251_output_' + str(i) + '.tif'
    array2raster(local_file_path, arrays[i], refl_md, ras_dir = Out_Dir, epsg = 32611)
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)































=======
# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190615_171836_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)

# Read in image #2
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata//NEON_D17_TEAK_DP1_20190615_171251_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_171251_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_171251_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32612 # for wgs 84, UTM 12N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax
>>>>>>> main

# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
print("Creating next set of arrays")
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190615_171251_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)
    
# Image 3
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata//NEON_D17_TEAK_DP1_20190617_185346_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_185346_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
#brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_185346_reflectance_brdf_coeffs_topo_brdf.json'
#neon.load_coeffs(brdf_coeffs,'brdf')

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32612 # for wgs 84, UTM 12N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax

# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
print("Creating next set of arrays")
arrays = [neon.get_wave(wave, corrections= ['topo']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190617_185346_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)
    
# Image 4
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata//NEON_D17_TEAK_DP1_20190617_172159_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_172159_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_172159_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32612 # for wgs 84, UTM 12N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax

# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
print("Creating next set of arrays")
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190617_172159_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)

# Image 5
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata//NEON_D17_TEAK_DP1_20190617_164750_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_164750_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190617_164750_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32612 # for wgs 84, UTM 12N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax

# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
print("Creating next set of arrays")
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190617_164750_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)
    
    
# Image 6
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata//NEON_D17_TEAK_DP1_20190615_181119_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_181119_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_TEAK/NEON_D17_TEAK_DP1_20190615_181119_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 32612 # for wgs 84, UTM 12N
refl_md['res'] = {}
refl_md['res']['pixelWidth'] = float(mapInfo[5])
refl_md['res']['pixelHeight'] = float(mapInfo[6])
# Extract the upper left-hand corner coordinates from mapInfo
xMin = float(mapInfo[3])  # convert from string to floating point number
yMax = float(mapInfo[4])
# Calculate the xMax and yMin values from the dimensions
xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
refl_md['ext_dict'] = {}
refl_md['ext_dict']['xMin'] = xMin
refl_md['ext_dict']['xMax'] = xMax
refl_md['ext_dict']['yMin'] = yMin
refl_md['ext_dict']['yMax'] = yMax

# Loop - requires refl_md be established for neon file
wavelength = header_dict['wavelength']
print("Creating next set of arrays")
arrays = [neon.get_wave(wave, corrections= ['topo','brdf']) for wave in wavelength]
for i in range(0,425):
    refl_md['wavelength'] = i
    newRaster = Out_Dir + '/TEAK_flightlines/20190615_181119_output_' + str(i) + '.tif'
    array2raster(newRaster, arrays[i], refl_md, Out_Dir, epsg = 32612)


# batch process arrays
wave_batches = [(a,neon) for a in np.array_split(wavelength, cpu_count() - 1) if a.any()]

# sequential processing for debugging
# for batch in window_batches:
#     window_calcs(batch)
process_map(
    wave_calcs,
    wave_batches,
    max_workers=cpu_count() - 1)

# view output
raster = rasterio.open(Out_Dir + "/output_2.tif")
plt.imshow(raster.read(1), cmap="BrBG")
plt.title("Temperature")
plt.show()


import earthpy.spatial as es
import earthpy.plot as ep

stack_band_paths = glob.glob(Out_Dir + "/output_*")
raster_out_path = os.path.join(Out_Dir, "raster.tiff")
with rasterio.Env(CHECK_DISK_FREE_SPACE="FALSE"):
    array, raster_prof = es.stack(stack_band_paths, out_path=raster_out_path)
    
# can also merge these in QGIS if necessary using One Click Raster Stacking plug-in
# super fast, creates lightweight stack
# then re-load rasters to mosaic??

extent = plotting_extent(array[0], raster_prof["transform"])

fig, ax = plt.subplots(figsize=(6, 6))
ep.plot_rgb(
    array,
    ax=ax,
    rgb=[3, 2, 1],
    stretch=True,
    extent=extent,
    str_clip=0.5,
    title="RGB Image of Un-cropped Raster",
)
plt.show()
