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
import re
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

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

def array2rastermb(newRaster, reflBandArray, reflArray_metadata, Out_Dir, epsg, bands):
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
    print(pwd)
    os.chdir(Out_Dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    print(cols,rows)
    pixelWidth = float(refl_md['res']['pixelWidth'])
    pixelHeight = -float(refl_md['res']['pixelHeight'])
    originX = refl_md['ext_dict']['xMin']
    originY = refl_md['ext_dict']['yMax']
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    print(outRaster)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    #outband = outRaster.GetRasterBand(1)
    #outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        print(band)
        outRaster.GetRasterBand(band + 1).WriteArray(reflBandArray[:, :, band])
    #outRaster.WriteArray(reflBandArray[:, :, :])
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

# Select files to pull from NEON

# Select files to pull from NEON
#SITECODE = 'SERC' # NEON site of interest
#PRODUCTCODE = 'DP1.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
#YEAR = '2019-05'
#file_paths = find_neon_files(SITECODE,
#                             PRODUCTCODE,
#                          YEAR)
#print(file_paths)

flights_SERC = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_145907_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_150424_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_150934_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_151446_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_152015_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_152538_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_153118_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_153658_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_154239_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_154747_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_155313_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_155853_reflectance.h5',                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_160426_reflectance.h5',
                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_161116_reflectance.h5']

#flights_SRER = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_164806_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_165553_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_183556_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175509_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175923_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_182626_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_180750_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_171112_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_182142_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_170324_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_164014_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_174137_reflectance.h5',
#               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_180334_reflectance.h5',
#                       'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175032_reflectance.h5',
#                        'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_172621_reflectance.h5',
#                         'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_184033_reflectance.h5',
#                          'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_183106_reflectance.h5',
#                           'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_181215_reflectance.h5',
#                            'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_162527_reflectance.h5',
#                             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_171856_reflectance.h5',
#                              'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_163317_reflectance.h5',
#                               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_160747_reflectance.h5',
#                                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_161733_reflectance.h5',
#                                 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_181654_reflectance.h5',
#                                  'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_180147_reflectance.h5',
#                                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_173359_reflectance.h5',
#                                    'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_164418_reflectance.h5',
#                                     'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_175358_reflectance.h5',
#                                      'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_172901_reflectance.h5',
#                                       'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_163624_reflectance.h5',
#                                        'https://storage.googleapis.com/neon-aop-products/2019/Fu

# Loop through all SERC files
for i,file in enumerate(flights_SERC):
    match = re.search(r'DP1_(.*?)_reflectance', file)
    if match:
        file_name = match.group(1)
        print(file_name)
    else:
        print("Pattern not found in the URL.")
    files = []
    files.append(file)
    retrieve_neon_files(files, Data_Dir)
    img = Data_Dir + "/NEON_D02_SERC_DP1_" + file_name + '_reflectance.h5'
    neon = ht.HyTools() 
    neon.read_file(img,'neon')
    print("file loaded")
    topo_file = "NEON BRDF-TOPO Corrections/2019_SERC/NEON_D02_SERC_DP1_" + file_name + "_reflectance_topo_coeffs_topo.json"
    print(topo_file)
    brdf_file = "NEON BRDF-TOPO Corrections/2019_SERC/NEON_D02_SERC_DP1_" + file_name + "_reflectance_brdf_coeffs_topo_brdf.json"
    try:
    # Attempt to download the file
        s3.download_file(bucket_name, topo_file, Data_Dir + '/topo.json')
        print("File downloaded successfully.")
    except FileNotFoundError:
        print("The file does not exist in the specified S3 bucket.")
    except NoCredentialsError:
        print("AWS credentials could not be found.")
    except PartialCredentialsError:
        print("AWS credentials are incomplete.")
    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            print("The object does not exist in the S3 bucket.")
        else:
            print("An error occurred:", e)
    except Exception as e:
        print("An unexpected error occurred:", e)
    #s3.meta.client.download_file(bucket_name,topo_file, Data_Dir + 'topo.json')
    s3.download_file(bucket_name,brdf_file, Data_Dir + '/brdf.json')
    topo_coeffs = Data_Dir + "/topo.json"
    brdf_coeffs = Data_Dir + "/brdf.json"
    neon.load_coeffs(topo_coeffs,'topo')
    neon.load_coeffs(brdf_coeffs, 'brdf')
    print("corrections loaded")
    # Store map info for raster
    mapInfo= neon.map_info
    header_dict = neon.get_header()
    refl_md = {}
    refl_md['mapInfo'] = header_dict['map info']
    refl_md['wavelength'] = header_dict['wavelength']
    refl_md['shape'] = [neon.lines, neon.columns, neon.bands]
    #Extract no data value & scale factor
    refl_md['noDataVal'] = float(header_dict['data ignore value'])
    refl_md['scaleFactor'] = float(0.996)
    refl_md['bad_band_window1'] = np.array([1340, 1445])
    refl_md['bad_band_window2'] = np.array([1790, 1955])
    refl_md['epsg'] = 32618 # for wgs 84, UTM 11N --> note that this changed by site!! 12N for SRER
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
    # Export with corrections
    wavelength = header_dict['wavelength']
    good_wl = np.where((wavelength < 1340) | (wavelength > 1955), wavelength, np.nan)
    good_wl_list = good_wl[~np.isnan(good_wl)]
    print("creating arrays")
    arrays = [neon.get_wave(wave, corrections= ['topo','brdf'], mask = None) for wave in good_wl_list]
    print("stacking arrays")
    fullarraystack = np.dstack(arrays)
    destination_s3_key = 'SERC_flightlines/'+ str(file_name)+'_output_' + '.tif'
    local_file_path = Out_Dir + '/output_fullarray_' + file_name + '.tif'
    print(local_file_path)
    array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)
    print("flightline complete")


# if the array does not have the right dimensions, use --->
#with rasterio.open('image.tif', 'w', **meta) as dst:
#    dst.write(fullarray)
#    # If array is in (x, y, z) order (cols, rows, bands)
#    dst.write(np.moveaxis(array, [0, 1, 2], [2, 1, 0]))

# mosaic flightlines together
# virutally? too large otherwise?
# gdalbuildvrt  mosaic.vrt *.bsq
# Try with rasterio
#tiff_files = glob.glob(###)
#src_files_to_mosaic = []
#for file in tiff_files:
#   src = rasterio.open(file)
#   src_files_to_mosaic.append(src)
#mosaic, out_trans = merge(src_files_to_mosaic)

# Then clip according to area around each field site
# feed into fric_iterate script
