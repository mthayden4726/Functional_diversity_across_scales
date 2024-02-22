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
from S01_Functions import * # add scripts folder to python path manager
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

flights_D19_HEAL = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_223544_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_222709_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_222054_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_221439_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_220801_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_220134_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_215503_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_214808_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_214138_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_213445_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_212815_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_212119_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_211430_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_210758_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_210138_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_205506_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_204851_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_204220_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_203553_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_202906_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_202234_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_201535_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_200910_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_200156_reflectance.h5',
                    'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_195508_reflectance.h5'] # only 08 ones for now, not 06

# Trying parallelization
volumes = process_map(
        coeff_process,
        flights_D19_HEAL,
        max_workers=cpu_count() - 1
    )
print("Finished parallel process")

# Loop through all HEAL files
for i,file in enumerate(flights_D19_HEAL):
    match = re.search(r'DP1_(.*?)_reflectance', file)
    if match:
        file_name = match.group(1)
        print(file_name)
    else:
        print("Pattern not found in the URL.")
    files = []
    files.append(file)
    retrieve_neon_files(files, Data_Dir)
    img = Data_Dir + "/NEON_D19_HEAL_DP1_" + file_name + '_reflectance.h5'
    neon = ht.HyTools() 
    neon.read_file(img,'neon')
    print("file loaded")
    topo_file = "NEON BRDF-TOPO Corrections/2019_HEAL/NEON_D19_HEAL_DP1_" + file_name + "_reflectance_topo_coeffs_topo.json"
    print(topo_file)
    brdf_file = "NEON BRDF-TOPO Corrections/2019_HEAL/NEON_D19_HEAL_DP1_" + file_name + "_reflectance_brdf_coeffs_topo_brdf.json"
    try:
    # Attempt to download the file
        s3.download_file(bucket_name, topo_file, Data_Dir + '/topo.json')
        s3.download_file(bucket_name,brdf_file, Data_Dir + '/brdf.json')
        print("Files downloaded successfully.")
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
    destination_s3_key = 'HEAL_flightlines/'+ str(file_name)+'_output_' + '.tif'
    local_file_path = Out_Dir + '/output_fullarray_' + file_name + '.tif'
    print(local_file_path)
    array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    os.remove(local_file_path)
    print("flightline complete")
