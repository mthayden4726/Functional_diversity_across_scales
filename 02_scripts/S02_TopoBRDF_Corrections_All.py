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

flights_D17_TEAK = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_171251_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_173632_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_173103_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_172405_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_175502_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_174859_reflectance.h5',
             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D17/2019_TEAK_4/L1/Spectrometer/ReflectanceH5/2019061515/NEON_D17_TEAK_DP1_20190615_174242_reflectance.h5']

flights_D02_SERC = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D02/2019_SERC_4/L1/Spectrometer/ReflectanceH5/2019051512/NEON_D02_SERC_DP1_20190515_145907_reflectance.h5',
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

flights_D14_SRER = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_164806_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_165553_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_183556_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175509_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175923_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_182626_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_180750_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_171112_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_182142_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_170324_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_164014_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_174137_reflectance.h5',
               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_180334_reflectance.h5',
                       'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_175032_reflectance.h5',
                        'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_172621_reflectance.h5',
                         'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_184033_reflectance.h5',
                          'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_183106_reflectance.h5',
                           'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_181215_reflectance.h5',
                            'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_162527_reflectance.h5',
                             'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_171856_reflectance.h5',
                              'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_163317_reflectance.h5',
                               'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_160747_reflectance.h5',
                                'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_161733_reflectance.h5',
                                 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_181654_reflectance.h5',
                                  'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_180147_reflectance.h5',
                                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090114/NEON_D14_SRER_DP1_20190901_173359_reflectance.h5',
                                    'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_164418_reflectance.h5',
                                     'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_175358_reflectance.h5',
                                      'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_172901_reflectance.h5',
                                       'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090314/NEON_D14_SRER_DP1_20190903_163624_reflectance.h5']

flights_D19_HEAL = ['https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_202759_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_223544_reflectance.h5',
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
                    'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019081918/NEON_D19_HEAL_DP1_20190819_195508_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_213520_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_212633_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_211857_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_211126_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_210409_reflectance.h5',
                   'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D19/2019_HEAL_3/L1/Spectrometer/ReflectanceH5/2019062617/NEON_D19_HEAL_DP1_20190626_205654_reflectance.h5'] # still need to add a few more
# Combine into list of lists
#all_flights = [flights_D17_TEAK, flights_D02_SERC, flights_D14_SRER]

def topo_brdf_correct(all_flights):
  # Loop through all SERC files
  for m, site in enumerate(all_flights):
    match = re.search(r'flights_(.*?)', site)
    if match:
        site_name = match.group(1)
        print(site_name)
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
        img = Data_Dir + "/NEON_" + site_name + "_DP1_" + file_name + '_reflectance.h5'
        neon = ht.HyTools() 
        neon.read_file(img,'neon')
        print("file loaded")
        topo_file = "NEON BRDF-TOPO Corrections/2019_" + site_name + "/NEON_" + site_name + "_DP1_" + file_name + "_reflectance_topo_coeffs_topo.json"
        print(topo_file)
        brdf_file = "NEON BRDF-TOPO Corrections/2019_" + site_name + "/NEON_" + site_name + "_DP1_" + file_name + "_reflectance_brdf_coeffs_topo_brdf.json"
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
        refl_md, header_dict = store_metadata(neon)
        # Export with corrections
        wavelength = header_dict['wavelength']
        good_wl = np.where((wavelength < 1340) | (wavelength > 1955), wavelength, np.nan)
        good_wl_list = good_wl[~np.isnan(good_wl)]
        print("creating arrays")
        arrays = [neon.get_wave(wave, corrections= ['topo','brdf'], mask = None) for wave in good_wl_list]
        print("stacking arrays")
        fullarraystack = np.dstack(arrays)
        destination_s3_key = site_name + '_flightlines/'+ str(file_name)+'_output_' + '.tif'
        local_file_path = Out_Dir + '/output_fullarray_' + file_name + '.tif'
        print(local_file_path)
        array2rastermb(local_file_path, fullarraystack, refl_md, Out_Dir, epsg = refl_md['epsg'], bands = fullarraystack.shape[2])
        upload_to_s3(bucket_name, local_file_path, destination_s3_key)
        os.remove(local_file_path)
        print("flightline complete")
