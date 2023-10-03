#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:22:47 2023

@author: meha3816
"""
sys.path.append("/Users/local/hdf5/lib/")
## Load necessary packages ##
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
print(sys.path)

Data_Dir = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata'

SITECODE = 'SRER' # NEON site of interest
PRODUCTCODE = 'DP1.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2019-09'

files = os.listdir(Data_Dir)
files_TEAK = glob.glob(Data_Dir + "/NEON_D17*")
files_SRER = glob.glob(Data_Dir + "/NEON_D14_SRER_DP1*")
files_TALL = glob.glob(Data_Dir + "/NEON_D08*")
files_KONZ = glob.glob(Data_Dir + "/NEON_D06*") # never downloaded, too slow
files_SERC = glob.glob(Data_Dir + "/NEON_D02*") # never downloaded, too slow



file_paths = find_neon_files(SITECODE,
                             PRODUCTCODE,
                          YEAR)

# count proportion of cells in scene with ndvi == 0
list = []
for file in files_SRER:
    neon = ht.HyTools() 
    neon.read_file(file,'neon')
    # show_rgb(neon, r=660,g=550,b=440, correct = [])
    ndvi = neon.ndi()
    zero_els = np.count_nonzero(ndvi==0)
    list.append(file)
    list.append((zero_els/(1000*1000))*100)
    #plt.imshow(ndvi)
    #plt.show()

# plot scenes 
fig, axs = plt.subplots(nrows=1, ncols=5)
for file, ax in zip(files_SRER, axs.ravel()):
    neon = ht.HyTools() 
    neon.read_file(file,'neon')
    # show_rgb(neon, r=660,g=550,b=440, correct = [])
    ndvi = neon.ndi()
    # filter df for ticker and plot on specified axes
    # show image
    shw = ax.imshow(ndvi)
      
    # make bar
    bar = plt.colorbar(shw)
      
    # show plot with labels
    plt.xlabel('Meters')
    plt.ylabel('Meters')
    
plt.imshow(ndvi)

# For a single plot
# make plot
#fig, ax = plt.subplots()
  
# show image
#shw = ax.imshow(ndvi)
  
# make bar
#bar = plt.colorbar(shw)
  
# show plot with labels
#plt.xlabel('Meters')
#plt.ylabel('Meters')
#bar.set_label('NDVI')
#plt.show()

# Looking at correction changes
fp = file_paths[6:7]
flights_SRER = retrieve_neon_files(file_paths, Data_Dir)

neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON_D14_SRER_DP1_20190901_165553_reflectance.h5'

neon = ht.HyTools()
neon.read_file(neon_image,'neon')

topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_SRER/NEON_D14_SRER_DP1_20190901_165553_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')

brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_SRER/NEON_D14_SRER_DP1_20190901_165553_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')

rgb = show_rgb(neon, correct = [])
rgbcorr = show_rgb(neon, correct = ['topo', 'brdf'])

difference = neon.get_wave(1150) - neon.get_wave(1150,corrections= ['topo',
                                                                  'brdf'])
diff_sub = difference[300:500, 300:500]
ndvi_sub = ndvi[300:500, 300:500]
plt.matshow(difference)

# Run next four lines together for plot
fig, ax = plt.subplots()
shw = ax.imshow(diff_sub)
bar = plt.colorbar(shw)
plt.show()

rgb_sub = rgb[300:500, 300:500, ]
plt.imshow(rgb_sub)

pixel = neon.get_chunk(0,100,0,100, corrections = ['topo','brdf'])
plt.imshow(pixel)

def show_rgb(hy_obj,r=660,g=550,b=440, correct= []):

    rgb=  np.stack([hy_obj.get_wave(r,corrections= correct),
                    hy_obj.get_wave(g,corrections= correct),
                    hy_obj.get_wave(b,corrections= correct)])
    rgb = np.moveaxis(rgb,0,-1).astype(float)
    rgb[rgb ==hy_obj.no_data] = np.nan

    bottom = np.nanpercentile(rgb,5,axis = (0,1))
    top = np.nanpercentile(rgb,95,axis = (0,1))
    rgb = np.clip(rgb,bottom,top)

    rgb = (rgb-np.nanmin(rgb,axis=(0,1)))/(np.nanmax(rgb,axis= (0,1))-np.nanmin(rgb,axis= (0,1)))

    height = int(hy_obj.lines/hy_obj.columns)

    fig  = plt.figure(figsize = (7,7) )
    plt.imshow(rgb)
    plt.show()
    plt.close()
    return(rgb)

## Mosaicing files together
files_to_mosaic = glob.glob(Data_Dir + '/NEON_D14_SRER_DP1_*')
files_to_mosaic

files_string = " ".join(files_to_mosaic)
print(files_string)

command = "gdal_merge.py -o /Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON_D14_SRER_DP1_mosaic.tiff" + files_string
print(os.popen(command).read())

# use neon functions
import numpy as np
import matplotlib.pyplot as plt
import h5py, gdal, osr, copy

srerRefl, srerRefl_md, array = h5refl2array('/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON_D14_SRER_DP1_20190901_165553_reflectance.h5')
array2raster('SRER.tif',srerRefl,srerRefl_md)


def h5refl2array(refl_filename):
    """h5refl2array reads in a NEON AOP reflectance hdf5 file and returns the
    reflectance array, select metadata, and wavelength dataset.
    --------
    Parameters
        refl_filename -- full or relative path and name of reflectance hdf5 file
    --------
    Returns 
    --------
    reflArray:
        array of reflectance values
    metadata:
        dictionary containing the following metadata:
            EPSG: coordinate reference system code (integer)
            *bad_band_window1: [1340 1445] range of wavelengths to ignore
            *bad_band_window2: [1790 1955] range of wavelengths to ignore 
            ext_dict: dictionary of spatial extent 
            extent: array of spatial extent (xMin, xMax, yMin, yMax)
            mapInfo: string of map information 
            *noDataVal: -9999.0
            projection: string of projection information
            *res: dictionary containing 'pixelWidth' and 'pixelHeight' values 
            *scaleFactor: 10000.0
            shape: tuple of reflectance shape (y, x, # of bands)
        * Asterixed values are the same for all NEON AOP hyperspectral 
        reflectance files processed 2016 & after.
    wavelengths:
        Wavelengths dataset. This is the same for all NEON AOP reflectance hdf5 files.
        wavelengths.value[n-1] gives the center wavelength for band n 
    --------
    This function applies to the NEON hdf5 format implemented in 2016, which 
    applies to data acquired in 2016 & 2017 as of June 2017. Data in earlier 
    NEON hdf5 format is will be re-processed after the 2017 flight season. 
    --------
    Example
    --------
    sercRefl, sercRefl_md, wavelengths = h5refl2array('NEON_D02_SERC_DP1_20160807_160559_reflectance.h5') """
    
    #Read in reflectance hdf5 file 
    #include full or relative path if data is located in a different directory
    hdf5_file = h5py.File(refl_filename,'r')

    #Get the site name
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]
    
    #Extract the reflectance & wavelength datasets
    refl = hdf5_file[sitename]['Reflectance']
    reflArray = refl['Reflectance_Data']
    refl_shape = reflArray.shape
    wavelengths = refl['Metadata']['Spectral_Data']['Wavelength']
    
    #Create dictionary containing relevant metadata information
    metadata = {}
    metadata['shape'] = reflArray.shape
    metadata['mapInfo'] = refl['Metadata']['Coordinate_System']['Map_Info'][()]

    #Extract no data value & set no data value to NaN
    metadata['noDataVal'] = float(reflArray.attrs['Data_Ignore_Value'])
    metadata['scaleFactor'] = float(reflArray.attrs['Scale_Factor'])
    
    #Extract bad band windows
    metadata['bad_band_window1'] = (refl.attrs['Band_Window_1_Nanometers'])
    metadata['bad_band_window2'] = (refl.attrs['Band_Window_2_Nanometers'])
    
    #Extract projection information
    metadata['projection'] = refl['Metadata']['Coordinate_System']['Proj4'][()]
    metadata['epsg'] = int(refl['Metadata']['Coordinate_System']['EPSG Code'][()])
    
    #Extract map information: spatial extent & resolution (pixel size)
    mapInfo = refl['Metadata']['Coordinate_System']['Map_Info'][()]
    mapInfo_string = str(mapInfo); 
    mapInfo_split = mapInfo_string.split(",")
    
    #Extract the resolution & convert to floating decimal number
    metadata['res'] = {}
    metadata['res']['pixelWidth'] = float(mapInfo_split[5])
    metadata['res']['pixelHeight'] = float(mapInfo_split[6])
    
    #Extract the upper left-hand corner coordinates from mapInfo
    xMin = float(mapInfo_split[3]) #convert from string to floating point number
    yMax = float(mapInfo_split[4])
    #Calculate the xMax and yMin values from the dimensions
    xMax = xMin + (refl_shape[1]*metadata['res']['pixelWidth']) #xMax = left edge + (# of columns * resolution)",
    yMin = yMax - (refl_shape[0]*metadata['res']['pixelHeight']) #yMin = top edge - (# of rows * resolution)",
    metadata['extent'] = (xMin,xMax,yMin,yMax) #useful format for plotting
    metadata['ext_dict'] = {}
    metadata['ext_dict']['xMin'] = xMin
    metadata['ext_dict']['xMax'] = xMax
    metadata['ext_dict']['yMin'] = yMin
    metadata['ext_dict']['yMax'] = yMax
    hdf5_file.close        
    
    return reflArray, metadata, wavelengths

def array2raster(newRaster,reflBandArray,reflArray_metadata): 
    
    '''array2raster reads in a reflectance array and associated metadata and returns a geotif raster named newRaster.tif
    --------
    Parameters
    --------
        newRaster: string, name of new geotif raster created
        reflBandArray: reflectance array to be converted to raster
        reflArray_metadata: reflectance metadata associated with reflectance array (generated by h5refl2array function)
    --------
    Returns 
    --------
        newRaster.tif: geotif raster created from reflectance array and associated metadata
    --------
    See Also
    --------
    h5refl2array: 
        reads in a NEON hdf5 reflectance file and returns the reflectance array, select metadata, and the wavelength dataset
    Example:
    --------
    array2raster('SERC_b56_clean.tif',SERC_b56_clean,sercRefl_md) ''' 
    
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = reflBandArray.shape[2]
    pixelWidth = float(reflArray_metadata['res']['pixelWidth'])
    pixelHeight = -float(reflArray_metadata['res']['pixelHeight'])
    originX = reflArray_metadata['ext_dict']['xMin']
    originY = reflArray_metadata['ext_dict']['yMax']
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create('hopb_b56.tif', cols, rows, bands, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(bands)
    outband.WriteArray(reflBandArray)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg']) 
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

reflBandArray = reflArray
reflArray_metadata = srerRefl_md

try :
  import gdal
except:
  from osgeo import gdal

def array_to_geotiff(array,hyObj,dstFile):
    """
    Export numpy array as geotiff.

    Parameters
    ----------
    array : Numpy array
    hyObj : HyTools objects corresponding to the input array
    dstFile : Output filename
    
    Returns
    -------
    None
    
    Geotiff saved to dstFile
    
    """
    
    if  hyObj.file_type == "ENVI":
        gdalFile = gdal.Open(hyObj.file_name)
        projection =gdalFile.GetProjection()
        transform = gdalFile.GetGeoTransform()
        
    elif hyObj.file_type == "HDF":
        print("HDF not supported yet.")
        return
    else:
        print("ERROR: File format not recognized.")
        return
    
    datatype_dict = {np.dtype('int16'): gdal.GDT_Int16,
                     np.dtype('int32'): gdal.GDT_Int32,
                     np.dtype('float32'): gdal.GDT_Float32,
                     np.dtype('float64'): gdal.GDT_Float64,
                     }
    
    datatype = datatype_dict[array.dtype]
    
    # Set the output raster transform and projection properties
    driver = gdal.GetDriverByName("GTIFF")
    tiff = driver.Create(dstFile,array.shape[1],array.shape[0],array.shape[2],datatype)
    tiff.SetGeoTransform(transform)
    tiff.SetProjection(projection)
        
    # Write bands to file
    for band in range(array.shape[2]):
        tiff.GetRasterBand(band +1).WriteArray(array[:,:,band])
        tiff.GetRasterBand(band +1).SetNoDataValue(hyObj.no_data)
        
    del tiff, driver
