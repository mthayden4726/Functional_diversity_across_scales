#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:28:57 2023

@author: meha3816
"""
Data_Dir = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/02_processed'

# Select files to pull from NEON
SITECODE = 'SRER' # NEON site of interest
PRODUCTCODE = 'DP1.30006.001' # NEON product of interest (DP3.30006.001 is orthorectified mosaic)
YEAR = '2019-09'
file_paths = find_neon_files(SITECODE,
                             PRODUCTCODE,
                          YEAR)
# Pull from NEON
flights_SRER = retrieve_neon_files(file_paths, Data_Dir)

# Read in image
neon_image= '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON_D14_SRER_DP1_20190901_165553_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')
# Load correction coefficients
topo_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_SRER/NEON_D14_SRER_DP1_20190901_165553_reflectance_topo_coeffs_topo.json'
neon.load_coeffs(topo_coeffs,'topo')
brdf_coeffs = '/Users/meha3816/Library/CloudStorage/OneDrive-UCB-O365/Desktop/BioSCape_across_scales/01_data/01_rawdata/NEON BRDF-TOPO Corrections/2019_SRER/NEON_D14_SRER_DP1_20190901_165553_reflectance_brdf_coeffs_topo_brdf.json'
neon.load_coeffs(brdf_coeffs,'brdf')


X = []
bad_bands = [[300,400],[1300,1450],[1780,2000],[2450,2600]]
neon.create_bad_bands(bad_bands)
ndvi = neon.ndi()
#neon.mask['veg'] = (neon.mask['no_data']) & (ndvi > ndvi_threshold)
X_chunk2 = neon.get_chunk(0, neon.columns, 0, neon.lines, corrections= ['topo','brdf'])
reflBandArray = X_chunk

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
header_dict = neon.get_header()
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = 
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
# Name of output raster
newRaster = 'output.tif'

# Array to raster
outRas = array2raster(newRaster, array, refl_md, Out_Dir, epsg)
def array2raster(newRaster, reflBandArray, reflArray_metadata, extent, ras_dir, epsg):
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
    bands = reflBandArray.shape[2]
    pixelWidth = float(refl_md['res']['pixelWidth'])
    pixelHeight = -float(refl_md['res']['pixelHeight'])
    originX = refl_md['ext_dict']['xMin']
    originY = refl_md['ext_dict']['yMax']
#
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[array.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        outRaster.GetRasterBand(band + 1).WriteArray(array[:, :, band])
#
    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)
    

# Create dictionary containing metadata for creation of raster
# data from NEON
mapInfo= neon.map_info
# empty dictionary
refl_md = {}
refl_md['mapInfo'] = header_dict['map info']
refl_md['wavelength'] = header_dict['wavelength']
refl_md['shape'] = array.shape
#Extract no data value & scale factor
refl_md['noDataVal'] = float(header_dict['data ignore value'])
refl_md['scaleFactor'] = float(0.9996)
refl_md['bad_band_window1'] = np.array([1340, 1445])
refl_md['bad_band_window2'] = np.array([1790, 1955])
refl_md['epsg'] = 4326
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

