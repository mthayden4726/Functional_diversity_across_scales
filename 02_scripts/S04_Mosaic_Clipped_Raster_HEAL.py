import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.warp import transform_bounds
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject
import glob
import os
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
import fiona
from fiona.crs import from_epsg
import pycrs
from osgeo import gdal
import numpy as np
from S01_Functions import * # add scripts folder to python path manager

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# Find files for mosaicing (define search terms)
#search_criteria1 = "20190515"
#search_criteria2 = "20190901_17"
#dirpath = "SERC_flightlines/"

# List objects in the S3 bucket in the matching directory
#objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
#files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria1 in obj['Key'])]
#print(files)
# Or select files based on QGIS identification
#files = ['TEAK_flightlines/20190901_164806_output_.tif']
# List shapefile prefices
#shapefiles = ['Site_boundaries/SERC/SERC_010_EPSG',
#              'Site_boundaries/SERC/SERC_009_EPSG',
#              'Site_boundaries/SERC/SERC_005_EPSG',
#              'Site_boundaries/SERC/SERC_004_EPSG',
#              'Site_boundaries/SERC/SERC_044_EPSG',
#              'Site_boundaries/SERC/SERC_012_EPSG',
#              'Site_boundaries/SERC/SERC_001_EPSG']

# Set config options
gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
gdal.SetConfigOption('CHECK_DISK_FREE_SPACE', 'FALSE')

#for i, file in enumerate(files):
#    print('Loading file from S3')
#    s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
#    flight = Out_Dir + '/file_' + str(i) + '.tif'
#    print(flight)
#    with rasterio.open(flight) as src:
#        clipped_data, transform = rasterio.mask.mask(src, clip_polygon.geometry, crop=True)
#        clipped_dataset = rasterio.open('clipped_file_' + str(i) + '.tif', 'w', **src.profile)
#        print(clipped_dataset)
#        clipped_dataset.write(clipped_data)
#        datasets.append(clipped_dataset)
#print(datasets)
# Open in read mode and add to file list

src_files_to_mosaic = []
#for j,shape in enumerate(shapefiles):
#    
#    print(shape)
#
#    # Download shapefile files
#    downloaded_files = download_shapefile(bucket_name, shape, Out_Dir)
#    shapefile_path = next(file for file in downloaded_files if file.endswith('.shp'))
#    
#    # Open shapefile and access geometry
#    with fiona.open(shapefile_path, "r") as shapefile:
#        shapes = [feature["geometry"] for feature in shapefile]
#    
#    #minx, miny, maxx, maxy = box(*shapes[0].bounds).bounds
#    
#    # Access the bounds of the entire GeoDataFrame
#    gdf = gpd.read_file(shapefile_path)
#    #minx, miny, maxx, maxy = gdf.total_bounds
#    #print(minx, miny, maxx, maxy)
#
file_ID = ['002',
           '004',
           '005',
           '013',
           '015',
           '018',
           '019',
           '024',
           '026']

for i,ID in enumerate(file_ID):
    
    # List files associated with a single buffer shape
    search_criteria = str(ID) + '_Clipped_file_'
    dirpath = "HEAL_flightlines/Site_boundaries/HEAL_"

    # List objects in the S3 bucket in the matching directory
    objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
    # Filter objects based on the search criteria
    files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria in obj['Key'])]
    print(files)
    for j,file in enumerate(files):
        flight  = Out_Dir + '/file_' + str(j) + '.tif'
        try:
            s3.download_file(bucket_name, file, flight)
            print(f"The file '{file}' exists.")
        except Exception as e:
            print(f"Error: {e}")
        src = rasterio.open(flight)
        src_files_to_mosaic.append(src)

    # Mosaic files
    print(src_files_to_mosaic)
    mosaic, out_trans = merge(src_files_to_mosaic, nodata = -9999)
    print('Merge complete')
    # Update metadata
    out_meta = src.meta.copy()
    print(out_meta)
    print(mosaic.shape)
    out_meta.update({"driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": "+init=epsg:32618 +units=m +no_defs "})
    print(out_meta)

    # Write to computer, send to S3
    mosaic_name = Out_Dir + "/mosaic_SRER.tif"
    with rasterio.open(mosaic_name, "w", **out_meta) as dest:
        dest.write(mosaic)
    print("File written")
    
    # Push to S3 bucket
    destination_s3_key = 'HEAL_flightlines/Mosaic_HEAL_'+str(ID)+'.tif'
    local_file_path = mosaic_name
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    print("File uploaded to S3")
    
    # Remove unneeded files (mosaic and shapefile)
    os.remove(local_file_path)
    #os.remove(shape)

