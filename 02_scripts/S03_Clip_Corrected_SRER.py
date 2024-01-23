
# Script to pull rasters from S3, clip them to the extent for the field sites, and then re-export them to S3.
# Next they will be pulled in to be merged. 
# Meghan Hayden - 10/11/2023

import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import fiona
import glob
import os
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
from osgeo import gdal

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

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

def download_shapefile(bucket, prefix, output_dir):
    # List all files with the given prefix
    files = s3.list_objects(Bucket=bucket, Prefix=prefix)['Contents']

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    downloaded_files = []
    
    # Download all files
    for file in files:
        key = file['Key']
        local_path = os.path.join(output_dir, os.path.basename(key))
        s3.download_file(bucket, key, local_path)
        downloaded_files.append(local_path)

    return downloaded_files

#files = ['TEAK_flightlines/20190615_171251_output_0.tif', 
#         'TEAK_flightlines/20190615_171251_output_1.tif',
#         'TEAK_flightlines/20190615_171251_output_2.tif',
#         'TEAK_flightlines/20190615_171251_output_3.tif',
#         'TEAK_flightlines/20190615_171251_output_4.tif',
#         'TEAK_flightlines/20190615_171251_output_5.tif',
#         'TEAK_flightlines/20190615_171251_output_6.tif',
#         'TEAK_flightlines/output_fullarray_170625.tif']

# Find files for mosaicing (define search terms)
search_criteria1 = "20190901"
#search_criteria2 = "20190901_17"
dirpath = "TEAK_flightlines/"

# List objects in the S3 bucket in the matching directory
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria1 in obj['Key'])]
print(files)
# Or select files based on QGIS identification
#files = ['TEAK_flightlines/20190901_164806_output_.tif']
# List shapefile prefices
shapefiles = ['Site_boundaries/SRER/SRER_023',
              'Site_boundaries/SRER/SRER_002',
              'Site_boundaries/SRER/SRER_027',
              'Site_boundaries/SRER/SRER_026',
              'Site_boundaries/SRER/SRER_021',
              'Site_boundaries/SRER/SRER_006',
              'Site_boundaries/SRER/SRER_003',
              'Site_boundaries/SRER/SRER_028',
              'Site_boundaries/SRER/SRER_014']
               #             'Site_boundaries/SERC/SERC_009_EPSG',
               #             'Site_boundaries/SERC/SERC_005_EPSG',
               #             'Site_boundaries/SERC/SERC_004_EPSG',
               #             'Site_boundaries/SERC/SERC_044_EPSG',
               #             'Site_boundaries/SERC/SERC_012_EPSG',
               #             'Site_boundaries/SERC/SERC_001_EPSG']

# Load the polygon for clipping ()
for j,shape in enumerate(shapefiles):
    print(shape)
    downloaded_files = download_shapefile(bucket_name, shape, Out_Dir)
    shapefile_path = next(file for file in downloaded_files if file.endswith('.shp'))
    
    # Open shapefile and access geometry
    with fiona.open(shapefile_path, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    print(shapes)
    #minx, miny, maxx, maxy = box(*shapes[0].bounds).bounds
    
    # Access the bounds of the entire GeoDataFrame
    #gdf = gpd.read_file(shapefile_path)
    #minx, miny, maxx, maxy = gdf.total_bounds
    
    for i, file in enumerate(files):
        print('Loading file from S3')
        s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
        flight = Out_Dir + '/file_' + str(i) + '.tif'
        print(flight)
        with rasterio.open(flight) as src:
            try:
                out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True, nodata = 0)
                out_meta = src.meta
                print("File Clipped")
                print(out_meta)
                out_meta.update({"driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                    "nodata": 0})
                local_file_path = Out_Dir + "/clip.tif"
                with rasterio.open(local_file_path, "w", **out_meta) as dest:
                    dest.write(out_image)
                    print("File Written")
                destination_s3_key = 'SRER_flightlines/Shape_' + str(j) + '_Clipped_file_' + str(i) + '.tif'
                upload_to_s3(bucket_name, local_file_path, destination_s3_key)
                print("File uploaded to S3")
                os.remove(local_file_path)
            except ValueError as e:
                # Handle the case where there is no overlap between the raster and the shapefiles
                print(f"Skipping file {i} as it does not overlap with the shapefile.")
        os.remove(flight)
        print("Cleared data files")
