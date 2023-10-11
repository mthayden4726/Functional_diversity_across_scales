
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

files = ['TEAK_flightlines/20190615_171251_output_0.tif', 
         'TEAK_flightlines/20190615_171251_output_1.tif',
         'TEAK_flightlines/20190615_171251_output_2.tif',
         'TEAK_flightlines/20190615_171251_output_3.tif',
         'TEAK_flightlines/20190615_171251_output_4.tif',
         'TEAK_flightlines/20190615_171251_output_5.tif',
         'TEAK_flightlines/20190615_171251_output_6.tif',
         'TEAK_flightlines/output_fullarray_170625.tif']

# Load the polygon for clipping ()
clip_file = 'TEAK_clip.shp'
s3.download_file(bucket_name, clip_file, Out_Dir + '/clip_polygon.shp')
clip_polygon = Out_Dir + '/clip_polygon.shp'
with fiona.open(clip_polygon, "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]
print(shapes)
for i, file in enumerate(files):
    print('Loading file from S3')
    s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
    flight = Out_Dir + '/file_' + str(i) + '.tif'
    print(flight)
    with rasterio.open(flight) as src:
        try:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta
            print("File Clipped")
            print(out_meta)
            out_meta.update({"driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform})
            local_file_path = Out_Dir + "/clip.tif"
            with rasterio.open(local_file_path, "w", **out_meta) as dest:
                dest.write(out_image)
            print("File Written")
            destination_s3_key = 'TEAK_flightlines/Clipped_file_' + str(i) + '.tif'
            upload_to_s3(bucket_name, local_file_path, destination_s3_key)
            print("File uploaded to S3")
            os.remove(local_file_path)
        except ValueError as e:
            # Handle the case where there is no overlap between the raster and the shapefiles
            print(f"Skipping file {i} as it does not overlap with the shapefile.")
    os.remove(flight)
    print("Cleared data files")
