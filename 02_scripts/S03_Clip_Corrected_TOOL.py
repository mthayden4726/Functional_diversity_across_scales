
# Script to pull rasters from S3, clip them to the extent for the field sites, and then re-export them to S3.
# Next they will be pulled in to be merged. 
# Meghan Hayden - 2/22/2024

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
from S01_Functions import *

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# Find files for mosaicing (define search terms)
search_criteria = "201907"
dirpath = "TOOL_flightlines/"

# List objects in the S3 bucket in the matching directory
objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
# Filter objects based on the search criteria
files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria in obj['Key'])]
print(files)

# List shapefile prefices
shapefiles = ['Site_boundaries/TOOL/TOOL_020',
              'Site_boundaries/TOOL/TOOL_018',
              'Site_boundaries/TOOL/TOOL_003',
              'Site_boundaries/TOOL/TOOL_026',
              'Site_boundaries/TOOL/TOOL_014',
              'Site_boundaries/TOOL/TOOL_010',
              'Site_boundaries/TOOL/TOOL_024',
              'Site_boundaries/TOOL/TOOL_071',
              'Site_boundaries/TOOL/TOOL_028',
              'Site_boundaries/TOOL/TOOL_043',
              'Site_boundaries/TOOL/TOOL_022',
              'Site_boundaries/TOOL/TOOL_023'
             ]

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
                out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True) # remove nodata = 0 to see if this helps with mosaic
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
                destination_s3_key = 'TOOL_flightlines/' + str(shape) + '_Clipped_file_' + str(i) + '.tif'
                upload_to_s3(bucket_name, local_file_path, destination_s3_key)
                print("File uploaded to S3")
                os.remove(local_file_path)
            except ValueError as e:
                # Handle the case where there is no overlap between the raster and the shapefiles
                print(f"Skipping file {i} as it does not overlap with the shapefile.")
            os.remove(flight)
        print("Cleared data files")
