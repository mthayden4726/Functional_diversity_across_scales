import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
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

# List clipped files in S3
files = ['TEAK_flightlines/Clipped_file_0.tif',
         'TEAK_flightlines/Clipped_file_1.tif',
         'TEAK_flightlines/Clipped_file_2.tif',
         'TEAK_flightlines/Clipped_file_3.tif',
         'TEAK_flightlines/Clipped_file_4.tif',
         'TEAK_flightlines/Clipped_file_5.tif',
         'TEAK_flightlines/Clipped_file_6.tif',
         'TEAK_flightlines/Clipped_file_7.tif']

src_files_to_mosaic = []
for i,file in enumerate(files):
    print("Loading file from S3")
    s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
    flight  = Out_Dir + '/file_' + str(i) + '.tif'
    print(flight)
    src = rasterio.open(flight)
    src_files_to_mosaic.append(src)

print(src_files_to_mosaic)

mosaic, out_trans = rasterio.merge.merge(src_files_to_mosaic)
print('Merge complete')

out_meta = src.meta.copy()
print(out_meta)
print(mosaic.shape)
out_meta.update({"driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_trans})
print(out_meta)
# Write to computer, send to S3
local_file_path = Out_Dir + "/mosaic.tif"
with rasterio.open(local_file_path, "w", **out_meta) as dest:
    dest.write(mosaic)
    print("File written")

# Push to S3 bucket
destination_s3_key = 'TEAK_flightlines/Mosaic_0-3.tif'
local_file_path = Out_Dir + '/mosaic.tif'
upload_to_s3(bucket_name, local_file_path, destination_s3_key)
print("File uploaded to S3")
os.remove(local_file_path)


