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
import numpy as np

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

def fix_no_data_value(input_file, output_file, no_data_value=0):
    with rasterio.open(input_file, "r+") as src:
        src.nodata = no_data_value
        with rasterio.open(output_file, 'w',  **src.profile) as dst:
            for i in range(1, src.count + 1):
                band = src.read(i)
                band = np.where(band==no_data_value,no_data_value,band)
                dst.write(band,i)

# List clipped files in S3
files = ['TEAK_flightlines/Clipped_file_update_0.tif',
         'TEAK_flightlines/Clipped_file_update_1.tif',
         'TEAK_flightlines/Clipped_file_update_2.tif',
         'TEAK_flightlines/Clipped_file_update_3.tif',
         'TEAK_flightlines/Clipped_file_update_4.tif',
         'TEAK_flightlines/Clipped_file_update_5.tif',
         'TEAK_flightlines/Clipped_file_update_6.tif',
         'TEAK_flightlines/Clipped_file_update_7.tif']

src_files_to_mosaic = []
for i,file in enumerate(files):
    print("Loading file from S3")
    #s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
    flight  = Out_Dir + '/file_' + str(i) + '.tif'
    print(flight)
    #flight_update = Out_Dir + '/update_' + str(i) + '.tif'
    #no_data_value = 0
    #fix_no_data_value(flight, flight_update,no_data_value) 
    #with rasterio.open(flight, 'r') as src:
    src = rasterio.open(flight, 'r') 
    src_files_to_mosaic.append(src)
    print(src.nodata)
    print(src.meta)

print(src_files_to_mosaic)

#def pre121_max(old_data, new_data, old_nodata, new_nodata, **kwargs):
#    mask = np.logical_and(~old_nodata, ~new_nodata)
#    old_data[mask] = np.maximum(old_data[mask], new_data[mask])
#    mask = np.logical_and(old_nodata, ~new_nodata)
#    old_data[mask] = new_data[mask]

mosaic, out_trans = rasterio.merge.merge(src_files_to_mosaic, method = "max")
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
destination_s3_key = 'TEAK_flightlines/Mosaic_TEAK.tif'
local_file_path = Out_Dir + '/mosaic.tif'
upload_to_s3(bucket_name, local_file_path, destination_s3_key)
print("File uploaded to S3")
os.remove(local_file_path)


