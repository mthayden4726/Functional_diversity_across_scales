# Documenation of Workflow for BioSCape Scale-Normalized Functional Diversity 

## Table of Contents
1. [Environment Setup for BioSCape: Biodiversity Across Scales](#environment-setup-for-bioscape)
2. [Function Library](#function-library)
3. [Image Correction: Topographic & BRDF Correction of Flightlines](#image-correction:-topographict-&-BRDF-Correction-of-Flightlines)
4. [Calculating Sun Angles](#calculating-sun-angles)
5. [Extracting Slope and Aspect for Drone Data using DEM](#extracting-slope-and-aspect-for-drone-data-using-dem)
6. [Topographic Correction using Methods for Drone Data](#topographic-correction-using-methods-for-drone-data)
7. [Topo and BRDF Correction using Hytools (Steps to Use It)](#topo-and-brdf-correction-using-hytools-steps-to-use-it)
8. [Resampling](#resampling)
9. [NEON Data Access using API](#neon-data-access-using-api)
10. [About Cyverse](#about-cyverse)

## Environment Setup for BioSCape
This section details the steps for setting up the environment required for the BioSCape project. Follow these steps to ensure a smooth and consistent development environment.
### Dependencies
Before starting, make sure you have the [environment.yml](https://github.com/mthayden4726/BioSCape_across_scales/blob/4fd0ba43a21235032893b45b26042cf354a5a22a/environment.yml) file. This file contains all the necessary dependencies.

If you are setting up a new instance, follow the [Step-by-Step Guide](#step-by-step-guide). If you have already set up an instance, skip ahead to [restarting an instance](#restarting-an-instance). 
### Step-by-Step Guide
1.  Open Command Line Interface
Depending on your operating system, use the following:
    * Windows Users: Anaconda Prompt or Command Prompt.
    * macOS/Linux Users: Terminal application.
2.  SSH into AWS instance
    * ``` ssh -i meha3816.pem ec2-user@ec2-44-238-190-104.us-west-2.compute.amazonaws.com ```
    * The path will vary depending on the user and instance ID. To find the instance ID, click on "InstanceID" and then select "Connect" in the AWS EC2 interface.
3.  Install anaconda, run the installer and make sure it is on your path.
    * ``` wget https://repo.anaconda.com/archive/Anaconda3-2023.03-1-Linux-x86_64.sh ```
    * ``` bash Anaconda3-2023.03-1-Linux-x86_64.sh ```
    * ``` source ~/.bashrc ```
4. Install git
    * ``` sudo yum install git ```
5. Clone the git repository
    * ``` git clone https://github.com/mthayden4726/BioSCape_across_scales ```
    * You will be asked for a password for which you can enter your personal token. The personal token for **mthayden4726** expires July 13, 2024 and is: **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**
6. Create and activate the environment
    * First, enter the project directory: ```cd BioSCape_across_scales ```
    * Next, create environment: ``` conda env create -n bioscape-env --file environment.yml ```
          * Note: If this environment does not build, I have additional steps listed in steps #9 & #11 in [AWS & Git Integration Steps](https://docs.google.com/document/d/1slMC_8aWb2bYjJKfHYCFtAxefhWc12xr_H9IauA0B7g/edit)
      * If the conda build is taking too long, use libmama instead. Follow these alternative steps to update conda, install libmama and then build the environment:
        1. ``` conda update conda ```
        2. ``` conda update conda-build ```
        3. ``` conda install -n base conda-libmamba-solver ```
        4. ``` conda config --set solver libmamba ```
        5. ``` conda env create -n bioscape-env --file environment.yml ```
    * Finally, activate the environment: ``` source activate bioscape-env ```
7. After these steps, your environment is set up and ready for project work. For all future work on this instance, start at [Re-starting an instance](#re-starting-an-instance).

### Re-starting an instance
Once your instance is launched on AWS, for all subsequent times you connect you can:
1. Open terminal
2. ssh into instance (e.g., ssh -i meha3816.pem ec2-user@ec2-44-238-190-104.us-west-2.compute.amazonaws.com)
3. Start a virtual screen
    * *This is not strictly necessary but I like to do this to ensure that the code continues to run if my local machine stops*
    * To start a new screen, ``` screen -S name ```
    * To resume an existing screen, ``` screen -r name ```
    * To exit a screen, Ctrl + A.
5. Change directory and activate environment.
    * ``` cd BioSCape_across_scales ```
    * ``` source activate bioscape-env ```

## Function Library
ADD TEXT HERE
## Image Correction: Topographic & BRDF Correction of Flightlines
The first step of the workflow is to implement topographic and BRDF corrections based on correction coefficients provided by Kyle Kovach. *If correction coefficients are not available, this step could be skipped*
This section outlines the steps for correcting the NEON data product, spectrometer orthorectified surface directional reflectance: [DP1.30006.001](https://data.neonscience.org/data-products/DP1.30006.001)). 

**Objective:** Correct variations in reflectance caused by BRDF & topographic effects like slope and aspect using the script [02_scripts/S02_TopoBRDF_Corrections.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/a73d3ea27cc2bff5dd30d4a6140351a09a150007/02_scripts/S02_TopoBRDF_Corrections_BART.py)

### Implementation:
This script requires the following input from users:
   1. Name of the NEON site (e.g., BART)
   2. Domain of the NEON site (e.g., D01)
   3. EPSG of NEON site (e.g., 32619)
   4. Date of desired flights (e.g., 20190825)
Global parameters include:
   * ndvi_threshold = 0.25
   * nir_band = 90
   * red_band = 58
   * Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
   * Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
   * bucket_name = 'bioscape.gra'
   * s3 = boto3.client('s3')

The script includes:
*   Find all flightlines with matching correction coefficients (located in S3 bucket "NEON BRDF-TOPO Corrections/")
*   Download those flightlines. For each flightline:
   *   Load matching BRDF & Topo correction coefficients.
   *   Remove bad bands.
   *   Apply corrections.
   *   Calculate NDVI and apply NDVI threshold.
   *   Rasterize the corrected array, update metadata and export as a .tif.
   *   Upload to S3 and delete file locally. 

For a detailed walkthrough, see the notebook: [NA]


### SCS+C (Sun-Canopy-Sensor + Cosine) Topographic Correction
**Objective:** Extend the SCS method by including a cosine correction factor, enhancing effectiveness in rugged terrains.

**Implementation:**

The notebook includes:

*   Parameters are derived from metadata within the NEON file, with solar zenith angles averaged over multiple flightlines.

*   Extract Parameters from Metadata: Similar to Notebook 1, detailing the extraction process of relevant parameters from the NEON data metadata.

*   Topographic Correction using SCS+C Method: Implementing the SCS+C method for topographic correction which includes an additional cosine correction factor.

*   Various Graphical Plots for Topographic Correction: Visualizing data and correction effects, especially focusing on NIR band 93.

*   Function for Plotting Aspect and Illumination: Similar to Notebook 1, but tailored for the SCS+C method.

*   Statistical Analysis of Pixel Values: Examining the impact of SCS+C correction on pixel values.

*   Correlation Analysis: Assessing correlations post-topographic correction using the SCS+C method.

*   NDVI Analysis and Correction: Detailed steps for applying NDVI on corrected data, including methods for choosing the nearest red and NIR bands.

*   Comparative NDVI Graphs: Visual comparison between reflectance NDVI and post-correction NDVI.

For an in-depth explanation and code, view: [Topo_Corr_SCS_C.ipynb](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Topo_Corr_using_Methods/Topo_Corr_SCS_C_Final.ipynb)

## Calculating Sun angles

This section of the documentation outlines the methodologies used in the Jupyter notebook for calculating the sun angles, which is a crucial step in topographic and radiometric corrections for remote sensing data.

### Overview
The notebook provides two distinct methods for calculating sun angles, leveraging geographical data from drone imagery. It includes functions to extract latitude and longitude coordinates from the drone images and to compute sun angles for each pixel.

### Detailed Steps
1.  Extracting Latitude and Longitude from Drone Images:
Functionality is developed to extract geographic coordinates (latitude and longitude) directly from the metadata of drone imagery.

2.  Method 1: Basic Sun Angle Calculation
This method calculates the sun angles using a simplified approach, ideal for scenarios where detailed topographic information is not critical.

3.  Method 2: Advanced Sun Angle Calculation
An advanced technique for calculating sun angles, providing more accuracy and detail. This method is particularly useful for rigorous topographic and radiometric analyses.

4.  Function to Calculate Sun Angles for Each Image Pixel:
A comprehensive function that calculates the sun angles for every pixel in the drone image, using the extracted latitude and longitude coordinates. This function is essential for detailed pixel-by-pixel analysis in remote sensing applications.

### Application
These methods are crucial for understanding the solar illumination conditions for each pixel, which significantly impacts the reflectance values in remote sensing data. Accurate sun angle calculation allows for more precise corrections and analyses in subsequent steps of the project.

For an in-depth explanation and code, view: [calculating_sun_angles.ipynb](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Topo_Corr_using_Methods/Calulating_Sun_angles.ipynb)



## Extracting Slope and Aspect for drone data using DEM

This section of the documentation details the process outlined in the Jupyter notebook for deriving slope and aspect data from Digital Elevation Models (DEM) applied to drone imagery. This procedure is crucial for understanding the topography of the area captured by drone data and getting slope and aspect for further corrections as they are the required parameters.

### Overview
The notebook demonstrates the procedure to calculate slope and aspect using DEM data, which are critical parameters in topographical analysis and can significantly impact the accuracy of remote sensing data interpretation.

### Detailed Steps
*   Loading the Digital Elevation Model (DEM) File: The first step involves loading the DEM file specified in the DEM_path variable. This file contains the elevation data necessary for calculating slope and aspect.
*   Calculating Slope and Aspect:The notebook provides code for calculating two key topographic parameters:
    *   Slope: This measures the steepness or degree of incline of the terrain. The slope is essential for understanding the terrain's gradient and is calculated in degrees.
    *   Aspect: This refers to the compass direction that the slope faces. Aspect is crucial for determining the direction of the sun's illumination on the terrain.
    
*   Output Path for Saving Results:
The notebook includes functionality to save the calculated slope and aspect data to specified output paths. This ensures that the results are stored for further analysis and use in the project.

For an in-depth explanation and code, view: [slope_aspect_drone_dtm.ipynb](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Topo_Corr_using_Methods/Slope_Aspect_Drone_Data.ipynb)

### Application
Slope and aspect data extracted from DEM are fundamental in environmental and geographical studies, especially in projects involving remote sensing and aerial imagery. They provide insights into the terrain characteristics, which are vital for various analyses, including ecological studies, land-use planning, and agricultural assessments.

## Topographic Correction using methods for Drone Data

### Overview
The notebook focuses on applying the SCS+C method for correcting topographic effects in drone imagery. This method is particularly effective for landscapes with varied terrain, as it accounts for differences in sunlight angles and terrain features.

### Detailed Steps
1.  Assigning Paths for Drone Data, Slope, and Aspect:
Define file paths for the drone data, slope, and aspect datasets necessary for topographic correction.

2.  Setting Parameters for SCS+C Topographic Correction:
Outlines the parameters needed for the SCS+C method, including the reading of specific bands from drone data.

3.  Function for Slope and Aspect Extraction:
Describes the process of using rasterio to extract slope and aspect data, essential components for the SCS+C correction.

4.  Calculating Sun Angles:
Details the calculation of sun angles for each pixel in the drone imagery, a key factor in the topographic correction process.

5.  Gathering Parameters for Topographic Correction:
Consolidates all necessary parameters (including illumination, slope, aspect, and sun angles) for executing the SCS+C correction.

6.  Illumination Plotting:
Visualizing the illumination component, crucial for understanding the light dynamics over the terrain.

7.  Topographic Correction using SCS+C Method:
Implementation of the SCS+C topographic correction, including a detailed function to perform this correction on the drone data.

8.  Visual Analysis for NIR Band:
Includes plots and analysis focusing on the NIR band to assess the effectiveness of the topographic correction.
9.  NDVI Calculation:
Demonstrates the calculation of the Normalized Difference Vegetation Index (NDVI), both with and without a threshold, to analyze vegetation health post-correction.

For an in-depth explanation and code, view: [Topo_Corr_drone_data.ipynb](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Topo_Corr_using_Methods/Topo_corr_Drone_data.ipynb)

#### Application
The application of topographic correction using the SCS+C method is vital in ensuring the accuracy of remote sensing data, especially in areas with complex terrain. This notebook provides a comprehensive guide for applying this correction to drone imagery, which can be crucial for environmental monitoring, agricultural assessments, and land-use studies.
## Topo and BRDF correction using Hytools (steps to use it)

### Overview
This documentation covers a comprehensive workflow for processing NEON data using Python scripts. The process involves converting NEON data to ENVI format, generating configuration files for topographic and BRDF corrections, and applying these corrections to the imagery.

### Python Scripts Description
1. neon2envi2.py: NEON to ENVI Conversion, code: [neon2envi.py](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/BRDF-Topo-HyTools/Topo%20and%20Brdf%20Corr/neon2envi2.py)
Purpose: Converts NEON AOP H5 data files to ENVI format.
Usage:
Run the script from the command line with the dataset path and output folder.
Optional flag -anc to export ancillary data.
Example:
```python neon2envi2.py <path-to-dataset_name> <path-to-output_folder> -anc```
2. config_generator.py: Configuration File Generation, code: [config_generator.py](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/BRDF-Topo-HyTools/Topo%20and%20Brdf%20Corr/config_generator.py)
Functionality: Generates JSON configuration files for applying topographic (TOPO) and Bidirectional Reflectance Distribution Function (BRDF) corrections.
Configuration Options: Includes settings for various correction types, wavelengths, and other parameters.
Running the Script: Edit the script according to the desired corrections and run it to create config_<iteration>.json files.
Example Command:
```python config_generator.py```
3. image_correct.py: Applying Corrections, code: [image_correct.py](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/BRDF-Topo-HyTools/Topo%20and%20Brdf%20Corr/image_correct.py)
Purpose: Reads the generated JSON configuration file and applies the specified TOPO and BRDF corrections to the imagery.
Execution: Run the script with the configuration file as a command-line argument.
```python image_correct.py <path-to-config-file>```


### Overview for config_generator.py
The config_generator.py script is designed to automate the creation of configuration files for topographic (TOPO) and Bidirectional Reflectance Distribution Function (BRDF) corrections of geospatial imagery. It allows customization to accommodate different correction methods and input formats.

#### TOPO Correction Methods
The script supports various TOPO correction methods, including SCS (Sun-Canopy-Sensor), SCS+C (Sun-Canopy-Sensor + Cosine), and C correction. For this project, the SCS+C method has been chosen due to its effectiveness in handling varied terrain by incorporating an additional cosine correction factor.

##### Key Features of SCS+C Method:
*   Accounts for solar zenith angle, slope, and aspect.
*   Adjusts reflectance values based on pixel-specific illumination conditions.
*   Particularly effective in landscapes with significant elevation changes.

#### BRDF Correction Methods
Two primary methods for BRDF correction are supported: the Universal method and the Flex method. In this project, the Flex method is used due to its adaptability and suitability for the specific requirements of NEON data.

##### Key Features of Flex Method:
*   Tailors BRDF corrections based on scene-specific characteristics.
*   Handles a wide range of surface and atmospheric conditions.
*   Note: Diagnostic plots are more challenging with the Flex method compared to the Universal method, as Flex returns different values that require extensive modifications to the HyTools library.

##### Customization and Preferences
Users can modify the config_generator.py script to choose their preferred methods for both TOPO and BRDF corrections. The script is structured to allow easy switching between different correction algorithms and settings.

#### Configuration File Generation:
*   The script generates JSON files for each ancillary file related to the main image.
*   Users can specify bad bands, file types, input files, and other settings.
*   The output includes detailed settings for export options, masks, coefficients, and correction-specific parameters.
#### Usage:
*   Ideal for workflows requiring specific topographic and BRDF corrections.
*   Users can edit the script to select desired correction methods and parameters.
*   The output JSON files serve as input for subsequent correction processes using tools like image_correct.py.

#### Flexibility and Extensibility:
The configuration generator script offers flexibility and extensibility, allowing users to adapt the correction process to their specific needs. By modifying the script, users can experiment with different correction methods and parameters, optimizing their workflow for the best possible results in geospatial imagery analysis.
### Steps to Run the Workflow ([readme](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/BRDF-Topo-HyTools/Topo%20and%20Brdf%20Corr/README.md))
*   Step 1: Convert NEON Data to ENVI
Ensure the output folder exists before running the conversion script.
Example:
```python neon2envi2.py neon.h5 output/ -anc```
*   Step 2: Generate Configuration JSON
Modify config_generator.py as needed for the specific corrections.
Run the script to generate the configuration file.
Example:
```python config_generator.py```
*   Step 3: Perform Correction
Use image_correct.py with the generated config file to apply corrections.
Example:
```python image_correct.py output/config_01.json```

### Applications of the Workflow
This workflow is ideal for remote sensing professionals and researchers working with NEON data who require precise spectral matching and corrections for their analysis. The streamlined process from data conversion to correction application ensures accuracy and efficiency in multispectral and ecological studies.


## Resampling
### Overview 
This tool facilitates the resampling of NEON and drone hyperspectral data to align with Landsat sensor specifications. It utilizes a JSON file for defining sensor parameters and a Python script incorporating the HyTools library for the resampling process.
### [Readme](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Resampling/README_RESAMPLING.md)
### Components
*   [landsat_band_parameters.json](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Resampling/landsat_band_parameters.json): This JSON file includes band parameters for various Landsat missions. It can be adapted to include specifications for NEON and drone sensors, allowing users to resample their data to match the Landsat spectral response.

*   [resampling_demo.py](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/Resampling/resampling_demo.py): A Python script that executes the resampling of NEON and drone data. Key components include:

    *   resampler_hy_obj Class: Initializes with the sensor type (Landsat, NEON, or drone), JSON file path, and HDR file path. Manages the reading and application of band parameters for resampling.

    *  create_header_info Method: Generates necessary header information from the HDR file for the resampling process.

    *   save_envi_data Method: Saves resampled data in the ENVI format.

    *   load_envi_data Function: Loads hyperspectral data from an ENVI binary file, preparing it for resampling.

### Prerequisites
*   Python 3.x
*   NumPy
*   Spectral Python (SPy)
*   HyTools library

### Installation

Install necessary Python libraries:
```pip install numpy spectral hytools```

### Usage
*   Update the landsat_band_parameters.json with NEON or drone sensor parameters if required.
*   Ensure both the JSON file and resampling_demo.py are in your working directory.
*   Execute the script with appropriate arguments for sensor type, file paths, and output specifications.

#### Example
```python
from resampling_deom import resampler_hy_obj

# Initialize resampler object for Landsat 8 OLI
resampler = resampler_hy_obj(sensor_type='Landsat 8 OLI', json_file='landsat_band_parameters.json')

# Apply resampling (add details based on your data and requirements)
```

### Application
This tool is crucial for researchers and professionals working with NEON and drone hyperspectral data who need to align their datasets with Landsat spectral characteristics. Such resampling is essential for comparative analysis across different sensors, enhancing the validity of environmental and geographical studies.


## NEON data access using API
### Overview
This guide provides a generalized method to access data from the National Ecological Observatory Network (NEON) by altering the site name and product ID. It is based on a Python script that utilizes the NEON API for data retrieval.
### [Code](https://github.com/earthlab/cross-sensor-cal/blob/janushi-main/neon-api.ipynb)
### Step-by-Step Guide
1. Defining NEON API Endpoint
Define the base URL for the NEON API. For example:
```NEON_API_ENDPOINT = "https://data.neonscience.org/api/v0/"```
2. Specifying Site and Product ID
Assign variables for the site name and product ID. These can be changed to access different datasets.
```site_name = "NEON_SITE_NAME"  # Replace with desired site name```
```product_id = "NEON_PRODUCT_ID"  # Replace with specific product ID```
3. Constructing the API Request
Build the request URL using the specified site and product ID, and make the API request:
```request_url = f"{NEON_API_ENDPOINT}data/{product_id}/{site_name}"```
```response = requests.get(request_url)```
4. Parsing the Response
Convert the response to a JSON format and extract relevant data:
```data = response.json()```
4. Handling Data
Depending on your requirements, process or analyze the retrieved data. This might involve data cleaning, analysis, visualization, etc.


### Application
This method is ideal for ecologists, environmental scientists, and data analysts who require access to NEON's vast ecological datasets. By simply changing the site name and product ID, a wide range of ecological data can be accessed and utilized for research and analysis.


## About Cyverse

*   I have downloaded NIWO and RMNP data from year 2020 and stored them inside commuitydata -> earthlab -> macrosystems -> NIWO_and_RMNP
*   I have also cloned gitHub repo into macrosystems environment
*   While following above steps we can run the files required.

