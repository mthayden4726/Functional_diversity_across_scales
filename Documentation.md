# Documenation of Workflow for BioSCape Scale-Normalized Functional Diversity 

## Table of Contents
1. [Workflow for BioSCape: Biodiversity Across Scales](#workflow-for-bioscape)
2. [Environment Setup for BioSCape: Biodiversity Across Scales](#environment-setup-for-bioscape)
3. [Function Library](#function-library)
4. [Image Correction: Topo/BRDF Correction and Radiometric Filtering of Flightlines](#image-correction)
5. [Image Clipping: Subsetting Data to NEON Plots](#image-clipping)
6. [Image Mosaic: Setting the Scene for Diversity Calculations](#image-mosaic)
7. [Functional Diversity Computation: Calculating Functional Richness and Functional Divergence across Window Sizes](#functional-diversity-computation)
8. [Environmental Covariates: Assessing Environmental Characteristics of NEON Sites](#environmental-covariates)
9. [Scaling Analysis: Creating a Scale-Normalized Functional Diversity Metric](#scaling-analysis)
   
## Workflow for BioSCape 
ADD TEXT

![image](https://github.com/mthayden4726/BioSCape_across_scales/assets/70178120/d83d8a23-c4c2-4a13-86ed-6b104902836b)


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
All of the supporting functions necessary to run this workflow are in the script [02_scripts/S01_Functions.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/608d9a2d1b883640afa6566f39fd39991e04b9ee/02_scripts/S01_Functions.py).

The functions to compute functional richness and functional divergence across neighborhood/window sizes have their own files:
   * [02_scripts/S01_Moving_Window_FRic.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/608d9a2d1b883640afa6566f39fd39991e04b9ee/02_scripts/S01_Moving_Window_FRIC.py)
   * [02_scripts/S01_Moving_Window_FDiv.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/608d9a2d1b883640afa6566f39fd39991e04b9ee/02_scripts/S01_Moving_Window_FDiv.py)

Supporting functions included in this script include functions for:
   * **Interacting with the NEON database**
      * ```find_neon_files()```: returns the file paths for NEON files according to site code, product code, and desired dates
      * ```retrieve_neon_files()```: downloads NEON files from list of file paths
   * **Interacting with AWS S3 Buckets**
      * ```download_shapefile()```: downloads the set of files constituting a shapefile from the specified bucket
      * ```upload_to_s3()```: uploads local file to S3
   * **Processing images**
      * ```arraytoraster()```: converts an array to a single band raster
      * ```array2rastermb()```: converts an array to a multi-band raster
      * ```clip_raster()```: clips a raster to bounds set by a shapefile
      * ```scale_transform()```: centers, scales, and fits a PCA to a raster
      * ```pca_steps()```: takes a raster through all the steps to a PCA (center, scale, fit and transform)
      * ```store_metadata()```: stores metadata of a raster
   * **Visualizing images**
      * ```rbg_show()```: plots hytools object in RGB

### Computing Functional Richness
Our computation of functional richness uses the ConvexHull() function implemented in scipy.spatial. We also wrote a function which subsets the input PCA into subarrays and calculates functional richness for each window using a moving window approach. This gives us a set of functional richness values across every window in the scene, from which we take the median to represent the functional richness for the plot at a given area. 

Specifically, the function ```window_calcs_fric()``` uses a moving window approach to loop across the input PCA, segmenting the PCA into subarrays of the appropriate window size and calculating functional richness as the volume of the convex hull for each subarray.
   * This function is parallelized such that the same approach is running simultaneously for multiple window sizes. 

### Computing Functional Divergence
There are two primary components to compute functional divergence:
   * Function for calculating functional divergence for a single subarray
   * Function for subsetting PCA into subarrays and calculating functional divergence for each using a moving window approach

First, the function ```calc_fdiv()``` computes functional divergence.
   * This function requires a PCA-transformed raster as an input.
   * The function calculates the mean pixel value in dimension 1, and then calculates the euclidean distance of each pixel from the mean.
   * Ultimately, functional divergence combines the mean of those distances and the deviations from the mean, following Villeger et al., 2008 as:
     
     ![image](https://github.com/mthayden4726/BioSCape_across_scales/assets/70178120/b51a62ec-0c03-4593-a95e-83ac7769e3e7)

Second, the function ```window_calcs_fdiv()``` uses a moving window approach to loop across the PCA, segmenting the PCA into subarrays of the appropriate window size and calculating functional divergence for each subarray.
   * This function is parallelized such that the same approach is running simultaneously for multiple window sizes. 


## Image Correction
The first step of the workflow is to implement topographic and BRDF corrections (as well as an NDVI threshold) based on correction coefficients provided by Kyle Kovach. *If correction coefficients are not available, this step could be skipped*
This section outlines the steps for correcting the NEON data product, spectrometer orthorectified surface directional reflectance: [DP1.30006.001](https://data.neonscience.org/data-products/DP1.30006.001)). 

**Objective:** Correct variations in reflectance caused by BRDF & topographic effects like slope and aspect, as well as remove non-vegetated pixels, using the script [02_scripts/S02_TopoBRDF_Corrections.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/a73d3ea27cc2bff5dd30d4a6140351a09a150007/02_scripts/S02_TopoBRDF_Corrections_BART.py)

### Implementation:

This script requires the following input from users:
   1. Name of the NEON site (e.g., BART)
   2. Domain of the NEON site (e.g., D01)
   3. ID of desired flights (e.g., 5)
   4. Date of desired flights (e.g., 20190825)
   5. Date and ID of desired flights (e.g., 2019082513)
   6. EPSG of NEON site (e.g., 32619)
   7. Desired NDVI threshold (e.g., 0.25)

Example command: ``` python 02_scripts/S02_TopoBRDF_Corrections.py --SITECODE BART --DOMAIN D01 --ID_NO 5 --DATE 20190825 --DATE_ID 2019082513 --EPSG 36219 --NDVI 0.25```

*For a list of parameters for all sites included in analysis, see the [2019 NEON Site List](https://docs.google.com/spreadsheets/d/17DJtV1BKtq0uLfcYCtM2kp7JjjjaWpuEV_jt95l_830/edit#gid=124418455)*

Global parameters include:
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

### Application
ADD TEXT

## Image Clipping
The second step of the workflow is to clip corrected flightlines to the desired boundaries. In this case, we want a 2.5 x 2.5 km2 box around each NEON plot of interest. For our implementation, there are two inputs:
1. The corrected flightlines (located in S3 bucket "SITENAME_flightlines/")
2. The shapefiles for each plot boundary (located in S3 bucket "Site_boundaries/")

**Objective:** Clip flightlines to areas of interest (to reduce processing time) using the script [02_scripts/S03_Clip_Corrected.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/4bee04f61cb48478cda5c2c340aaf6c690e3a29c/02_scripts/S03_Clip_Corrected_BART.py)

### Implementation:

This script requires the following input from users:
   1. Name of the NEON site (e.g., BART)
   2. Date of desired files (e.g., 201908). Some sites had flights occur at multiple times throughout the year - for mosaicking purposes, we want to combine flights that were part of the same survey (or occurred at around the same time). 

Example command: ``` python 02_scripts/S03_Clip_Corrected.py --SITECODE BART --YEAR 201908```

*For a list of all sites included in analysis, see the [2019 NEON Site List](https://docs.google.com/spreadsheets/d/17DJtV1BKtq0uLfcYCtM2kp7JjjjaWpuEV_jt95l_830/edit#gid=124418455)*

Global parameters include:
   * Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
   * Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
   * bucket_name = 'bioscape.gra'
   * s3 = boto3.client('s3')

The script includes:
*   Find all corrected flightlines (located in S3 bucket "SITENAME_flightlines/")
*   Load all shapefiles associated with site. For each shapefile:
   * Load a flightline.
   * If there is overlap, clip flightline to shapefile.
   * Update flightline metadata and upload to S3.
   * Repeat with all flightlines.

For a detailed walkthrough, see the notebook: [NA]

## Image Mosaic
The third step of the workflow is to mosaic the clipped and corrected flightlines to create a 2.5 km x 2.5 km scene for analyses. For our implementation, there is one input:
1. The clipped and corrected flightlines (located in S3 bucket "SITENAME_flightlines/Site_boundaries/SITENAME/")

**Objective:** Mosaic clipped and corrected flightlines into scenes for analyses [02_scripts/S04_Mosaic_Clipped_Raster.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/4bee04f61cb48478cda5c2c340aaf6c690e3a29c/02_scripts/S04_Mosaic_Clipped_Raster_BART.py)

### Implementation:

This script requires the following input from users:
   1. Name of the NEON site (e.g., BART)
   2. EPSG of the NEON site (e.g., 32619)

Example command: ``` python 02_scripts/S04_Mosaic_Clipped_Raster.py --SITECODE BART --EPSG 36219```

*For a list of the EPSG of all sites included in analysis, see the [2019 NEON Site List](https://docs.google.com/spreadsheets/d/17DJtV1BKtq0uLfcYCtM2kp7JjjjaWpuEV_jt95l_830/edit#gid=124418455)*

Global parameters include:
   * Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
   * Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
   * bucket_name = 'bioscape.gra'
   * s3 = boto3.client('s3')

The script includes:
*   Find all clipped and corrected flightlines (located in S3 bucket "SITENAME_flightlines/Site_boundaries/SITENAME/")
*   Load all rasters associated with a plot, site. For each file:
   * Load file.
   * Change nodata values to 0.
   * Append to list of files for each plot, site to merge.
*   Merge files using method = "max".
*   Update metadata for mosaicked file.
*   Upload mosaic to S3 and delete local file. 

For a detailed walkthrough, see the notebook: [NA]

## Functional Diversity Computation
The fourth step of the workflow is to compute functional richness and divergence across a set of window sizes for each mosaic produced in the previous step (representing a plot nested within a NEON site). For our implementation, there is one input:
1. The mosaics (located in S3 bucket "SITENAME_flightlines/")

**Objective:** Calculate functional richness and divergence across a set of neighborhood/window sizes for each scene using [02_scripts/S05_Compute_FRic_FDiv.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/4bee04f61cb48478cda5c2c340aaf6c690e3a29c/02_scripts/S05_Compute_FRic_FDiv.py)

### Implementation:

This script requires the following input from users:
   1. Name of the NEON site (e.g., BART)

Example command: ``` python 02_scripts/S05_Compute_FRic_FDiv.py --SITECODE BART```

*For a list of all sites included in analysis, see the [2019 NEON Site List](https://docs.google.com/spreadsheets/d/17DJtV1BKtq0uLfcYCtM2kp7JjjjaWpuEV_jt95l_830/edit#gid=124418455)*

Global parameters include:
   * Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
   * Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
   * bucket_name = 'bioscape.gra'
   * s3 = boto3.client('s3')
   * window_sizes = [60, 120, 240, 480, 960, 1200, 1500, 2000, 2200] (List of window sizes for computation of functional diversity metrics)
   * comps = 3 (Number of components to retain from principal component analysis)

The script includes:
*   Find all mosaics associated with a site (located in S3 bucket "SITENAME_flightlines/")
*   For each mosaic:
   * Load mosaic.
   * Prepare mosaic for dimensionality reduction:
        * Convert to numpy array.
        * Reshape and rescale data.
        * Impute values for NAs using "mean".
        * Scale and standarize according to mean and standard deviation.
        * Fit PCA on a subset of the data.
   * Transform with PCA.
   * Using a moving window, calculate functional richness and divergence for each window.
      * Functional richness is computed as the volume of the minimum convex hull that encompasses all PCA pixels in the window. (REF)
      * Functional divergence is computed as the distance of all PCA pixels in the window from the center of gravity. (REF)
   * Export FRic and FDiv for each window as a .csv.
   * Upload to S3 and delete local files.  

For a detailed walkthrough, see the notebook: [NA]

## Environmental Covariates
The fifth step of the workflow is to process the environmental covariates for each site to get summary characteristics for each plot (within a NEON site). For our implementation, there is one input:
1. The plot shapefiles (located in S3 bucket "Site_boundaries/")

**Objective:** Extract summaries of elevation, slope and canopy height from NEON data at the plot-level using [02_scripts/S06_Process_Covariates.py](https://github.com/mthayden4726/BioSCape_across_scales/blob/5411a8ab1337e47b79a0dbe6e3bc7141dde12c09/02_scripts/S06_Process_Covariates.py)

### Implementation:

This script requires the following input from users:
   1. Name of the NEON site in all caps (e.g., BART)
   2. Domain of the NEON site (e.g., D01)
   4. ID of NEON site (e.g., 4)
   5. Date of desired flights as YYYY-MM (e.g., 201908)
   6. Environmental raster of interest (DTM or CHM)

Example command: ```python 02_scripts/S06_Process_Covariates.py --SITECODE BART --DOMAIN D01 --ID_NO 5 --YEAR 2019-08 --ENV CHM```

*For a list of parameters associated with sites included in analysis, see the [2019 NEON Site List](https://docs.google.com/spreadsheets/d/17DJtV1BKtq0uLfcYCtM2kp7JjjjaWpuEV_jt95l_830/edit#gid=124418455)*

Global parameters include:
   * Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
   * Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
   * bucket_name = 'bioscape.gra'
   * s3 = boto3.client('s3')

The script includes:
*   Find all flightlines associated with the environmental covariate of interest and the area of interest (shapefile).
*   For each flightline, clip to the area of interest.
*   Mosaic together all clipped flightlines that align with the same area of interest.
*   Extract summary metrics for each mosaic, including mean, median, min, max, std, and var.
*   Export as .csv and upload to S3. 

For a detailed walkthrough, see the notebook: [NA]

## Scaling Analysis 
The final step of the workflow is to fit scaling relationships to the functional richness and divergence outputs and extract parameters (exponent and coefficient) from the model fits. For our implementation, there is one input:
1. The FRic and FDiv results files (located in S3 bucket "/")

**Objective:** Fit power law functions to the functional diversity metrics and extract parameters of interest using [02_scripts/S07_Scaling_Functions.py](NA)

### Implementation:
Example command: ``` python 02_scripts/S07_Scaling_Functions.py ```

SCRIPT IN PROGRESS - TRANSLATING FROM R AND FOR USE WITH S3. 
