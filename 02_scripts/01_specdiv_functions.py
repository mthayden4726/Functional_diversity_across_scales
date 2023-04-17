#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script containing functions for calculating spectral diversity metrics.
- Alpha diversity
- Beta diversity
- Functional richness

Additionally contains processing functions to:
    - plot an RGB map of multi-band raster.
    - 
    

Author: M. Hayden
Date: 4/12/2023
"""
def find_neon_files(SITECODE = site, PRODUCTCODE = product, 
                    YEAR = year):
    """Function to create API URLs for NEON images of interest.
    
    Parameters:
    -----------
    site: 
        String...
    product: 
        String...
    year:
        String...
        
    Returns:
    -----------
    file_paths: list of file paths for desired NEON images.
        URL for API where image is located.
    
    
    """
    # Server's URL will remain the same
    SERVER = 'http://data.neonscience.org/api/v0/'

    # Build so that we can loop through reading files in
    url = SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/'+YEAR
    
    # Request the url
    data_request = requests.get(url)

    # Convert the request to Python JSON object
    data_json = data_request.json()

    # Create list of file paths of interest for given site, product, year
    file_paths = []
    for file in data_json['data']['files'][:20]:
      if 'reflectance.h5' in file['name']:
            file_paths.insert(1, file['url'])
            print(file['url'])

def show_rgb(x,r=660,g=550,b=440, correct= []):
    """Display raster in RGB.
    
    Parameters:
    -----------
    x: hytools object
        Hytools object of image of interest.
    r: int
        Wavelength/band corresponding to red.
    g: int
        Wavelength/band corresponding to green.
    b: int
        Wavelength/band corresponding to blue.
    Returns:
    -----------
    rgb : raster
        RGB raster for display in plot viewer.
    
    Plots rgb raster.
    
    """
    gb=  np.stack([hy_obj.get_wave(r,corrections= correct),
                    hy_obj.get_wave(g,corrections= correct),
                    hy_obj.get_wave(b,corrections= correct)])
    rgb = np.moveaxis(rgb,0,-1).astype(float)
    rgb[rgb ==hy_obj.no_data] = np.nan

    bottom = np.nanpercentile(rgb,5,axis = (0,1))
    top = np.nanpercentile(rgb,95,axis = (0,1))
    rgb = np.clip(rgb,bottom,top)

    rgb = (rgb-np.nanmin(rgb,axis=(0,1)))/(np.nanmax(rgb,axis= (0,1))-np.nanmin(rgb,axis= (0,1)))

   # height = int(hy_obj.lines/hy_obj.columns)

    # fig  = plt.figure(figsize = (7,7) )
    plt.imshow(rgb)
    plt.show()
    plt.close()

def subsample():
def prog_bar():
    
def calc_shannon(arr):
    """ Calculate shannon index for an image.
    
    Parameters:
    -----------
    arr: array
    
    Returns:
    -----------
    value : shannon index for given window size and image.
    
    """
    p = np.bincount(arr.ravel()) / float(arr.size) # proportion of community made of species i 
    return -np.sum(p * np.log2(p), axis=0, where=(p > 0)) # sum of proportions times their natural logs

def calc_alpha(arr, windows):
    
    """ Calculate shannon index for different window sizes in an image.
    
    Parameters:
    -----------
    arr: array
    
    windows: list of window sizes (or single value)
    
    Returns:
    -----------
    results: dict of shannon values per window size
    
    """
    results = {} # storing results as a dict
    
    for window in windows:
        half_window = window // 2
        shannon = np.zeros(arr.shape)
        for i in range(half_window, arr.shape[0]-half_window):
            for j in range(half_window, arr.shape[1]-half_window):
                sub_arr = arr[i-half_window:i+half_window+1, j-half_window:j+half_window+1]
                shannon[i, j] = calc_shannon(sub_arr)
        results[window] = shannon
    return results

def calc_chv(arr):
    """ Calculate convex hull volume for an array.
    
    Parameters:
    -----------
    arr: 
        a numpy array of points, where each row is a pixel and each 
        column is a trait.
    Returns:
    -----------
    value : shannon index for given window size and image.
    
    """
    hull = ConvexHull(points = arr)
    volume = hull.volume
    return volume
    
def calc_fun_rich(neon, window_sizes):
    """ Calculate convex hull volume for an array at a variety of window sizes.
    
    Parameters:
    -----------
    neon: hytools image
    window_sizes: list/array of integers
    
    Returns:
    -----------
    volume_mean: functional richness for given window size and image.
    
    """
    volumes = {}
    results_FR = {}
    iterator = neon.iterate(by = 'chunk',chunk_size = (500,500))
    while not iterator.complete:
        chunk = iterator.read_next()
        X_chunk = chunk[:,:,~neon.bad_bands].astype(np.float32)
        X_chunk = X_chunk.reshape((X_chunk.shape[0]*X_chunk.shape[1],X_chunk.shape[2]))
        X_chunk -=x_mean
        X_chunk /=x_std
        X_chunk[np.isnan(X_chunk) | np.isinf(X_chunk)] = 0
        pca_chunk=  pca.transform(X_chunk)
        pca_chunk = pca_chunk.reshape((chunk.shape[0],chunk.shape[1],comps))
        pca_chunk[chunk[:,:,0] == neon.no_data] =0
        for window in window_sizes:
            half_window = window // 2
            fric = np.zeros(pca_chunk.shape)
            for i in range(half_window, pca_chunk.shape[0]-half_window):
                for j in range(half_window, pca_chunk.shape[1]-half_window):
                    sub_arr = pca_chunk[i-half_window:i+half_window+1, j-half_window:j+half_window+1, :]
            sub_arr = sub_arr.reshape((sub_arr.shape[0]*sub_arr.shape[1],comps))
            if np.nanmean(sub_arr) == 0.0:
                continue
            hull = ConvexHull(sub_arr)
            fric = hull.volume
            results_FR[iterator.current_line, window] = fric
        volumes[iterator.current_line] = np.array(list(results_FR.values())).mean()
    volume_mean = np.array(list(volumes.values())).mean()
    return volume_mean

    
#### In progress below this point ####



window = 15 # desired window size (# of pixels)
classes = classes # array of k-means cluster value assigned to each pixel
nclusters = nclusters # number of clusters used in k-means
def calc_alpha(classes,window,nclusters):
    shannon = np.zeros(classes.shape) # makes blank matrix the size of the array
    simpson = np.zeros(classes.shape) # makes blank matrix the size of the array
    windows_pixels = window**2 # stores shape of window

    lines = classes.shape[0] # stores number of lines in image
    columns = classes.shape[1] # stores number of columns in image

    half_window = int(window/2)

    for line in range(half_window,lines-half_window): # from a half window down into image to a half window from bottom of image
        for col in range(half_window,columns-half_window): # and from a half into into image to a half window from end of side of image
            nbhd = classes[line-half_window:line+half_window,col-half_window:col+half_window].flatten()
            cover = [0 for x in range(nclusters)]
            nbhd_size = 0

            for a in nbhd:
                if a !=nclusters:
                    cover[a]+=1
                    nbhd_size+=1

            if nbhd_size/windows_pixels < .75:
                continue

            shn,smp = 0,0

            for c in range(nclusters):
                if cover[c] !=0:
                    p = cover[c]/nbhd_size
                    shn += p * np.log(p)
                    smp += p**2
            shannon[line,col] = -shn
            simpson[line,col] = 1/smp
    return shannon,simpson


def calc_alpha(classes, windows, nclusters):
    shannon = np.zeros((len(windows), *classes.shape))
    simpson = np.zeros((len(windows), *classes.shape))
    windows_pixels = windows[:, np.newaxis, np.newaxis] ** 2

    lines = classes.shape[0]
    columns = classes.shape[1]

    half_windows = windows // 2

    for i, window in enumerate(windows):
        for line in range(half_windows[i], lines - half_windows[i]):
            for col in range(half_windows[i], columns - half_windows[i]):
                nbhd = classes[:, line - half_windows[i]:line + half_windows[i] + 1, col - half_windows[i]:col + half_windows[i] + 1].reshape((classes.shape[0], -1))
                cover = np.zeros((classes.shape[0], nclusters))
                nbhd_size = np.zeros(classes.shape[0])

                for a in range(classes.shape[0]):
                    mask = nbhd[a] != nclusters
                    cover[a][nbhd[a][mask]] = np.bincount(nbhd[a][mask], minlength=nclusters)
                    nbhd_size[a] = np.sum(mask)

                if np.any(nbhd_size / windows_pixels[i] < .75):
                    continue

                shn, smp = 0, 0

                for c in range(nclusters):
                    p = cover[:, c] / nbhd_size
                    mask = p > 0
                    shn += np.sum(p[mask] * np.log(p[mask]))
                    smp += np.sum(p[mask] ** 2)

                shannon[i][line][col] = -shn
                simpson[i][line][col] = 1 / smp

    return shannon, simpson

def calc_alpha(classes, windows, nclusters):
    shannon = np.zeros((len(windows), *classes.shape))
    simpson = np.zeros((len(windows), *classes.shape))
    windows_pixels = windows[:, np.newaxis, np.newaxis] ** 2

    lines = classes.shape[0]
    columns = classes.shape[1]

    half_windows = windows // 2

    # create an array of indices for each window size
    window_indices = np.indices((len(windows), lines - np.max(half_windows), columns - np.max(half_windows)))

    # offset the indices by the appropriate half window size
    window_indices[:, :] += half_windows[:, np.newaxis, np.newaxis]

    # create an array of neighborhood indices for each window size
    nbhd_indices = np.indices((len(windows), *windows_pixels.shape))

    # offset the neighborhood indices by the appropriate half window size
    nbhd_indices[1:] -= half_windows[:, np.newaxis, np.newaxis]

    # reshape the neighborhood indices to a 2D array
    nbhd_indices = nbhd_indices.reshape((len(windows), -1, nbhd_indices.shape[-1]))

    # create a mask for pixels with the right neighborhood size
    mask = (nbhd_indices.shape[-1] - np.sum(nbhd_indices == 0, axis=-1)) / windows_pixels >= 0.75

    # apply the mask to the neighborhood indices
    nbhd_indices = nbhd_indices[:, mask, :]

    # create a flat index for the neighborhood
    nbhd_flat_index = np.ravel_multi_index(nbhd_indices, classes.shape)

    # index the classes array with the neighborhood index
    nbhd = classes.ravel()[nbhd_flat_index]

    # reshape the neighborhood array to a 2D array
    nbhd = nbhd.reshape((len(windows), -1, nbhd.shape[-1]))

    # count the number of occurrences of each cluster in the neighborhood
    cover = np.zeros((len(windows), nbhd.shape[1], nclusters))
    nbhd_size = np.sum(nbhd != nclusters, axis=-1)
    mask = nbhd_size > 0
    for c in range(nclusters):
        cover[:, mask, c] = np.sum(nbhd == c, axis=-1)[:, mask] / nbhd_size[:, mask]

    # calculate Shannon entropy and Simpson diversity
    shn = -np.sum(cover * np.log(cover), axis=-1)
    smp = np.sum(cover ** 2, axis=-1)

    # fill in the output arrays with the calculated values
    shannon[:, half_windows[0]:lines-half_windows[0], half_windows[0]:columns-half_windows[0]][mask] = shn[mask]
    simpson[:, half_windows[0]:lines-half_windows[0], half_windows[0]:columns-half_windows[0]][mask] = 1 / smp[mask]

    return shannon, simpson


# CHATGPT solution for efficent loading of imagery and conversion to array
import requests
import numpy as np

# Make a GET request to the API to retrieve the image data
response = requests.get('https://example.com/hyperspectral_image')

# Convert the response content to a numpy array
image_data = np.frombuffer(response.content, dtype=np.uint16)

# Reshape the array to the correct dimensions (1389, 8691)
image_array = np.reshape(image_data, (1389, 8691))
