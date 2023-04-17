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


def show_rgb(hy_obj,r=660,g=550,b=440, correct= []):
    """Display raster in RGB.
    
    Parameters:
    -----------
    hy_obj: hytools object
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
    rgb=  np.stack([hy_obj.get_wave(r,corrections= correct),
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

    
# Trying the ChatGPT version
def shannon_index(arr):
    p = np.bincount(arr.ravel()) / float(arr.size)
    return -np.sum(p * np.log2(p), axis=0, where=(p > 0))

window_sizes = [25, 35, 50, 75, 100, 150]  # list of window sizes to test
results = {} # storing results as a dict

arr = classes
for window in window_sizes:
    half_window = window // 2
    shannon = np.zeros(arr.shape)

    for i in range(half_window, arr.shape[0]-half_window):
        for j in range(half_window, arr.shape[1]-half_window):
            sub_arr = arr[i-half_window:i+half_window+1, j-half_window:j+half_window+1]
            shannon[i, j] = shannon_index(sub_arr)

    results[f"window_{window}"] = shannon

print(results)

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
