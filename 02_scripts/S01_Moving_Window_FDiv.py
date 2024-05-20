import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
import csv
from tqdm import tqdm

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

# Calculate functional divergence from PCA
def calculate_FDiv(pca):
    """Function to calculate functional divergence from a PCA.
    
    Parameters:
    -----------
    pca: PCA object.
        
    Returns:
    -----------
    FDiv: value of FDiv for the given PCA chunk.
    
    """
    # Calculate mean of PCA data along the first axis
    mean_pca = np.mean(pca, axis=0)
    # Euclidean distances to mean along the first axis
    dist_to_mean = np.linalg.norm(pca - mean_pca, axis=1)
    # Mean of distances
    meandB = np.mean(dist_to_mean)
    # Deviations to mean
    devdB = dist_to_mean - meandB
    # Computation of FDiv
    FDiv = (np.sum(devdB) / len(dist_to_mean) + meandB) / (np.sum(np.abs(devdB)) / len(dist_to_mean) + meandB)
    return FDiv


def window_calcs_fdiv(args):
    
    """ Calculate functional divergence using a PCA chunk and window size.
    FOR USE IN PARALLEL PROCESSING OF FUNCTIONAL EVENNESS.
    
    Parameters:
    -----------
    pca_chunk: PCA chunk from NEON image.
    window_sizes: list/array of integers
    
    Returns:
    -----------
    FDiv: functional divergence for given window size and image.
    
    """
    windows, pca_x, results_FD, local_file_path  = args
    for window in tqdm(windows, desc='Processing window for batch'):
        print(window)
        half_window = window // 2
        FDiv_values = []
        for i in tqdm(range(half_window, pca_x.shape[0] - half_window, 15), desc='Processing window index'):
            for j in range(half_window, pca_x.shape[1] - half_window, 15):
                sub_arr = pca_x[i - half_window:i + half_window + 1, j - half_window:j + half_window + 1, :]
                if sub_arr is None:
                    print(f"sub_arr is None at position ({i}, {j})")
                    continue  # Skip this iteration if sub_arr is None
                #print(f"sub_arr shape: {sub_arr.shape}")
                FDiv = calculate_FDiv(sub_arr)
                print(FDiv)
                if isinstance(FDiv, (int, float)):  # Check if FDiv is a single value
                    FDiv_values.append(FDiv)  # Append the single value to FDiv_values
                else:
                    FDiv_values.extend(FDiv)  # Append multiple values if FDiv is iterable
                
        with open(local_file_path, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['Window_Size', 'FDiv'])  # Write header
            for FDiv_value in FDiv_values:
                csvwriter.writerow([window, FDiv_value])
