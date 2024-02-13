import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

# Calculate functional divergence from PCA
def calculate_FDiv_from_PCA(pca):
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

# Calculate functional evenness from PCA

# First create vertex information (assign each row in PCA to a vertex)
vertex_info = pca # pca must be an array
nb = pca.shape[0]

# Compute pairwise Euclidean distances between vertices
distances = pdist(vertex_info)

# Convert the condensed distance matrix to a square distance matrix
dist_matrix = squareform(distances)

# Compute minimum spanning tree
mst = minimum_spanning_tree(dist_matrix)
mstvect = mst.toarray()

def calculate_FEve(mstvect, dist_matrix, nb):
    # computation of EW for the (nb - 1) segments to link the nb points
    EW = np.zeros(nb - 1)
    flag = 0
    for m in range((nb - 1) * nb // 2):
        if mstvect[m] != 0:
            EW[flag] = dist_matrix[m]
            flag += 1
    return EW
