# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 18:42:38 2024

@author: tzar
"""
from sklearn.cluster import DBSCAN 
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np

def double_clusterisation(selected_X, MEDIAN, SD, METRIC, MIN_SIZE):
    """
    Searching of unproperly mapped reads' clusters, that might be associated with deletion

    Parameters
    ----------
    selected_X : np.array
        For each read selected for clusterisation its coordinate and its paired read coordinate
    MEDIAN : float
        Median of size of insertion between reads in pair
    SD : float
        Standard deviation of size of insertion between reads in pair
    METRIC : str
        Metric from sklearn.metrics.pairwise.pairwise_distances
    MIN_SIZE : int
        Minimum amout of reads in cluster

    Returns
    -------
    labels : np.array
        Cluster label for each read
    sublabels : np.array
        Subcluster label for each read
    
    Note that cluster divided to subclusters untill each subcluster size <= MEDIAN + 2*SD
    Noize is labeld as (-1, 0)

    """

    def calculate_cluster_size(X):
        """
        Calculate cluster's size using METRIC

        Parameters
        ----------
        X : np.array
            cluster

        Returns
        -------
        size : float
            Cluster's size

        """
        distances = pairwise_distances(X, X, metric=METRIC)
        return np.max(distances)
    
    db = DBSCAN(eps=MEDIAN, min_samples=MIN_SIZE, metric=METRIC).fit(selected_X)
    labels = db.labels_
    
    unique_labels = set(labels)
    
    sublabels = np.zeros_like(labels)
    
    for c in unique_labels:
        if c == -1:
            continue
        mask = labels == c
        cluster = selected_X[mask]
        distance = calculate_cluster_size(cluster)
        if distance > MEDIAN + 2*SD:
            max_distance = distance
            i = 2
            while max_distance > MEDIAN + 2*SD:
                kmeans = KMeans(n_clusters=i, random_state=0).fit(cluster)
                sublabels[mask] = kmeans.labels_
                unique_sublabels = set(sublabels[mask])
                max_temp_cl_size = 0
                for temp_c in unique_sublabels:
                    submask = mask & (sublabels == temp_c)
                    subcluster = selected_X[submask]
                    temp_cl_size = calculate_cluster_size(subcluster)
                    if temp_cl_size > max_temp_cl_size:
                        max_temp_cl_size = temp_cl_size 
                    
                i += 1
                max_distance = max_temp_cl_size
    
    for c in np.unique(labels):
        sc = 0
        while True:
            mask = labels == c
            submask = mask & (sublabels == sc)
            size = sum(submask)
            if size == 0:
                break
            if size < MIN_SIZE:
                labels[submask] = -1
                sublabels[submask] = 0
                sublabels -= mask & (sublabels > sc)
                sc -= 1
            sc += 1
            
        
    return labels, sublabels
    


