# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 23:39:37 2024

@author: tzar
"""
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.utils._param_validation import InvalidParameterError
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import copy

def _calculate_cluster_size(X, metric):
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
    distances = pairwise_distances(X, X, metric=metric)
    return np.max(distances)

def _constant(C):
    def func(x, y):
        return C
    return func


class modifyed_DBSCAN(DBSCAN, KMeans):
    def __init__(self, **kwargs):            
        super().__init__(**kwargs)
        #self.set_params(**DBSCAN().get_params())
        self._parameter_constraints = DBSCAN._parameter_constraints.copy()
        del self._parameter_constraints["min_samples"]
        
    def get_params(self, deep=True):
        params = super().get_params(deep)
        # Hack to make get_params return base class params...
        cp = copy.copy(self)
        cp.__class__ = DBSCAN
        params.update(DBSCAN.get_params(cp, deep))
        return params
    
    def fit(self, X, y=None, sample_weight=None):
        self.X = X 
        if type(self.min_samples) is int:
            self.min_samples_func = _constant(self.min_samples)
        elif type(self.min_samples) is type(_calculate_cluster_size):
            self.min_samples_func = self.min_samples
            self.min_samples = np.array([self.min_samples_func(x, y) for x, y in self.X])
        elif type(self.min_samples) is np.ndarray:
            if self.X.shape[0] == self.min_samples.shape[0] == self.min_samples.size:
              self.min_samples_func = None  
            else:
               raise ValueError("min_samples is not correct") 
        else:
            raise ValueError("min_samples is not correct")
        super().fit(X, y, sample_weight)
        return self
    
    def split(self, max_cluster_size):
        labels = self.labels_
        unique_labels = set(labels)
        sublabels = np.zeros_like(labels)
        for c in unique_labels:
            if c == -1:
                continue
            mask = labels == c
            cluster = self.X[mask]
            distance = _calculate_cluster_size(cluster, self.metric)
            if distance > max_cluster_size:
                max_distance = distance
                i = 2
                while max_distance > max_cluster_size:
                    kmeans = KMeans(n_clusters=i, random_state=0).fit(cluster)
                    sublabels[mask] = kmeans.labels_
                    unique_sublabels = set(sublabels[mask])
                    max_temp_cl_size = 0
                    for temp_c in unique_sublabels:
                        submask = mask & (sublabels == temp_c)
                        subcluster = self.X[submask]
                        temp_cl_size = _calculate_cluster_size(subcluster, self.metric)
                        if temp_cl_size > max_temp_cl_size:
                            max_temp_cl_size = temp_cl_size 
                        
                    i += 1
                    max_distance = max_temp_cl_size
            sc = 0
            while True:
                mask = labels == c
                submask = mask & (sublabels == sc)
                size = sum(submask)
                try:
                    MIN_SIZE = min(self.min_samples[submask])
                except:
                    MIN_SIZE = self.min_samples
                if size == 0:
                    break
                if size < MIN_SIZE:
                    labels[submask] = -1
                    sublabels[submask] = 0
                    sublabels -= mask & (sublabels > sc)
                    sc -= 1
                sc += 1
            self.labels_ = labels
            self.sublabels_ = sublabels
        return self
        