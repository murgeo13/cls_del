# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 13:18:46 2024

@author: tzar
"""
import numpy as np
from disjoint_set import DisjointSet

def MakeFeatures(average_covers, average_all_covers, d_clusters, 
                 min_deletion_size = 0,
                 FEATURETSV = "features.tsv"):
    """
    Make a table with found deletions

    Parameters
    ----------
    average_covers : np.array
        Cover by reads used in clusterisation
    average_all_covers : np.array
        Cover by all reads
    d_clusters : dict
        Dictionary of dots in cluster
    min_deletion_size : int
        Minimal size of deletion for features filtration. By default there is no filtration.
    FEATURETSV : str, optional
        Path to output file with features. The default is "features.tsv".

    Returns
    -------
    finds : np.array
        row for each cluster includes
        First_mean First_std Second_mean Second_std Supporting_reads Total_reads Expected_count

    """
    
    ans = None
    for key, val in d_clusters.items():
        d_l, d_subl = key
        if d_l == -1:
            continue
        n_dim = len(val.shape)
        if n_dim == 1:
            val = np.array([val])
        elif n_dim == 2:
            val = val
        else:
            val = []
            print(f"Cluster {key} is wierd")  
        cls_X = val[:,:2]
        cls_flags = val[:,2]
        total = len(cls_flags)
        left = sum(np.isin(cls_flags, np.array([97, 161])))
        right = sum(np.isin(cls_flags, np.array([81, 145])))
        same = sum(np.isin(cls_flags, np.array([65, 113, 129, 177])))
        if left > 0:
            first_mean = np.mean(cls_X[:,0], dtype="int32")
            first_std = np.std(cls_X[:,0], dtype="int32")
            second_mean = np.mean(cls_X[:,1], dtype="int32")
            second_std = np.std(cls_X[:,1], dtype="int32")
            expected = average_covers[1, average_covers[0,:] == first_mean][0]
            row = np.array((d_l, d_subl, first_mean, first_std, second_mean, second_std, left, total, expected))
            if ans is None:
                ans = row
            else:
                ans = np.row_stack((ans, row))
        if right > 0:
            first_mean = np.mean(cls_X[:,1], dtype="int32")
            first_std = np.std(cls_X[:,1], dtype="int32")
            second_mean = np.mean(cls_X[:,0], dtype="int32")
            second_std = np.std(cls_X[:,0], dtype="int32")
            expected = average_covers[1, average_covers[0,:] == first_mean][0]
            row = np.array((d_l, d_subl, first_mean, first_std, second_mean, second_std, right, total, expected))
            if ans is None:
                ans = row
            else:
                ans = np.row_stack((ans, row))        
                
    ds = DisjointSet()
    for i in range(ans.shape[0]):
        ds.find(i)
        for j in range(ans.shape[0]):
                if ans[i,2]-ans[i,3] < ans[j,2] < ans[i,2]+ans[i,3] and ans[i,4]-ans[i,5] < ans[j,4] < ans[i,4]+ans[i,5]:
                    ds.union(i,j)
    
    finds = None   
    list_labels = []             
    for i in ds.itersets():
        group = ans[list(i)]
        labels = list(map(tuple, group[:,:2].astype(int)))
        first_mean = np.mean(group[:,2], dtype="int32")
        first_std = np.linalg.norm(group[:,3]).astype(int)
        second_mean = np.mean(group[:,4], dtype="int32")
        second_std = np.linalg.norm(group[:,5]).astype(int)
        reads = np.sum(group[:,6], dtype="int32")
        total = np.sum(group[:,7], dtype="int32")
        expected = np.sum(group[:,8], dtype="int32")
        row = np.array((first_mean, first_std, second_mean, second_std, reads, total, expected))
        deletion_size = abs(second_mean - first_mean)
        if deletion_size > min_deletion_size:
            list_labels.append(labels)
            if finds is None:
                finds = row
            else:
                finds = np.row_stack((finds, row))
              
    with open(FEATURETSV, "w") as out:
        print("Labels\tFirst_mean\tFirst_std\tSecond_mean\tSecond_std\tSupporting_reads\tTotal_reads\tExpected_count", file = out)
        for label, row in zip(list_labels, finds):
            print(label, *row, sep="\t", file=out) 
    
    return list_labels, finds