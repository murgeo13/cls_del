# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:52:58 2024

@author: tzar
"""
import numpy as np
#from double_clusterisation import double_clusterisation
from modifyed_DBSCAN  import modifyed_DBSCAN 


def moving_average(X, Y, window, circle=False):
    """
    Calculate average in moving window

    Parameters
    ----------
    X : np.array
        array of x_s
    Y : np.array
        array of y_s
    window : int
        size of window
    circle : bool, optional
        Shoud be date interpreated as circle. The default is False.

    Returns
    -------
    np.array
        (X, averge_Y).

    """
    X = np.array(X)
    Y = np.array(Y)
    new_X = np.convolve(X, np.ones(window), "valid") // window
    if circle:
        new_X = np.concatenate((X[window//2:], X[:window//2]))
        Y = np.concatenate((Y, Y[:window-1]))
    
    new_Y = np.convolve(Y, np.ones(window), "valid") / window  
    
    average = np.array([new_X, new_Y])
    average = average[:, average[0, :].argsort()]
    
    return average

def DoubleClusterDots(cover, dots, all_cover,
                      MEDIAN, 
                      SD, 
                      METRIC = "manhattan", 
                      MIN_SIZE = "auto", 
                      SELECTION = ["all"],
                      OUTTSV = "clusters.tsv",
                      max_cluster_size = "auto"):
    """
    

    Parameters
    ----------
    cover : np.array
        cover for each position for reads for clusterisation
    dots : np.array
        row for each read is 
        read's coordinate/paired read's coordinate/bam-flag
    all_cover : np.array
        cover for each position
    MEDIAN : float
        Median of size of insertion between reads in pair
    SD : float
        Standard deviation of size of insertion between reads in pair
    METRIC : str, optional
        Metric from sklearn.metrics.pairwise.pairwise_distances. The default is "manhattan".
    MIN_SIZE : int or function(x, y)
        Minimum amout of reads in cluster. For each read elected for clusterisation it might be
        a function of two coordinates (read and its pair). The defaut is local average cover.
    SELECTION : list, optional
        ["select", list_of_flags] or ["exclude", list_of_flags]. The default is ["all"].
    OUTTSV : str, optional
        Path to otput file. The default is "clusters.tsv".

    Returns
    -------
    average_covers : np.array
        Computed in window 2*3*SD
    average_all_covers : np.array
        Computed in window 2*3*SD
    d_clusters : dict
        Dictionary of dots in cluster

    """
    
    window = 2*3*SD
    average_covers = moving_average(cover[0], cover[1], window, circle=True)  
    average_all_covers = moving_average(all_cover[0], all_cover[1], window, circle=True) 
    
    X = dots[:,:2]
    flags = dots[:,2]
    
    if SELECTION[0] == "all":
        mask = np.ones_like(flags, dtype=bool) #select all
        
    if SELECTION[0] == "select":
        #select only [97, 161, 81, 145] flags
        mask = np.zeros_like(flags, dtype=bool)
        for i in SELECTION[1]:
            mask = mask | (flags == i)
    
    if SELECTION[0] == "exclude":
        #select all except properly mapped  
        mask = np.ones_like(flags, dtype=bool)
        for i in SELECTION[1]:
            mask = mask & (flags != i)   
    
    selected_X = X[mask]
    selected_flags = flags[mask]
    
    # make nupay array of clusters' and subclusters labels
    def func(x, y):
        return average_covers[1, x-1]
    if MIN_SIZE == "auto":
        MIN_SIZE = func
    if max_cluster_size == "auto":
        max_cluster_size=MEDIAN+2*SD
    
    db = modifyed_DBSCAN(eps=MEDIAN, min_samples=MIN_SIZE, metric=METRIC).fit(selected_X).split(max_cluster_size=max_cluster_size)
    d_l = db.labels_
    d_subl = db.sublabels_
    
    #delete clusters without discordant reads
    correct_flags = np.array([99, 147, 83, 163])
    for c in np.unique(d_l):
        sc = 0
        while True:
            mask = d_l == c
            submask = mask & (d_subl == sc)
            cls_flags = np.unique(selected_flags[submask])
            size = len(np.unique(np.concatenate((cls_flags, correct_flags))))
            if len(cls_flags) == 0:
                break
            if size == len(correct_flags):
                d_l[submask] = -1
                d_subl[submask] = 0
                d_subl -= mask & (d_subl > sc)
                sc -= 1
            sc += 1
    
    # make tsv with result of clusterisation
    ans = np.column_stack((selected_X, selected_flags, d_l, d_subl))
    with open(OUTTSV, "w") as out:
        print("First\tSecond\tFlag\tCluster\tSubcluster", file = out)
        for row in ans:
            print(*row, sep="\t", file=out) 
    
    # make dictionary of all subclusters including noise-cluster (-1, 0)
    d_clusters = dict()
    for dot, flag, cl, scl in zip(selected_X, selected_flags, d_l, d_subl):
        key = tuple([cl, scl])
        if key in d_clusters.keys():
            d_clusters[key] = np.row_stack((d_clusters[key], np.insert(dot, 2, flag)))
        else:
            d_clusters[key] = np.insert(dot, 2, flag)
    
       
    
    return average_covers, average_all_covers, d_clusters