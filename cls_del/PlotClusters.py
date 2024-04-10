# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 14:08:18 2024

@author: tzar
"""
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

def PlotCluster(d_clusters,
                OUTDIR = "./cls_plots"):
    """
    Plot each cluster

    Parameters
    ----------
    d_clusters : dict
        Dictionary of dots in cluster
    OUTDIR : str, optional
        Path to output directory for plots. The default is "./cls_plots".

    Returns
    -------
    None.

    """
          
    for i in d_clusters.keys():
        if i[0] == -1:
            continue
        cl = d_clusters[i]
        for dot in cl:
            x, y, flag = dot
            if flag in [65, 113, 129, 177]:
                purple = plt.scatter(x, y, s=5, c="purple", label='same direction') 
            elif flag in [83, 99, 147, 163]:
                green = plt.scatter(x, y, s=5, c="green", label='properly mapped') 
            elif flag in [97, 161]:
                red = plt.scatter(x, y, s=5, c="red", label='large TLEN -->') 
            elif flag in [81, 145]:
                blue = plt.scatter(x, y, s=5, c="blue", label='large TLEN <--')
            else:
                black = plt.scatter(x, y, s=5, c="black", label='other')
        plt.plot([np.min(cl[:,:2])-10, np.max(cl[:,:2])+10],
                 [np.min(cl[:,:2])-10, np.max(cl[:,:2])+10], linestyle='dashed')
        all_colours = []
        try:
            all_colours.append(purple)
            del purple
        except:
            pass
        try:
            all_colours.append(green)
            del green
        except:
            pass
        try:
            all_colours.append(blue)
            del blue
        except:
            pass
        try:
            all_colours.append(red)
            del red
        except:
            pass
        try:
            all_colours.append(black)
            del black
        except:
            pass
        plt.legend(handles=all_colours, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(f"{OUTDIR}/Cluster_{i}.png", bbox_inches="tight")
        plt.close()
        
    
    
