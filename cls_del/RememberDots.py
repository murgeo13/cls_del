# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 22:44:40 2023

@author: tzar
"""
import numpy as np

def RememberDots(CHRCOV,
                 CHRNAME):
    """
    Preprocessing of cover file

    Parameters
    ----------
    CHRCOV : str
        path to cover created with
        "samtools mpileup ${bam} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${CHRCOV}"
    CHRNAME : str
        Name of chromosome for deletion search

    Returns
    -------
    cover : np.array
        cover for each position
    dots : np.array
        row for each read is 
        read's coordinate/paired read's coordinate/bam-flag

    """

    with open(CHRCOV) as inp:
        dots = set()
        poses = []
        covers = []
        for line in inp:
                line = line.split("\t")
                poses.append(int(line[1]))
                covers.append(int(line[3]))
                flags = line[6].strip().split(",")
                might_xs = line[7].strip().split(",")
                chrs = line[8].strip().split(",")
                might_ys = line[9].strip().split(",")
                for flag, x, chr_, y in zip(flags, might_xs, chrs, might_ys):
                    if chr_ == CHRNAME and y != "*":
                        dots.add(tuple([int(x), int(y), int(flag)]))
                        
    cover = np.array([poses, covers])                       
    dots = np.array(list(dots))
    
    return cover, dots