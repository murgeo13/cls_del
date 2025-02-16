# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 23:56:28 2024

@author: tzar
"""

import numpy as np
import argparse
import pickle
from matplotlib.lines import Line2D
from RememberDots import RememberDots
from DoubleClusterDots import moving_average
from PlotCircle import PlotCircle
import os

parser = argparse.ArgumentParser(
    description="If you used delly call -o ${delly_out_bcf} -g ${ref} ${bam_on_chr} and bcftools view ${delly_out_bcf} > ${delly_out_vcf}")
parser.add_argument("--delly-vcf", type=str, required=True,
                    help="Path to .vcf")
parser.add_argument("--out", type=str, required=False, default="./",
                    help='Path to output default="./"')
parser.add_argument("--chrname", type=str, required=False, default="mt",
                    help='Path to output default="./"')
parser.add_argument("--chrcov-total", type=str, required=True,
                    help="Path to total chromosom cover")
parser.add_argument("--window", type=int, required=False, default=1000,
                    help="Size of window for moving average, "
                    "default=1000")

parser.add_argument("--save-temp-files", action="store_true", required=False, default=False,
                    help="Saving of temporary .pikle files")
parser.add_argument("--debug", action="store_true", required=False, default=False,
                    help="In debugging mode existing temporary .pickle files will be used")

args = parser.parse_args()

window = args.window
CHRNAME = args.chrname
CHRCOV_TOTAL = args.chrcov_total
CLSPICKLE = f"{args.out}/clusters.pickle"
DOTSPICKLE_TOTAL = f"{args.out}/total_dots.pickle"

if args.debug:
    print("Loading variables from existing .pikle")
    try:
        with open(DOTSPICKLE_TOTAL, 'rb') as inp:
            all_cover = pickle.load(inp)
        with open(CLSPICKLE, 'rb') as inp:
            average_covers = pickle.load(inp)
            average_all_covers = pickle.load(inp)
    except:
        pass

try:
    all_cover
except: 
    print("RememberDots(CHRCOV_TOTAL, CHRNAME) is started")
    all_cover, all_dots = RememberDots(CHRCOV_TOTAL, CHRNAME)
    with open(DOTSPICKLE_TOTAL, 'wb') as out:
        pickle.dump(all_cover, out)
        pickle.dump(all_dots, out)
    print("RememberDots(CHRCOV_TOTAL, CHRNAME) is done")
try:
    average_covers, average_all_covers
except:    
    average_all_covers = moving_average(all_cover[0], all_cover[1], window, circle=True)



list_labels = list()
finds = None

with open(args.delly_vcf, "r") as inp:
    line = inp.readline()
    while line:
        if not line.startswith("#"):
            line = line.strip().split("\t")
            start = int(line[1])
            name = line[2]
            end = line[7]
            if name.startswith("BND"):
                end = None
            else:
                end = end.strip().split(";")
                end = int(end[4].split("=")[1])
            list_labels.append(name)
            row = np.array([start, 0, end, 0])
            if finds is None:
                finds = row
            else:
                finds = np.row_stack((finds, row))
        line = inp.readline()

with open(f"{args.out}/delly_features.tsv", "a") as out:
    print("Labels\tFirst_mean\tFirst_std\tSecond_mean\tSecond_std\tSupporting_reads\tTotal_reads\tExpected_count", file = out)
    for label, row in zip(list_labels, finds):
        print(label, *row, sep="\t", file=out) 


def delly_func(labels, finds, circos):
    names = labels
    
    colour_map = ["green" if name.startswith("DUP") else
     "yellow" if name.startswith("INV") else
     "red" if name.startswith("DEL") else None
     for name in names]
    
    line_handles = [
        Line2D([], [], color="green", label="DUP"),
        Line2D([], [], color="yellow", label="INV"),
        Line2D([], [], color="red", label="DEL"),
    ]
    return colour_map, line_handles
    
PlotCircle(average_covers, average_all_covers, list_labels, finds,
               CHRNAME,
               OUTDIR = f"{args.out}",
               plot_name ="delly",
               chr_colour = "red",
               chr_name_colour = "white",
               cover_colour = "blue",
               all_cover_colour = "green",
               func=delly_func,
               legend_title="Delly feature")

if not args.save_temp_files:
    print("Removing all temportary files")
    os.remove(DOTSPICKLE_TOTAL)
    os.remove(CLSPICKLE)