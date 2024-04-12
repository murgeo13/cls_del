# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 23:56:28 2024

@author: tzar
"""

import numpy as np
import argparse
from pycirclize import Circos

parser = argparse.ArgumentParser(
    description="If you used delly call -o ${delly_out_bcf} -g ${ref} ${bam_on_chr} and bcftools view ${delly_out_bcf} > ${delly_out_vcf}")
parser.add_argument("--delly-vcf", type=str, required=True,
                    help="Path to .vcf")
parser.add_argument("--out", type=str, required=False, default="./",
                    help='Path to output default="./"')
parser.add_argument("--chrcov-total", type=str, required=False, default=None,
                    help="Path to total chromosom cover")
parser.add_argument("--median", type=int, required=True,
                    help="Median insert size (TLEN)")
parser.add_argument("--sd", type=int, required=True,
                    help="SD of insert size (TLEN)")
parser.add_argument("--metric", type=str, required=False, default="manhattan",
                    help="see valid values in sklearn.metrics.pairwise.pairwise_distances, default=manhattan")
parser.add_argument("--minsize", type=int, required=False, default=5,
                    help="Minimal size of cluster, default=5")

group1 = parser.add_argument_group("Seletion of flags",
                                  "If both --select and --exlude are unused, all flags will selected").add_mutually_exclusive_group()
group1.add_argument("--select", type=list, required=False, default=[],
                   help="Select ONLY flags in LIST")
group1.add_argument("--exclude", type=list, required=False, default=[],
                   help="Select all exept flags in LIST")
parser.add_argument("--save-temp-files", action="store_true", required=False, default=False,
                    help="Saving of temporary .pikle files")
parser.add_argument("--debug", action="store_true", required=False, default=False,
                    help="In debugging mode existing temporary .pickle files will be used")

args = parser.parse_args()

list_labels = list()
finds = None

with open(args.delly_vcf, "r") as inp:
    line = inp.readline()
    while line:
        if not line.startswith("#"):
            line = line.strip().split("\t")
            start = line[1]
            name = line[2]
            end = line[7]
            if name.startswith("BND"):
                end = None
            else:
                end = end.strip().split(";")
                end = end[4].split("=")[1]
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
    circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        title="Delly feature",
        handlelength=2,
    )
    return colour_map


PlotCircle(average_covers, average_all_covers, labels, finds,
               CHRNAME,
               OUTDIR = "./cls_plots",
               chr_colour = "red",
               chr_name_colour = "white",
               cover_colour = "blue",
               all_cover_colour = "green",
               percentile = 10,
               cmap = "RdPu",
               func=None)