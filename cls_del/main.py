# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 23:15:55 2024

@author: tzar
"""
from RememberDots import RememberDots
from DoubleClusterDots import DoubleClusterDots
from PlotClusters import PlotCluster
from MakeFeatures import MakeFeatures
from PlotCircle import PlotCircle
import pickle
import argparse
import os


parser = argparse.ArgumentParser(
    description="If you used samtools mpileup --output-extra FLAG,POS,RNEXT,PNEXT and picard CollectInsertSizeMetrics, you can find deletions")
parser.add_argument("--chr", type=str, required=True,
                    help="Name of chromosome")
parser.add_argument("--chrcov", type=str, required=True,
                    help="Path to chromosom cover")
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
parser.add_argument("--out", type=str, required=False, default="./",
                    help='Path to output default="./"')
group1 = parser.add_argument_group("Seletion of flags",
                                  "If both --select and --exlude are unused, all flags will selected").add_mutually_exclusive_group()
group1.add_argument("--select", type=int, required=False, default=[], nargs="+",
                   help="Select ONLY flags in LIST")
group1.add_argument("--exclude", type=int, required=False, default=[], nargs="+",
                   help="Select all exept flags in LIST")
parser.add_argument("--save-temp-files", action="store_true", required=False, default=False,
                    help="Saving of temporary .pikle files")
parser.add_argument("--debug", action="store_true", required=False, default=False,
                    help="In debugging mode existing temporary .pickle files will be used")

args = parser.parse_args()
os.makedirs(args.out, exist_ok=True)

CHRCOV = args.chrcov
CHRNAME = args.chr
CHRCOV_TOTAL = args.chrcov_total
if CHRCOV_TOTAL is None:
    CHRCOV_TOTAL = CHRCOV
DOTSPICKLE = f"{args.out}/dots.pickle"
DOTSPICKLE_TOTAL = f"{args.out}/total_dots.pickle"
SELECTION = ["all"]
if args.select: SELECTION = ["select", args.select]
if args.exclude: SELECTION = ["exclude", args.exclude]
MEDIAN = args.median
SD = args.sd
METRIC = args.metric
MIN_SIZE = args.minsize
OUTTSV = f"{args.out}/clusters.tsv"
CLSPICKLE = f"{args.out}/clusters.pickle"
OUTDIR = f"{args.out}/cls_plots"
FEATURETSV = f"{args.out}/features.tsv"
FEATURESPICKLE = f"{args.out}/features.pikle"

os.makedirs(OUTDIR, exist_ok=True)

if args.debug:
    print("Loading variables from existing .pikle")
    try:
        with open(DOTSPICKLE, 'rb') as inp:
            cover = pickle.load(inp)
            dots = pickle.load(inp)  
        with open(DOTSPICKLE_TOTAL, 'rb') as inp:
            all_cover = pickle.load(inp)
            all_dots = pickle.load(inp)
        with open(CLSPICKLE, 'rb') as inp:
            average_covers = pickle.load(inp)
            average_all_covers = pickle.load(inp)
            d_clusters = pickle.load(inp)
        with open(FEATURESPICKLE, 'rb') as inp:
            labels = pickle.load(inp)
            finds = pickle.load(inp)
    except:
        pass


try:
    cover, dots
except:
    print("RememberDots(CHRCOV, CHRNAME) is started")            
    cover, dots = RememberDots(CHRCOV, CHRNAME)
    with open(DOTSPICKLE, 'wb') as out:
        pickle.dump(cover, out)
        pickle.dump(dots, out)
    print("RememberDots(CHRCOV, CHRNAME) is done")

try:
    all_cover, all_dots
except: 
    print("RememberDots(CHRCOV_TOTAL, CHRNAME) is started")
    all_cover, all_dots = RememberDots(CHRCOV_TOTAL, CHRNAME)
    with open(DOTSPICKLE_TOTAL, 'wb') as out:
        pickle.dump(all_cover, out)
        pickle.dump(all_dots, out)
    print("RememberDots(CHRCOV_TOTAL, CHRNAME) is done")
    
try:
    average_covers, average_all_covers, d_clusters
except:
    print("DoubleClusterDots(*args) is started")
    average_covers, average_all_covers, d_clusters = DoubleClusterDots(cover, dots, all_cover,
                                                                       MEDIAN, SD, METRIC, MIN_SIZE,
                                                                       SELECTION, OUTTSV)
    with open(CLSPICKLE, 'wb') as out:
        pickle.dump(average_covers, out)
        pickle.dump(average_all_covers, out)
        pickle.dump(d_clusters, out)
    print("DoubleClusterDots(*args) is done")

PlotCluster(d_clusters, OUTDIR)

try:
    labels, finds
except:
    print("MakeFeatures(*args) is started")
    labels, finds = MakeFeatures(average_covers, average_all_covers, d_clusters, FEATURETSV)
    with open(FEATURESPICKLE, 'wb') as out:
        pickle.dump(labels, out)
        pickle.dump(finds, out)
    print("MakeFeatures(*args) is done")

print("PlotCircle(*args) is started")    
PlotCircle(average_covers, average_all_covers, labels, finds, CHRNAME, OUTDIR)
print("PlotCircle(*args) is started")

if not args.save_temp_files:
    print("Removing all temportary files")
    os.remove(DOTSPICKLE)
    os.remove(CLSPICKLE)
    os.remove(FEATURESPICKLE)