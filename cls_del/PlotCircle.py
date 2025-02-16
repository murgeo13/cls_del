# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:28:57 2024

@author: tzar
"""
import numpy as np
from pycirclize import Circos
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import scipy.stats as ss

def PlotCircle(average_covers, average_all_covers, labels, finds,
               CHRNAME,
               OUTDIR = "./cls_plots",
               plot_name ="Circle",
               chr_colour = "red",
               chr_name_colour = "white",
               cover_colour = "blue",
               all_cover_colour = "green",
               func=None,
               legend_title="Averege cover",
               **funckwargs):
    """
    Plot chromosome map with deletions coloured by liklihood.

    Parameters
    ----------
    average_covers : np.array
        Cover by reads used in clusterisation
    average_all_covers : np.array
        Cover by all reads
    finds : np.array
        row for each cluster includes
        First_mean First_std Second_mean Second_std Supporting_reads Total_reads Expected_count
    CHRNAME : str
        Name of chromosome.
    OUTDIR : str, optional
        Path to output directory for plots. The default is "./cls_plots".
    plot_name : str, optional
        Name for file with plot. The default is "Circle".
    chr_colour : str, optional
        DESCRIPTION. The default is "red".
    chr_name_colour : str, optional
        DESCRIPTION. The default is "white".
    cover_colour : str, optional
        DESCRIPTION. The default is "blue".
    all_cover_colour : str, optional
        DESCRIPTION. The default is "green".
    percentile : int, optional
        Percenteile to cut extremly high and low likliehood in colormap. The default is 10.
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap for gradient by finds' liklihoods. The default is "RdPu".
    func : function(labels, finds, circos)
        Function generating colours for each find and adding special legend to plot.
        Parameters
        ----------
        labels : np.array
            name of each cluster
        finds : list or other iterable
            row for each cluster includes
            First_mean First_std Second_mean Second_std Supporting_reads Total_reads Expected_count
        circos : pycirclize.circos.Circos
            plot of chromosome

        Returns
        -------
        colour_map : np.array
            colour for each find
        line_handles : list of Line2D
            being used for legend creation

    Returns
    -------
    None.

    """
    def default_func(labels, finds, circos, percentile=10, cmap="RdPu"):
        """
        Function generating colours for each find and adding special legend to plot.

        Parameters
        ----------
        labels : np.array
            name of each cluster
        finds : list or other iterable
            row for each cluster includes
            First_mean First_std Second_mean Second_std Supporting_reads Total_reads Expected_count
        circos : pycirclize.circos.Circos
            plot of chromosome

        Returns
        -------
        colour_map : np.array
            colour for each find
        line_handles : list of Line2D
            being used for legend creation

        """
        
        Ls = ss.binom.pmf(finds[:, 4], finds[:, 6], finds[:, 4]/finds[:, 6])
        Ls[np.isnan(Ls)] = 1
        logLs = np.log(Ls)
        
        sm = ScalarMappable(norm=Normalize(vmin=np.percentile(logLs, percentile),
                                           vmax=np.percentile(logLs, 100 - percentile    )),
                        cmap=cmap)
        colour_map = sm.to_rgba(logLs)
        circos.colorbar(vmin=np.percentile(logLs, percentile),
                        vmax=np.percentile(logLs, 100 - percentile),
                        cmap=cmap,
                        colorbar_kws=dict(label="log-liklihood"))
        line_handles = [
            Line2D([], [], linewidth=1, color=all_cover_colour, label="with all mapped reads"),
            Line2D([], [], linewidth=1, color=cover_colour, label="with reads used in clusterisation"),
        ]
        return colour_map, line_handles
    
    if func is None:
        func = default_func
    
    LEN = average_covers.shape[1]
      
    sectors = {f"{CHRNAME}" : LEN}
    circos = Circos(sectors, space=0)
    for sector in circos.sectors:
        track = sector.add_track((60, 70))
        track.axis(fc=chr_colour)
        track.text(sector.name, color=chr_name_colour, size=12)
        track.xticks_by_interval(10000)
    
    starts = finds[:,0]
    ends = finds[:,2]
    colours, line_handles = func(labels, finds, circos, **funckwargs)
    
    for t1, t2, colour in zip(starts, ends, colours):
        circos.link_line((f"{CHRNAME}", t1), (f"{CHRNAME}", t2),
                         lw=1, color=colour)
    track2 = circos.sectors[0].add_track((80, 100), r_pad_ratio=0.1)
    track2.axis()
    track2.line(average_all_covers[0,:], average_all_covers[1,:], lw=1, color=all_cover_colour)
    track2.line(average_covers[0,:], average_covers[1,:], lw=1, color=cover_colour)
    fig = circos.plotfig()    
    circos.ax.legend(
        handles=line_handles,
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        title=legend_title,
        handlelength=2,
        frameon=True)
    fig.savefig(f"{OUTDIR}/{plot_name}.svg", bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/{plot_name}.png", bbox_inches="tight")