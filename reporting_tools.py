#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)  # developmental hack, remove later!
import introSpect

__name__, __package__, invoked_directly = introSpect.cmdSupport(
    __name__, __package__, __file__
)
hint = introSpect.hint

import numpy as np
import scipy as sci
import pandas as pd
import seaborn as sns
from matplotlib import colors as mplCols
from typing import Union
from . import par_examples, plotting_tools

def recipe(*, table: Union[None, str] = None, verbose: bool = True,) -> dict:

    """
    Shortcuts to generate summaries and reports

    Parameters
    ----------
    table
        Path to a heatmap table.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    The heatmap image
    """

    ### Read the heatmap of p-values
    if table is None:
        ptable = pd.DataFrame.from_dict(par_examples.ptable)
    else:
        ptable = get_ptable(table)
    ptable = ptable.set_index('cohort')
    hint(verbose, "The -log10 p-values:\n", ptable.head())

    ### Create linkage matrix and reorder the table based on hierarchy
    row_linkage, col_linkage = tree_linkages(ptable)
    ptable, left_tree, bottom_tree = linkorder_table(ptable, row_linkage, col_linkage)
    hint(verbose, "Reordered table:\n", ptable.head())

    ### Create a latex table from the heatmap
    latex_hm = latex_heatmap(ptable)
    hint(verbose, "The LaTeX version:\n", latex_hm)

    ### Create an image of the heatmap
    image_hm = heatmap_image(ptable)

    return image_hm


def get_ptable(
    fn: str,
    *,
    idx: str = 'cohort',
    ) -> pd.DataFrame:

    """
    Shortcut to read the table of p-values into Pandas.

    Parameters
    ----------
    fn
        Path to the table.
    idx
        Column name of the row index.

    Returns
    -------
    Data frame with p-values.
    """
    
    ptable = pd.read_csv(fn, sep='\t')
    ptable = ptable.set_index(idx)
    return ptable

def heatmap_image(
    ptable: pd.DataFrame,
    *,
    row_linkage: Union[None, str] = None,
    col_linkage: Union[None, str] = None,
    ) -> sns.matrix.ClusterGrid:

    """
    Shortcut to ctreate a clustered heatmap with Seaborn.

    Parameters
    ----------
    ptable
        Table of P-values.
    row_linkage
        Row linkage.
    col_linkage
        Column linkage.

    Returns
    -------
    Data frame with p-values.
    """
    
    g = sns.clustermap(ptable, row_linkage=row_linkage, col_linkage=col_linkage, method='single', figsize=(8, 5.8),  xticklabels=True,  yticklabels=True)
    return g

def tree_linkages(
    ptable: pd.DataFrame,
    *,
    method: str = 'single',
    ) -> sns.matrix.ClusterGrid:

    """
    Shortcut to ctreate a clustered heatmap with Seaborn.

    Parameters
    ----------
    ptable
        Table of P-values.
    method
        Linkage method.

    Returns
    -------
    Linkage matrices for rows and columns.
    """

    row_linkage = sci.cluster.hierarchy.linkage(sci.spatial.distance.pdist(ptable), method=method)
    col_linkage = sci.cluster.hierarchy.linkage(sci.spatial.distance.pdist(ptable.T), method=method)
    return row_linkage, col_linkage

def linkorder_table(
    ptable: pd.DataFrame,
    row_linkage: Union[None, str] = None,
    col_linkage: Union[None, str] = None,
    ) -> sns.matrix.ClusterGrid:

    """
    Reorder the table based on hierarchical clustering.

    Parameters
    ----------
    ptable
        Table of P-values.
    row_linkage
        Row linkage.
    col_linkage
        Column linkage.

    Returns
    -------
    Reordered data frame.
    """
    
    ngenes = ptable.columns.values.tolist()
    tns = ptable.index.values.tolist()

    fig, ax1 = plotting_tools.plt.subplots(figsize=(0.3, len(tns)*0.263))
    row_order = sci.cluster.hierarchy.dendrogram(row_linkage, orientation='left', ax=ax1, link_color_func=lambda k: 'k', labels=tns)
    ax1.axis('off')

    fig, ax2 = plotting_tools.plt.subplots(figsize=(len(ngenes)*0.35, 0.3))
    col_order = sci.cluster.hierarchy.dendrogram(col_linkage, orientation='bottom', ax=ax2, link_color_func=lambda k: 'k', labels=ngenes)
    ax2.axis('off')

    ptable = ptable.iloc[row_order['leaves'],col_order['leaves']]

    return ptable, ax1, ax2

def latex_heatmap(
    ptable: pd.DataFrame,
    *,
    rCm: mplCols.LinearSegmentedColormap = plotting_tools.plt.get_cmap('Reds'),
    bCm: mplCols.LinearSegmentedColormap = plotting_tools.plt.get_cmap('Blues'),
    plot_extent: int = 1,
    ) -> str:

    """
    Reorder the table based on hierarchical clustering.

    Parameters
    ----------
    ptable
        Table of P-values.
    rCm
        Colormap for the background of the cell.
    bCm
        Colormap for the text in the cell.
    plot_extent
        How many pages a single group of figures span (how to increment link in heatmap).

    Returns
    -------
    LaTeX string for the heatmap table.
    """
    
    ngenes = ptable.columns.values.tolist()
    tns = ptable.index.values.tolist()
    t = '\\begin{tabular}{' + ' | '.join(['c'] * (1+len(ngenes))) + '}\n'
    t += ' & ' + ' & '.join(['\\rotatebox{90}{' + x +'}' for x in ngenes]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(plot_extent*i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    for i in range(len(tns)):
        row = ptable.iloc[i, :]
        t += tns[i] + ' & ' + ' & '.join(['\\color[rgb]{'+textPainter(x, bCm)+'} \\cellcolor[rgb]{'+cellPainter(x, rCm)+'} \\hyperlink{page.'+str(i+2)+'}{'+'{:0.2f}'.format(1*x) + '}' for x in row]) + '\\\\ \n'
    t += '\\end{tabular}'
    return t

def cellPainter(
    x: float,
    cm: mplCols.LinearSegmentedColormap,
    ) -> str:

    """
    Helper function producing a scaled background color.

    Parameters
    ----------
    x
        Log10 p-value.
    cm
        Colormap.

    Returns
    -------
    RGB color string.

    """
    
    c = cm(-0.2*(x+10**0.05))
    return ','.join([str(r) for r in c[:3]])

def textPainter(
    x: float,
    cm: mplCols.LinearSegmentedColormap,
    ) -> str:

    """
    Helper function producing a scaled text color.

    Parameters
    ----------
    x
        Log10 p-value.
    cm
        Colormap.

    Returns
    -------
    RGB color string.

    """
    
    c = cm(1-(-0.2*(x+10**0.05)))
    return ','.join([str(r) for r in c[:3]])

def main():

    """
    Define what should be executed when invoked from command line.
    """
    modified_kws = {
        "verbose": (
            0,
            "-v",
            "--verbose",
            {"dest": "verbose", "default": False, "action": "store_true"},
        ),
        "outFile": (
            1,
            "-o",
            "--outFile",
            {
                "dest": "outFile",
                "help": "Location where results should be saved. If not specified, STDOUT will be used.",
            },
        ),
    }
    mainFunction = introSpect.commandLines.cmdConnect(recipe, modified_kws)
    mainFunction.eval()
    mainFunction.save()
    return


__doc__ = recipe.__doc__
if invoked_directly:
    main()
