#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)  # developmental hack, remove later!
import introSpect

request_for_special_plotting_backend = os.environ.get("MPLBACK", None)
__name__, __package__, invoked_directly = introSpect.cmdSupport(
    __name__, __package__, __file__
)
hint = introSpect.hint

import matplotlib

if request_for_special_plotting_backend is not None:
    matplotlib.use(request_for_special_plotting_backend)
from matplotlib import pyplot as plt
import seaborn as sns
from typing import Union, Tuple


def recipe(
    *,
    fontsizes: Union[None, dict] = None,
    style_kws: Union[None, dict] = None,
    figsize: Union[None, Tuple[float, float]] = None,
    remove_spines: bool = True,
    verbose: bool = True,
) -> dict:

    """
    Helper functions to facilitate plotting and styling of plots.

    Parameters
    ----------
    fontsizes
        Dictionary of font size for small, medium and large text.
    style_kws
        Settings that will be passed on to the Pyplot rc params.
    figsize
        Dimensions of the figure in inches. Does not override value defined in style_kws.
    remove_spines
        Right and top spines are removed by default.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    The parameters that were changed from Matplotlib defaults.
    
    """

    ### Set x tick color to red
    style_kws = set_figure_rc(style_kws={"xtick.color": "red"})
    hint(verbose, "Examples supplied with this package:\n", style_kws)

    ### Show style of plot
    plt.scatter([7, 3, 2, 5, 3], [2, 6, 8, 1, 3])
    plt.savefig("tmp.png")

    ### Make a plot with legend only
    legend_only()
    return style_kws


def legend_only(
    *,
    labels: Union[None, Tuple[str, str]] = ("low expression", "high expression"),
    ax: Union[None, plt.Axes] = None,
) -> plt.Axes:

    """
    Creates an empyt plot with the legend only.

    Parameters
    ----------
    labels
        Legend label for the low and the high expression.
    ax
        The matplotlib axis object for the plot.

    Returns
    -------
    The matplotlib axis object with the legend.
    """

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(range(-4, -1), range(-4, -1), label=labels[0])
    ax.plot(range(-4, -1), range(-4, -1), label=labels[1])
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.legend(title="", loc="upper left", frameon=False)
    return ax


def set_figure_rc(
    *,
    fontsizes: Union[None, dict] = None,
    style_kws: Union[None, dict] = None,
    figsize: Union[None, Tuple[float, float]] = None,
    remove_spines: bool = True,
) -> dict:

    """
    Sets rcParams to match a custom style.

    Parameters
    ----------
    fontsizes
        Dictionary of font size for small, medium and large text.
    style_kws
        Settings that will be passed on to the Pyplot rc params.
    figsize
        Dimensions of the figure in inches. Does not override value defined in style_kws.
    remove_spines
        Right and top spines are removed by default.
    isnotebook
        The pgf backend does not work well with Jupyter; check if called from notebook.
    
    Returns
    -------
    The parameters that were changed from Matplotlib defaults.
    """

    if fontsizes is None:
        fontsizes = dict()
    _fontsizes = fontsizes.copy()
    fontsizes = {"SMALL_SIZE": 6, "MEDIUM_SIZE": 7, "BIGGER_SIZE": 8}
    fontsizes.update(_fontsizes)

    default_style = {
        "legend.frameon": False,
        "savefig.transparent": True,
        "figure.subplot.hspace": 0.6,
        "figure.subplot.wspace": 0.3,
        "figure.max_open_warning": -1,
        "figure.dpi": 600,
        "axes.labelpad": -0.15 * fontsizes["SMALL_SIZE"],
        "figure.titlesize": fontsizes["BIGGER_SIZE"],
        "axes.titlesize": fontsizes["MEDIUM_SIZE"],
        "axes.labelsize": fontsizes["SMALL_SIZE"],
        "xtick.labelsize": fontsizes["SMALL_SIZE"],
        "ytick.labelsize": fontsizes["SMALL_SIZE"],
        "legend.fontsize": fontsizes["MEDIUM_SIZE"],
        "font.size": fontsizes["SMALL_SIZE"],
    }

    pgf_style = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
    }

    if os.environ.get("MPLBACK", None) == "pgf":
        print("Plotting backend seems to be PGF, adding adequate keywords to rcParams")
        default_style.update(pgf_style)

    if figsize is None:
        # figsize = (6.4, 4.8)
        figsize = (7.2, 5.4)
    default_style["figure.figsize"] = figsize

    if remove_spines:
        default_style["axes.spines.top"] = False
        default_style["axes.spines.right"] = False

    if style_kws is None:
        style_kws = dict()
    default_style.update(style_kws)

    plt.rcParams.update(default_style)
    return style_kws


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
