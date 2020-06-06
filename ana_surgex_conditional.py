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
nextflowProcess = introSpect.flowNodes.nextflowProcess

import scipy
import numpy as np
from typing import Union, Tuple
from cancerGenoMiner import (
    environment_definiton,
    par_examples,
    gdc_features,
    plotting_tools,
    survival_tools,
    xena_tools,
    gex_tools,
    reporting_tools,
    ana_surgex_single,
)
import connectDrive

cgm_folder = os.path.dirname(os.path.realpath(gex_tools.__file__))
intro_folder = os.path.dirname(os.path.realpath(introSpect.__file__))
drive_folder = os.path.dirname(os.path.realpath(connectDrive.__file__))


def create_pipeline(
    *,
    location: str = os.getcwd(),
    nodes: Union[None, str] = None,
    main_kws: Union[None, dict] = None,
    comment_on_methods: Union[None, str] = None,
    conda: Union[None, str] = None,
    general_configs: Union[None, dict] = None,
    container_paths: Union[None, dict] = None,
    verbose: bool = True,
) -> str:
    """
    Create a Nextflow pipeline that take a list of genes and a list of cohorts, compares
    survival in high and low expression groups conditionally, creating four groups based
    on gene axpression and an additional factor. The overview heatmap will show the
    gained logP when applying the additional factor.

    Parameters
    ----------
    location
        Path where the pipeline directory should be created.
    nodes
        Objects that define processes as nodes linked by Nextflow.
    main_kws
        Initial pipeline parameters.
    comment_on_methods
        Description of methods.
    conda
        Path to a yaml file to be used during environment creation or the conda dir.
    general_configs
        Settings in the nextflow config file that apply to the whole pipeline in general.
    container_paths
        Paths to containers defined by their Dockerhub link. Gets automatically updated
        if an image is built, but adding a path here avoids rebuilding existing images.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    Body of the nextflow script
    """

    ### Define the main parameters for the Nextflow script
    comment_location = location + "/pipeline/comments.tex"
    default_main_kws = {
        "cohorts": [par_examples.cohort] + par_examples.lung_cohorts,
        "genes": par_examples.prognosis_genes,
        "genedict": {"NaN": "NaN"},  # par_examples.prognosis_genecodes,
        "survtab": "None",
        "report_title": '"Title of report"',
        "author_name": '"Author"',
        "lab_name": '"Laboratory or PI"',
        "bibliography": "biblio.bib",
        "slidesfile": "figures.pptx",
        "comments": comment_location,
        "mpl_backend": "pgf",
        "fn": "survival_table.tsv",
    }

    if main_kws is None:
        main_kws = dict()
    default_main_kws.update(main_kws)
    main_kws = default_main_kws.copy()

    if container_paths is None:
        container_paths = dict()

    ### Set the Conda environment
    if conda is None:
        os.makedirs(location + "/pipeline", exist_ok=True)
        conda = location + "/pipeline/environment.yml"
        with open(conda, "w") as f:
            f.write(environment_definiton.environment)
        conda = "'" + conda + "'"

    ### Add processes to the Nextflow pipeline
    if nodes is None:
        nodes = [
            ana_surgex_single.fetchClinicalFile(
                inchannels=["survtab"], outchannels=["local_survtab"], conda=conda
            ),
            ana_surgex_single.getSurvival(
                inchannels=["cohorts", "genes", "genedict", "nicer_survtab"],
                outchannels=["survivals", "gene_nd"],
                conda=conda,
            ),
            ana_surgex_single.plotSurvival(
                inchannels=[
                    "survivals",
                    "plotgenes",
                    "symdict",
                    "numplot",
                    "mpl_backend",
                ],
                outchannels=[
                    "plotnames",
                    "gexnames",
                    "images",
                    "plots",
                    "stats",
                    "titles",
                    "notebooks",
                ],
                conda=conda,
                capture=True,
            ),
            ana_surgex_single.makeHeatmap(
                inchannels=["stat", "heatmap_f"],
                outchannels=["heatmap", "heatimage", "cohort_order"],
                conda=conda,
            ),
            ana_surgex_single.pptFromFigures(conda=conda),
            ana_surgex_single.getRefs(
                inchannels=["bibliography", "drive_key"],
                outchannels=["bibtex"],
                conda=conda,
            ),
            ana_surgex_single.compileReport(
                inchannels=[
                    "plotsr",
                    "heatmap",
                    "ordered_cohorts",
                    "page_titles",
                    "comments",
                    "bibtex",
                    "report_title",
                    "author_name",
                    "lab_name",
                ],
                outchannels=["reportex"],
                conda=conda,
            ),
            ana_surgex_single.pdfFromLatex(),
        ]

    ### Compile the pipeline into a temporary folder
    introSpect.flowNodes.channelNodes(
        *nodes,
        main_kws=main_kws,
        location=location + "/pipeline",
        generalSettings=general_configs,
        containerPaths=container_paths,
        verbose=verbose,
    )
    hint(verbose, "Pipeline compiled to folder")

    ### Save comments to a file
    if comment_on_methods is None:
        comment_on_methods = ""
    with open(comment_location, "w") as f:
        f.write(comment_on_methods)
    return os.path.dirname(location + "/pipeline/main.nf")

recipe = create_pipeline

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
    }
    mainFunction = introSpect.commandLines.cmdConnect(recipe, modified_kws)
    mainFunction.eval()
    mainFunction.save()
    return


if invoked_directly:
    main()
