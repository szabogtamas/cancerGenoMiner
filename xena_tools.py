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

import os, io, requests
import numpy as np
import pandas as pd
import xenaPython as xena
from typing import Union, List
from . import par_examples


def recipe(
    *,
    xena_hub: str = par_examples.xena_hub,
    dataset: str = par_examples.dataset,
    probemap: str = par_examples.probemap,
    survcol: str = par_examples.survcol,
    verbose: bool = True,
) -> pd.DataFrame:

    """
    Shortcut functions to manage data access with xenaPython.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset we want to check on.
    probemap
        A probemap file for testing direct table download.
    survcol
        Column containing the critical event (survival).
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    A table with phenotype data for the given dataset
    """

    ### Check what datasets (cohorts) are available on the hub
    cohorts = check_cohort_names(xena_hub)
    hint(verbose, "Available cohorts on", xena_hub, ":\n", cohorts)

    ### Check what column names we have in a given dataset
    metanames = check_metafields(xena_hub, dataset)
    hint(verbose, "Available fields in", dataset, ":\n", metanames)

    ### Check sample names in the dataset
    samples = check_allsamples(xena_hub, dataset)
    hint(verbose, "First 5 samples in", dataset, ":\n", samples[:5])

    ### Read a probe mapping directly from the hub
    probes = read_xena_table(probemap, hubPrefix=xena_hub).head()
    hint(verbose, "Snapshot of", probemap, ":\n", probes)

    ### Read a clinical phenotype table
    clinicals = download_gdc_clinicals(xena_hub, dataset)
    hint(verbose, "Snapshot of", dataset, ":\n", clinicals.head())

    ### Changes factor levels in survival column to event codes
    clinicals = fix_phenotype_factorlevels(
        clinicals,
        xena_hub,
        dataset,
        survcol,
        leveldict={"NaN": 0},
        renamedict={"Dead": 1},
        renameremain=0,
    )
    hint(verbose, dataset, "after categories fixed:\n", clinicals.head())

    return clinicals


def check_cohort_names(
    xena_hub: str, *, blacklist: Union[None, List] = None, TCGAonly: bool = True,
) -> List:

    """
    List what cohorts are available at a given hub.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    blacklist
        A list of cohort names we want to exclude.
    TCGAonly
        Sets if cohort names should start with "TCGA". True by default.

    Returns
    -------
    List of cohort names.
    """

    prefixes = []
    if blacklist is None:
        blacklist = []
    for x in xena.all_datasets(xena_hub):
        x = x["name"].split(".")[0]
        if x[:4] == "TCGA" or not TCGAonly:
            if x not in blacklist:
                prefixes.append(x)
    prefixes = list(set(prefixes))
    return prefixes


def check_allsamples(
    xena_hub: str, dataset: str, *, n: Union[None, int] = None,
) -> List:

    """
    Check what samples are available in the dataset.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset we want to check on.
    n
        Limit the number of samples returned. Set to None if all samples are needed.

    Returns
    -------
    A list of all sample names.
    """

    return xena.dataset_samples(xena_hub, dataset, n)


def check_metafields(xena_hub: str, dataset: str,) -> List:

    """
    Show what phenotype information is available.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset we want to check on.

    Returns
    -------
    List of available phenotype information columns.
    """

    metafields = []
    for m in xena.feature_list(xena_hub, dataset):
        metafields.append(m["name"])
    return metafields


def read_xena_table(
    linkOrFile: str, *, hubPrefix: Union[None, str] = None,
) -> pd.DataFrame:

    """
    Read clinical information or probe names, preferentially directly form Xena.

    Parameters
    ----------
    linkOrFile
        The fastq file resulted from sequencing of the samples.
    hubPrefix
        A prefix for the link, typically the Xena hub url.

    Returns
    -------
    Data table with clinical information or gene and probe names.
    """

    if hubPrefix not in [None, ""]:
        linkOrFile = hubPrefix + "/download/" + linkOrFile
    if os.path.isfile(linkOrFile):
        f = linkOrFile
    else:
        response = requests.get(linkOrFile)
        f = io.StringIO(response.text)
    genemap = pd.read_csv(f, sep="\t")
    return genemap


def download_gdc_clinicals(
    xena_hub: str,
    dataset: str,
    *,
    rowKey: str = "SampleID",
    headers: List = [
        "sample_type.samples",
        "days_to_death.demographic",
        "days_to_last_follow_up.diagnoses",
        "vital_status.demographic",
    ],
) -> pd.DataFrame:

    """
    Clinical information contains many factor columns, where factor names are not
    readily accessible. This function decodes factor numbers into category names.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset containing survival information on the repository hub.
    rowKey
        The column containing (unique) identifiers - typically the sample IDs.
    headers
        A list of column names that we want to retrieve.

    Returns
    -------
    Dataframe with clinical information, including survival.
    """

    samples = xena.dataset_samples(xena_hub, dataset, None)
    pos, mat = xena.dataset_probe_values(xena_hub, dataset, samples, headers)
    clinicals = dict()
    clinicals[rowKey] = samples
    for i in range(len(mat)):
        clinicals[headers[i]] = mat[i]
    clinicals = pd.DataFrame(clinicals)
    clinicals = clinicals.set_index(rowKey)
    return clinicals


def fix_phenotype_factorlevels(
    clinicals: pd.DataFrame,
    xena_hub: str,
    dataset: str,
    targetcol: str,
    *,
    leveldict: dict = {"NaN": 0},
    renamedict: Union[None, dict] = None,
    renameremain: Union[None, str, int, float] = None,
) -> pd.DataFrame:

    """
    Clinical information contains many factor columns, where factor names are not
    readily accessible. This function decodes factor numbers into category names.

    Parameters
    ----------
    clinicals
        A dataframe where the column is to be fixed.
    xena_hub
        Url of the data repository hub.
    ds
        Name of the dataset on the repository hub.
    targetcol
        The column we have to fix.
    leveldict
        A dictionary to override some level naming if needed. Must have a NaN key.
    renamedict
        A dictionary to set an alternative name for categories (e.g. 1 for Deceased).
    renameremain
        If the category name is not in the rename dictionary, it keeps its name
        by default. If renameremain is not None, then remaning items get this name.

    Returns
    -------
    Dataframe with factor level names instead of level numbers.
    """

    codes = xena.field_codes(xena_hub, dataset, [targetcol])
    codes = codes[0]
    codes = codes["code"].split("\t")
    if renamedict is None:
        renamedict = dict()
    for i in range(len(codes)):
        e = codes[i]
        if e in renamedict:
            e = renamedict[e]
        else:
            if renameremain is not None:
                e = renameremain
        if i not in leveldict:
            leveldict[i] = e
    clinicals[targetcol] = clinicals[targetcol].map(leveldict)
    clinicals[targetcol] = clinicals[targetcol].map(renamedict)
    return clinicals


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
