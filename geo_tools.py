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

import urllib, re
import pandas as pd
from io import StringIO
from bs4 import BeautifulSoup
from typing import Union, Tuple, List, Sequence
from . import par_examples, xena_tools


def recipe(
    *,
    dataset: str = "GDS4053",
    genes_of_interest: str = par_examples.prognosis_genecodes,
    verbose: bool = True,
) -> pd.DataFrame:

    """
    Functions facilitating the interaction with Gene Expression Ominbus.

    Parameters
    ----------
    dataset
        Dataset GEO ID.
    genes_of_interest
        Gene symbols to be submitted.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    Long format data table featuring gene expression for genes of interest, in each
    sample of the dataset.
    """

    ### Download gene expression from GEO
    gex = agex_from_GEO_dataset(dataset, genes_of_interest)
    hint(
        verbose, "Gene expression in", dataset, gex.head(),
    )

    ### Explore grouping of samples
    gex["GroupName"] = gex.apply(explore_sample_grouping, axis=1)
    hint(
        verbose, "Groups found:", gex["GroupName"].unique(),
    )

    return gex


def gex_from_GEO_dataset(dataset: str, genes_of_interest: list) -> pd.DataFrame:

    """
    Retrieve gene expression data from GEO for a given dataset.

    Parameters
    ----------
    dataset
        Dataset GEO ID.
    genes_of_interest
        Gene symbols to be submitted.

    Returns
    -------
    A data table with gene expression for all genes in all samples of the dataset.
    """

    extracted_data = []
    for genesym in genes_of_interest:
        refs = resolve_genesym_to_RefId(dataset, genesym)
        for gene_ref in refs:
            url = (
                "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?ID="
                + dataset
                + ":"
                + gene_ref
            )
            try:
                contents = urllib.request.urlopen(url).read()
                contents = contents.decode("utf-8")
                parsed_html = BeautifulSoup(contents, "html.parser")
                title_data = parsed_html.body.find("table", attrs={"id": "titleTable"})
                tdd = []
                rows = title_data.find_all("tr")
                for row in rows:
                    tdd.append([col.text for col in row.find_all("td")])
                t = tdd[1][1]
                table_data = parsed_html.body.find("table", attrs={"id": "sampleTable"})
                extract = ""
                rows = table_data.find_all("tr")
                for row in rows:
                    cols = row.find_all("td")
                    temporary = [dataset, t, genesym, gene_ref]
                    for col in cols:
                        temporary.append(col.text)
                    if len(temporary) == 8:
                        extracted_data.append(temporary)
            except:
                pass
    extracted_data = pd.DataFrame(
        extracted_data,
        columns=[
            "Dataset",
            "Description",
            "Gene",
            "ProbeID",
            "SampleID",
            "SampleName",
            "GEX",
            "PercentileRank",
        ],
    )
    return extracted_data


def resolve_genesym_to_RefId(dataset: str, genesym: str) -> list:

    """
    Maps between probe names (IDs) and gene symbols.

    Parameters
    ----------
    dataset
        Dataset GEO ID.
    genesym
        Gene symbol for gene to be mapped.

    Returns
    -------
    A list of probes for a given gene.
    """

    refs = []
    url = (
        "https://www.ncbi.nlm.nih.gov/geoprofiles/?term="
        + dataset
        + "%5BACCN%5D+"
        + genesym
        + "%5BGene+Symbol%5D&report=docsum&format=text"
    )
    with urllib.request.urlopen(url) as site:
        contents = site.read()
        contents = contents.decode("utf-8")
        for line in contents.split("\n"):
            if line.find("Reporter: ") == 0:
                for i in line[10:].split(", "):
                    if i.find(" (ID_REF)") > -1:
                        refs.append(i.replace(" (ID_REF)", ""))
    return refs


def explore_sample_grouping(s: pd.Series) -> str:

    """
    Creates group names by stripping numbers off replicate samples. To be used with
    `df.apply`.

    Parameters
    ----------
    s
        Rows of GEO data table as a Pandas series.

    Returns
    -------
    Pandas data series with group name for each record.
    """

    REPLICATES_REGEX = re.compile(
        "(biological replicate|biological rep|replicate|rep|Rep)\s*\d*", re.S
    )
    ENDNUMBERS_REGEX = re.compile("\d*$", re.S)
    group = str(s["SampleName"])
    truncated = False
    match = REPLICATES_REGEX.search(group)
    if match:
        groupname_end = match.start()
        group = group[:groupname_end].strip()
        truncated = True
    if not truncated:
        match = ENDNUMBERS_REGEX.search(group)
        if match:
            groupname_end = match.start()
            group = group[:groupname_end].strip()
    return group


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
