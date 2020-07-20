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

import pandas as pd
import numpy as np
import xenaPython as xena
from typing import Union, Tuple, List, Sequence
from . import par_examples, xena_tools


def recipe(
    *,
    xena_hub: str = par_examples.xena_hub,
    dataset: str = par_examples.gexdata,
    probemap: str = par_examples.probemap,
    genes: str = par_examples.prognosis_genes,
    gene_symbols: str = par_examples.prognosis_genecodes,
    verbose: bool = True,
) -> pd.DataFrame:

    """
    Functions facilitating the interaction with UCSC Xena datahub to download gene
    expression.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset we want to check on.
    probemap
        A probemap file for testing direct table download.
    genes
        A list of genes to be queried.
    gene_symbols
        Gene symbols of the query genes.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    Gene expression matrix for the supplied genes and samples.
    
    """

    ### Retrieve the mapping between genes and probes (ENSEMBL IDs)
    probes = xena_tools.read_xena_table(probemap, hubPrefix=xena_hub)
    probedict = parse_gene_mapping(probes)
    for g in genes:
        if g in gene_symbols:
            n = gene_symbols[g]
        else:
            n = g
        if n in probedict:
            hint(verbose, g, n, probedict[n])
        else:
            hint(verbose, g, n, "not found")

    ### Test if every gene could be mapped correctly
    unmapped = check_probedict(genes, probedict, gene_symbols=gene_symbols)
    if len(unmapped) < 1:
        hint(verbose, "All genes could be mapped to probes.")
    else:
        hint(verbose, "Following genes could not be mapped:\n", unmapped)

    ### Search for a specific ENSEMBL ID in the mapping
    hint(verbose, probes.loc[probes["id"].str.contains("ENSG00000082175"), :])

    ### Retrieve a list of samples
    samples = xena_tools.check_allsamples(xena_hub, dataset)
    hint(verbose, "First 5 samples in", dataset, ":\n", samples[:5])

    ### Convert sample list to dataframe
    df = pd.DataFrame(samples, index=samples, columns=["SampleID"])
    hint(verbose, df.head())

    ### Get gene expression by genes
    gex = add_gene_expression_by_genes(
        genes, df, xena_hub, dataset, gene_names=gene_symbols
    )
    hint(
        verbose,
        "Gene expression in",
        dataset,
        "retrived directly on genes:\n",
        gex.head(),
    )

    ### Get gene expression by probes
    gex = add_gene_expression_by_probes(
        genes, df, xena_hub, dataset, probedict, gene_names=gene_symbols
    )
    hint(
        verbose,
        "Gene expression in",
        dataset,
        "retrived after mapping to probes:\n",
        gex.head(),
    )

    ### Mark samples where expression is in the lower or upper quartile
    for gene in genes:
        gex = split_by_gex_quartile(gex, gene)
    hint(verbose, "Gene expression categories based on quartiles:\n", gex.head())

    ### Split samples by median expression
    for gene in genes:
        gex = split_by_gex_median(gex, gene)
    hint(verbose, "Gene expression categories based on the median:\n", gex.head())

    return gex


def split_by_gex_median(
    gex: pd.DataFrame,
    gene: str,
    *,
    colname: str = "gex_",
    catname: str = "cat_",
    labels: Tuple[str, str] = ("low", "high"),
) -> pd.DataFrame:

    """
    Categorize gene expression into low (less than median) and high (median or greater).

    Parameters
    ----------
    gex
        A dataframe containing gene expression values.
    gene
        The gene whose expression we are looking for.
    colname
        The column with values we want to discretize.
    catname
        Name of the categorical column this function creates.
    labels
        Category names for low and high expression.

    Returns
    -------
    Dataframe with gene expression discretized into low and high categories.
    """

    colname += gene
    catname += gene
    gex = gex.loc[~pd.isnull(gex[colname]), :]
    gex[colname] = gex[colname].astype(float)
    gexmed = gex[colname].median()
    low, high = labels
    gex[catname] = gex[colname].apply(lambda x: low if x < gexmed else high)
    return gex


def split_by_gex_quartile(
    gex: pd.DataFrame,
    gene: str,
    *,
    colname: str = "gex_",
    catname: str = "cat_",
    labels: Tuple[str, str] = ("low", "high"),
) -> pd.DataFrame:

    """
    Categorize gene expression into low (1st quartile) and high (4th quartile).

    Parameters
    ----------
    gex
        A dataframe containing gene expression values.
    gene
        The gene whose expression we are looking for.
    colname
        The column with values we want to discretize.
    catname
        Name of the categorical column this function creates.
    labels
        Category names for low and high expression.

    Returns
    -------
    Dataframe with gene expression discretized into low and high categories.
    """

    colname += gene
    catname += gene
    gex = gex.loc[~pd.isnull(gex[colname]), :]
    gex[colname] = gex[colname].astype(float)
    lowquart, highquart = gex[colname].quantile([0.25, 0.75]).tolist()
    low, high = labels
    gex[catname] = gex[colname].apply(lambda x: low if x <= lowquart else "")
    gex[catname] = gex[catname].astype(str) + gex[colname].apply(
        lambda x: high if x >= highquart else ""
    )
    return gex


def add_gene_expression_by_probes(
    target_genes: Sequence,
    clinicals: pd.DataFrame,
    xena_hub: str,
    ds: str,
    *,
    probemap: Union[None, dict] = None,
    gene_names: Union[None, dict] = None,
    colprefix: str = "gex_",
) -> pd.DataFrame:

    """
    Downloads gene expression of a set of genes for all samples in the cohort.

    Parameters
    ----------
    target_genes
        A list of genes whose expression we are interested in.
    clinicals
        A dataframe containing sample information. A column for each gene will be added.
    xena_hub
        Url of the data repository hub.
    ds
        Name of the dataset on the repository hub.
    probemap
        A dictionary to convert gene symbols to ENSEBML IDs.
    gene_names
        If gene names are not identical to gene symbols (or the genome version is
        different), this maps genes to their expression data.
    colprefix
        Added to the gene name to label the column with expression data.

    Returns
    -------
    Dataframe with gene expression data added to the previous content.
    """

    target_probes = []
    if gene_names is None:
        gene_names = dict()
    if probemap is None:
        probemap = dict()
    for g in target_genes:
        if g in gene_names:
            g = gene_names[g]
        if g in probemap:
            g = probemap[g]
        target_probes.append(g)
    position, expression_matrix = xena.dataset_probe_values(
        xena_hub, ds, clinicals.index.values.tolist(), target_probes
    )
    for i in range(len(target_genes)):
        colname = colprefix + target_genes[i]
        clinicals[colname] = expression_matrix[i]
    return clinicals


def add_gene_expression_by_genes(
    target_genes: Sequence,
    clinicals: pd.DataFrame,
    xena_hub: str,
    ds: str,
    *,
    gene_names: Union[None, dict] = None,
    colprefix: str = "gex_",
) -> pd.DataFrame:

    """
    Downloads gene expression of a set of genes for all samples in the cohort.

    Parameters
    ----------
    target_genes
        A list of genes whose expression we are interested in.
    clinicals
        A dataframe containing sample information. A column for each gene will be added.
    xena_hub
        Url of the data repository hub.
    ds
        Name of the dataset on the repository hub.
    gene_names
        If gene names are not identical to gene symbols (or the genome version is
        different), this maps genes to their expression data.
    colprefix
        Added to the gene name to label the column with expression data.

    Returns
    -------
    Dataframe with gene expression data added to the previous content.
    """

    genes = []
    if gene_names is None:
        gene_names = dict()
    for g in target_genes:
        if g in gene_names:
            g = gene_names[g]
        genes.append(g)

    expression_matrix = xena.dataset_gene_probe_avg(
        xena_hub, ds, clinicals.index.values.tolist(), genes
    )
    for i in range(len(target_genes)):
        colname = colprefix + target_genes[i]
        if len(expression_matrix[i]["scores"][0]) < 1:
            print(
                colname,
                "not found. Are you sure you provided the gene symbol corresponding to the genome version?",
            )
            clinicals[colname] = "NaN"
        else:
            clinicals[colname] = expression_matrix[i]["scores"][0]
    return clinicals


def create_gene_chunks(
    cohort: str,
    *,
    xena_hub: str = par_examples.xena_hub,
    gex_prefix: str = par_examples.gextag,
    chunk_size: int = 500,
) -> list:

    """
    Retrieve all genes in the gex datasets in chunks.

    Parameters
    ----------
    cohort
        The TCGA cohort to check.
    xena_hub
        Url of the data repository hub.
    gex_prefix
        Constant part of the gene expression dataset name.
    chunk_size
        Number of genes to be grouped together.

    Returns
    -------
    List of evenly sized gene lists.
    """

    allgenes = np.array(
        xena_tools.xena.dataset_field_examples(xena_hub, cohort + gex_prefix, None)
    )
    N_genes = allgenes.shape[0]
    gsn = int(allgenes.shape[0] / chunk_size)
    rest = allgenes[chunk_size * gsn :]
    rest = rest.tolist()
    if rest[-1] == "sampleID":
        rest.pop()
    geneslices = allgenes[: chunk_size * gsn].reshape(-1, chunk_size).tolist()
    geneslices.append(rest)
    return geneslices


def check_probedict(
    genes: Sequence, genemap: str, *, gene_symbols: Union[None, dict] = None,
) -> dict:

    """
    Check if genes could be mapped to probes. List unmapped genes.

    Parameters
    ----------
    genes
        List of genes we want to map.
    genemap
        The dictionary created form the gene-to-probe mapping.
    gene_symbols
        If gene names are not identical to gene symbols (or the genome version is
        different), this maps genes to their expression data.

    Returns
    -------
    A list of genes that could not be mapped to probes.
    """

    if gene_symbols is None:
        gene_symbols = dict()

    unmapped = []
    for g in genes:
        if g in gene_symbols:
            n = gene_symbols[g]
        else:
            n = g
        if n not in genemap:
            unmapped.append(n)

    return unmapped


def parse_gene_mapping(
    genemap: str, *, probecol: str = "id", genecol: str = "gene",
) -> dict:

    """
    Read clinical information or probe names, preferentially directly form Xena.

    Parameters
    ----------
    genemap
        The table created form the gene map file specified in the dataset description.
    probecol
        The column with probe names.
    genecol
        The column with gene names.

    Returns
    -------
    A mapping from genes to probe names.
    """

    genedict = genemap.set_index(genecol)
    genedict = genedict.loc[:, probecol].to_dict()
    return genedict


def map_genenames(genes: str, genedictfile: str, *, idx: str = "name",) -> list:

    """
    Maps between gene names and gene symbols.

    Parameters
    ----------
    genes
        List of genes to be mapped.
    genedictfile
        Path to the mapping table.
    idx
        Index column (name or symbol).

    Returns
    -------
    A mapped list.
    """

    if genedictfile is None:
        return genes
    genedict = pd.read_csv(genedictfile, sep="\t", names=["name", "symbol"])
    genedict = genedict.set_index(idx)
    genedict = genedict.T.to_dict("records")
    genedict = genedict[0]
    ngenes = []
    for gene in genes:
        gene = gene.replace('"', "")
        if gene in genedict:
            ngenes.append(genedict[gene])
        else:
            ngenes.append(gene)
    return ngenes


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
