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

xena_hub = "https://gdc.xenahubs.net"
dataset = "TCGA-BRCA.GDC_phenotype.tsv"
gexdata = "TCGA-BRCA.htseq_fpkm-uq.tsv"
phenotypetag = ".GDC_phenotype.tsv"
gextag = ".htseq_fpkm-uq.tsv"
probemap = "gencode.v22.annotation.gene.probeMap"
cohort = "TCGA-BRCA"
lung_cohorts = ["TCGA-LUSC", "TCGA-LUAD"]
survcol = "vital_status.demographic"
target_gene = "TP53"
target_genes = [target_gene]
prognosis_genes = ["p53", "KI-67", "Estrogen receptor", "Progesterone receptor"]
prognosis_genecodes = {
    "p53": "TP53",
    "KI-67": "MKI67",
    "Estrogen receptor": "ESR1",
    "Progesterone receptor": "PGR",
}
ptable = {
    "cohort": ["TCGA-LIHC", "TCGA-BRCA"],
    "TP53": [0.48725, 0.435489],
    "MKI67": [0.1292, 3.375975],
}
curated_survivals = "clinical_endpoints.tsv"
pancan_sampletypes = "PANCAN_xena_sample_types.tsv"
gdrive_secret = os.path.expanduser("~") + "/.key"
quadKMlabels = [
    "WT low",
    "WT high",
    "mut low",
    "mut high",
    "low",
    "high",
]
quadKMcolors = ("#2ca02c", "#ff7f0e", "#1f77b4", "#d62728", "#1b9e77", "#d95f02")
hazardColors = {"WT": "#1b9e77", "mut": "#d95f02"}


def recipe(*, verbose: bool = True,) -> dict:

    """
    Example parameters that can be used as default for testing purposes.

    Parameters
    ----------
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    Examples supplied with this package
    
    """

    ### List variables that can be used to run tests
    all_def = globals()
    not_an_example = [
        "__name__",
        "__doc__",
        "__package__",
        "__loader__",
        "__spec__",
        "__file__",
        "__cached__",
        "__builtins__",
        "__annotations__",
        "os",
        "sys",
        "introSpect",
        "invoked_directly",
        "hint",
        "main",
        "recipe",
    ]
    example_variables = {
        key: value for key, value in all_def.items() if key not in not_an_example
    }
    hint(verbose, "Examples supplied with this package:\n", example_variables)

    return example_variables


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
