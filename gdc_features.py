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

gdc_any_tumor = [
    "Primary Tumor",
    "Human Tumor Original Cells",
    "Additional Metastatic",
    "Additional - New Primary",
    "Primary Blood Derived Cancer - Bone Marrow",
    "Recurrent Tumor",
    "Metastatic",
    "Recurrent Blood Derived Cancer - Peripheral Blood",
    "Recurrent Blood Derived Cancer - Bone Marrow",
    "Primary Blood Derived Cancer - Peripheral Blood",
]

tcga_cohort_diagnosis = {
    "TCGA-ACC": "Adrenocortical cc.",
    "TCGA-BLCA": "Bladder cc.",
    "TCGA-BRCA": "Breast cancer",
    "TCGA-CESC": "Cervical cancer",
    "TCGA-CHOL": "Cholangiocc.",
    "TCGA-COAD": "Colon cancer",
    "TCGA-DLBC": "Lymphoma",
    "TCGA-ESCA": "Esophageal cc.",
    "TCGA-GBM": "Glioblastoma multiforme",
    "TCGA-HNSC": "Head and neck cc.",
    "TCGA-KICH": "Kidney chromophobe cc.",
    "TCGA-KIRC": "Kidney clear cell cc.",
    "TCGA-KIRP": "Kidney papillary cc.",
    "TCGA-LGG": "Lower grade glioma",
    "TCGA-LIHC": "Hepatocellular cc.",
    "TCGA-LUAD": "Lung adenocc.",
    "TCGA-LUSC": "Lung squamous cell cc.",
    "TCGA-MESO": "Mesothelioma",
    "TCGA-OV": "Ovarian serous cc.",
    "TCGA-PAAD": "Pancreatic adenocc.",
    "TCGA-PCPG": "Pheochromocytoma",
    "TCGA-PRAD": "Prostate adenocc.",
    "TCGA-READ": "Rectal cancer",
    "TCGA-SARC": "Sarcoma",
    "TCGA-SKCM": "Cutaneous melanoma",
    "TCGA-STAD": "Gastric cancer",
    "TCGA-TGCT": "Testicular cancer",
    "TCGA-THYM": "Thymoma",
    "TCGA-THCA": "Thyroid cc.",
    "TCGA-UCS": "Uterine carcinosarcoma",
    "TCGA-UCEC": "Endometrial cc.",
    "TCGA-UVM": "Uveal melanoma",
}


def recipe(*, verbose: bool = True,) -> dict:

    """
    Common definitions for GDC terminology.

    Parameters
    ----------
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    Definitions supplied with this package
    
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
    hint(verbose, "Definitions supplied with this package:\n", example_variables)

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
