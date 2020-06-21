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


sampleFile = basedir+'/sample_types.tsv'
survivalFile = basedir+'/clinical_endpoints.tsv'

# Read clinical and sample data
clinicals = pd.read_csv(survivalFile, sep='\t')
clinicals = clinicals.set_index('bcr_patient_barcode')
sampleTypes = pd.read_csv(sampleFile, sep='\t')
sampleTypes = sampleTypes[~pd.isnull(sampleTypes['sample'])]
sampleTypes['patient'] = sampleTypes['sample'].apply(lambda x: '-'.join(x.split('-')[:3]))
sampleTypes = sampleTypes.set_index('sample')

rightSample = [
    'Primary Tumor',
    'Human Tumor Original Cells',
    'Additional Metastatic',
    'Additional - New Primary',
    'Primary Blood Derived Cancer - Bone Marrow',
    'Recurrent Tumor',
    'Metastatic',
    'Recurrent Blood Derived Cancer - Peripheral Blood',
    'Recurrent Blood Derived Cancer - Bone Marrow',
    'Primary Blood Derived Cancer - Peripheral Blood'
]

# Create a dictionary of gene symbols

response = requests.get(probeMapLink)
f = io.StringIO(response.text)
genedict = pd.read_csv(f, sep='\t')
genedict = genedict.set_index('id')
genedict = genedict.loc[:, 'gene'].to_dict()


# Retrieve a list of all genes

geneslice = 500
fl_allgenes = set()
for prf in prefixes:
    g = xena.dataset_field_examples(xena_hub, prf+gex_dataset, xena.dataset_field_n(xena_hub, prf+gex_dataset))
    fl_allgenes.update(g)
fl_allgenes = list(fl_allgenes)
allgenes = np.array(fl_allgenes)
N_genes = allgenes.shape[0]
gsn = int(allgenes.shape[0]/geneslice)
rest = allgenes[geneslice*gsn:]
allgenes = allgenes[:geneslice*gsn].reshape(-1, geneslice).tolist()
allgenes.append(rest.tolist())

# Get a list of samples
patients = clinicals.loc[:, ['OS', 'OS.time']]
sampleTab = sampleTypes.loc[sampleTypes['patient'].isin(patients.index.values) & sampleTypes['sample_type'].isin(rightSample),:]
samples = sampleTab.index.values.tolist()

# Add PIDD expression category
for prf in prefixes:
    position, expression_matrix = xena.dataset_probe_values(xena_hub, prf+gex_dataset, samples, basegenes)
    expression_matrix = np.array(expression_matrix)
    sampleTab['PIDD_gex'] = expression_matrix[0]

sampleTab = sampleTab.loc[sampleTab['PIDD_gex'] != 'NaN', :]
sampleTab['PIDD_gex'] = pd.to_numeric(sampleTab['PIDD_gex'])
p_mn = sampleTab['PIDD_gex'].median()
sampleTab['PIDD_cat'] = sampleTab['PIDD_gex'].apply(lambda x: 'low' if x < p_mn else'high')
sampleTab['survival'] = patients.loc[sampleTab['patient'],['OS']].values
sampleTab['time'] = patients.loc[sampleTab['patient'],['OS.time']].values
samples = sampleTab.index.values.tolist()
sampleTab.head()

# Assess interference of genes on PIDD survival impact
pRows = []
for prf in prefixes:
    for chunk in allgenes:
        position, expression_matrix = xena.dataset_probe_values(xena_hub, prf+gex_dataset, samples, chunk)
        expression_matrix = np.array(expression_matrix)
        for i in range(len(chunk)):
            tdf = sampleTab.copy()
            tdf['gex'] = expression_matrix[i]
            tdf = tdf.loc[tdf['gex'] != 'NaN', :]
            tdf['gex'] = pd.to_numeric(tdf['gex'])
            t_mn = tdf['gex'].median()
            tdf['gex'] = tdf['gex'].apply(lambda x: 'low' if x < t_mn else'high')
            pmask = (tdf['PIDD_cat'] == 'low')
            gmask = (tdf['gex'] == 'low')
            T = tdf['time']
            E = tdf['survival']
            s1 = statistics.logrank_test(T[pmask & gmask], T[pmask & ~gmask], event_observed_A=E[pmask & gmask], event_observed_B=E[pmask & ~gmask])
            s2 = statistics.logrank_test(T[~pmask & gmask], T[~pmask & ~gmask], event_observed_A=E[~pmask & gmask], event_observed_B=E[~pmask & ~gmask])
            pRows.append((chunk[i], np.log(s1.p_value), np.log(s2.p_value)))
    pRows = pd.DataFrame.from_records(pRows)
    pRows = pRows.fillna(0)
    pRows.columns = ['gene', 'logP1', 'logP2']
    pRows['gain'] = pRows['logP1']-pRows['logP2']
    pRows = pRows.sort_values(by='gain') 
    pRows.reset_index().to_csv(prf+'_PIDD1_synleth.tsv', sep='\t', index=False)

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
