#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)  # developmental hack, remove later!
import introSpect
import connectDrive

__name__, __package__, invoked_directly = introSpect.cmdSupport(
    __name__, __package__, __file__
)
hint = introSpect.hint

import pandas as pd
from lifelines import statistics, KaplanMeierFitter
from typing import Union, Tuple, List, Sequence
from . import par_examples, gdc_features, xena_tools, gex_tools, plotting_tools
plt = plotting_tools.plt

def recipe(
    *,
    xena_hub: str = par_examples.xena_hub,
    dataset: str = par_examples.dataset,
    gex_dataset: str = par_examples.gexdata,
    probemap: str = par_examples.probemap,
    gene: str = par_examples.target_gene,
    survcol: str = par_examples.survcol,
    survfile: str = par_examples.curated_survivals,
    samplefile: str = par_examples.pancan_sampletypes,
    gdrive_secret: str = par_examples.gdrive_secret,
    verbose: bool = True,
) -> pd.DataFrame:

    """
    Convenience functions to plot and compare Kaplan–Meier curves of individuals split
    up in two groups.

    Parameters
    ----------
    xena_hub
        Url of the data repository hub.
    dataset
        The dataset we want to check on.
    gex_dataset
        Dataset containing gene expression information.
    gene
        The gene to be queried.
    survcol
        Column containing the critical event (survival).
    survfile
        Data table on Google Drive containig curated survival information.
    gdrive_secret
        Location of the credential file for Google Drive access.
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    The Kaplan–Meier plot.
    
    """

    ### Read a clinical phenotype table
    clinicals = xena_tools.download_gdc_clinicals(xena_hub, dataset)
    hint(verbose, "Snapshot of", dataset, ":\n", clinicals.head())

    ### Fix survival column to contain event codes
    clinicals = fix_gdc_survival_categories(clinicals, xena_hub, dataset)
    hint(verbose, "Survival after decoding events:\n", clinicals.head())

    ### Login to Drive and explore directory structure
    drive = connectDrive.io.login_to_drive(gdrive_secret)
    gd_projectFolderID, datasetDir = connectDrive.io.explore_drive_content(drive)
    hint(verbose, "Data ID of the curated survival table:\n", datasetDir[survfile])

    ### Acces the curated set of survival information
    clinicals = curatedSurvival(survfile, samplefile, drive, datasetDir)
    hint(verbose, "Data in the curated survival table:\n", clinicals.head())

    ### Add gene expression to the table
    clinicals = gex_tools.add_gene_expression_by_genes([gene], clinicals, xena_hub, gex_dataset)
    hint(
        verbose,
        "Gene expression in",
        gex_dataset,
        "retrived directly on genes:\n",
        clinicals.head(),
    )

    ### Split samples by median expression
    clinicals = gex_tools.split_by_gex_median(clinicals, gene)
    hint(verbose, "Gene expression categories based on the median:\n", clinicals.head())

    ### Create a mask that defines our groups
    mask = (clinicals["cat_"+gene] == 'low')

    ### Calculate logRank statistics to compare survival in the two groups
    stat = logRankSurvival(clinicals['time'], clinicals['event'], mask)
    hint(verbose, "The probability for the two groups to have the same survival:\n", stat.p_value)

    ### Plot survival
    plotting_tools.set_figure_rc()
    ax = plotKMpair(clinicals, mask)

    return ax

def plotKMpair(
    df: pd.DataFrame,
    mask: pd.Series,
    *,
    alternative_mask: Union[None, pd.Series] = None,
    timeline: Union[None, Sequence] = None,
    labels: Union[None, Tuple[str, str]] = ('low expression', 'high expression'),
    xlabel: str = 'Overall survival (months)',
    title: str = '',
    calculate_stat: bool = True,
    make_legend: bool = True,
    ax: Union[None, plt.Axes] = None,
    ) -> plt.Axes:
        
    """
    Plots two Kaplan-Meier curves to compare survival in two groups.

    Parameters
    ----------
    df
        A data frame with two compulsory columns: time and event.
    mask
        A Pandas mask for the main group.
    alternative_mask
        The second group is the negation of the first mask
        by default. This parameter sets a custom mask.
    timeline
        The range of survival that has to be shown.
    labels
        Legend label for the low and the high expression.
    xlabel
        Label for the x axis.
    title
        Title of the (sub)plot.
    calculate_stat
        If a logrank statistic should be calculated.
    legend
        If a legend should be added to the plot.
    ax
        The matplotlib axis object for the plot.

    Returns
    -------
    The matplotlib axis object with the Kaplan-Meier curves.
    """
    
    if ax is None:
        fig, ax = plt.subplots()
    if alternative_mask is None:
        alternative_mask = ~mask
    l1, l2 = labels

    T = df['time']
    E = df['event']
    kmf = KaplanMeierFitter()
    kmf.fit(T[mask], E[mask], timeline=timeline, label=l1 + '(n=' + str(len(E[mask])) + ')')
    ax = kmf.plot(ax=ax, legend=make_legend)
    kmf.fit(T[alternative_mask], E[alternative_mask], timeline=timeline, label=l2 + '(n=' + str(len(E[alternative_mask])) + ')')
    ax = kmf.plot(ax=ax, legend=make_legend)
    if calculate_stat:
        s = logRankSurvival(T, E, mask, alternative_mask=alternative_mask)
        title += '\n(p={:1.6f})'.format(s.p_value)
    ax.set_title(title, y=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylim(0, 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if make_legend:
        ax.legend(title='', loc='lower left', frameon=False)
    else:
        try:
            ax.get_legend().remove()
        except:
            pass
    return ax

def curatedSurvival(
    cfn: str,
    sfn: str,
    drive: pd.Series,
    drivedir: dict,
    *,
    survivaltype: Tuple[str, str] = ('OS', 'OS.time'),
    clinical_indexcol: str = 'bcr_patient_barcode',
    sample_indexcol: str = 'sample',
    sampletypecol: str = 'sample_type',
    sampletype: Union[None, List] = gdc_features.gdc_any_tumor,
    timefactor: int = 30,
    ) -> pd.DataFrame:
        
    """
    Read clinical information from the dataset curated by Liu et al. in a 2018 Cell
    paper (https://www.ncbi.nlm.nih.gov/pubmed/29625055). Additionally, filter for tumor
    samples.

    Parameters
    ----------
        cfn
        Name of the clinical information file.
    sfn
        Name of the sample type information file.
    drive
        A Google Drive object after authentication.
    drivedir
        The directory structure of the drive (filne names mapped to IDs).
    survivaltype
        If overall disease-free or other survival is investigated. A tuple of column
        names for the event and the time.
    clinical_indexcol
        The index column with patient IDs.
    sample_indexcol
        Index of the sample type table, with sample IDs.
    sampletypecol
        Column describing sample type.
    sampletype
        Sample types accepted for the analysis. If None, all samples are considered.
    timefactor
        Scaling factor for time. Converts days to months by default.

    Returns
    -------
    A curated table of patient survival
    """

    clinicals = connectDrive.io.read_gd_table(drive, cfn, drivedir)
    clinicals = clinicals.set_index(clinical_indexcol)
    clinicals['event'] = clinicals[survivaltype[0]]
    clinicals['time'] = clinicals[survivaltype[1]]
    clinicals = clinicals.loc[clinicals['time'] != 'NaN', :]
    clinicals = clinicals.loc[clinicals['time'].notna(), :]
    clinicals = clinicals.loc[clinicals['event'] != 'NaN', :]
    clinicals = clinicals.loc[clinicals['event'].notna(), :]
    clinicals['time'] = clinicals['time'].astype(float)
    clinicals['time'] = clinicals['time']/timefactor
    clinicals = clinicals.iloc[:,1:]
    saved_columns = clinicals.columns.values.tolist()
    clinicals['tmp_col'] = clinicals.values.tolist()
    clinicals = clinicals.loc[:, 'tmp_col'].to_dict()
    sampleTypes = connectDrive.io.read_gd_table(drive, sfn, drivedir)
    sampleTypes = sampleTypes[~pd.isnull(sampleTypes[sample_indexcol])]
    sampleTypes['patient'] = sampleTypes[sample_indexcol].apply(lambda x: '-'.join(x.split('-')[:3]))
    if sampletype is not None:
        sampleTypes = sampleTypes.loc[sampleTypes[sampletypecol].isin(sampletype),:]
    sampleTypes = sampleTypes.loc[sampleTypes['patient'].isin(clinicals.keys()),:]
    sampleTypes['tmp_col'] = sampleTypes['patient'].map(clinicals)
    sampleTypes[saved_columns] = pd.DataFrame(sampleTypes['tmp_col'].tolist(), index=sampleTypes.index)
    sampleTypes = sampleTypes.set_index(sample_indexcol)
    sampleTypes = sampleTypes.loc[:,saved_columns+[sampletypecol, 'patient']]
    return sampleTypes

def logRankSurvival(
    T: pd.Series,
    E: pd.Series,
    mask: pd.Series,
    *,
    alternative_mask: Union[None, pd.Series] = None,
    ) -> None:
        
    """
    Calculated a logRank p-value of survival in two conditions
    being different.

    Parameters
    ----------
    T
        A series of survival times.
    E
        A series of events, where 1 is the event (death).
    mask
        A Pandas mask for the main group.
    alternative_mask
        The second group is the negation of the first mask
        by default. This parameter sets a custom mask.

    Returns
    -------
    A logrank statistics result.
    """
    
    if alternative_mask is None:
        alternative_mask = ~mask
    return statistics.logrank_test(T[mask], T[alternative_mask], event_observed_A=E[mask], event_observed_B=E[alternative_mask])

def fix_gdc_survival_categories(
    clinicals: pd.DataFrame,
    xena_hub: str,
    ds: str,
    *,
    leveldict: dict = {'NaN': 0},
    renamedict: Union[None, dict] = {'Dead': 1},
    renameremain: Union[None, str, int, float] = 0,
    survcol: str = 'vital_status.demographic',
    timecol: str = 'days_to_death.demographic',
    alttimecol: str = 'days_to_last_follow_up.diagnoses',
    timefactor: int = 30,
    ) -> pd.DataFrame:
        
    """
    Recodes survival in 0-1 events and time scaled by given factor (months by default).
    Drops rows where either time or event information is missing.

    Parameters
    ----------
    clinicals
        A dataframe where the survival column is to be recoded.
    xena_hub
        Url of the data repository hub.
    ds
        Name of the dataset on the repository hub.
    leveldict
        A dictionary to override some level naming if needed. Must have a NaN key.
    renamedict
        A dictionary to set an alternative name for categories (e.g. 1 for Deceased).
    renameremain
        If the category name is not in the rename dictionary, it keeps its name
        by default. If renameremain is not None, then remaning items get this name.
    survcol
        Column with survival data. Transformed and copied into the event column.
    timecol
        Primary survival time column. Usually Nan for surviving patients.
    alttimecol
        Secondary survival time column that contains information also for surviving
        patients. Typically the 
    timefactor
        Scaling factor for time. Converts days to months by default.

    Returns
    -------
    Dataframe with time and event columns, dropping missing value rows.
    """

    clinicals = xena_tools.fix_phenotype_factorlevels(clinicals, xena_hub, ds, survcol, leveldict=leveldict, renamedict=renamedict, renameremain=renameremain)
    clinicals['event'] = clinicals[survcol]
    clinicals = clinicals.loc[clinicals['event'] != 'NaN', :]
    clinicals['time'] = clinicals.apply(lambda x: x[timecol] if x[alttimecol] != 'NaN' else x[timecol], axis=1)
    clinicals = clinicals.loc[clinicals['time'] != 'NaN', :]
    clinicals['time'] = clinicals['time'].astype(float)
    clinicals['time'] = clinicals['time']/timefactor
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