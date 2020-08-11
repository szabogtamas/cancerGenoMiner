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

default_main_kws = {
    "cohorts": [par_examples.cohort] + par_examples.lung_cohorts,
    "genes": par_examples.prognosis_genes,
    "genedict": {"NaN": "NaN"},  # par_examples.prognosis_genecodes,
    "survtab": "None",
    "conditiontab": "None",
    "report_title": '"Title of report"',
    "author_name": '"Author"',
    "lab_name": '"Laboratory or PI"',
    "bibliography": "biblio.bib",
    "slidesfile": "figures.pptx",
    "comments": os.getcwd() + "/pipeline/comments.tex",
    "mpl_backend": "pgf",
    "fn": "survival_table.tsv",
}


def enlist_process_nodes(nodes, conda):
    nodes = [introSpect.flowNodes.helloWorld(inchannels=["cheers"])]
    return nodes

#TODO: It is not the best strategy to override the parent function. Need to come up with something better. Maybe exec?
#ana_surgex_single.enlist_process_nodes = enlist_process_nodes
ana_surgex_single.default_main_kws = default_main_kws
create_pipeline = ana_surgex_single.create_pipeline

recipe = create_pipeline


class rankSurvivalImpacts(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import numpy as np",
                "from typing import Union, Tuple",
                "from cancerGenoMiner import par_examples, gdc_features, survival_tools, xena_tools, gex_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder, drive_folder],
        }

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from(params.genedict)",
                "map{it.join('\\t')}",
                "collectFile(name: 'genedict.tsv', newLine: true)",
                "set{genedict}",
            ],
            ["local_survtab", "map{it.trim()}", "set{nicer_survtab}",],
        ]

    def channel_specifications(self):
        return {
            "genes": ("each", "gene", "gene", None, True),
            "genedict": ("each", "*genedict", "genedict", None, False),
            "nicer_survtab": ("each", "survtab", "survival_table", None, False),
            "cohorts": ("val", "cohort", "cohort", None, True,),
            "interactions": ("file", '"${cohort}_data.tsv"', "outFile", None, False,),
            "gene_nd": ("file", "'genedict.tsv'", "gd", None, False),
        }

    def customize_features(self):
        self.modified_kws = {
            "outFile": (
                1,
                "-o",
                "--outFile",
                {
                    "dest": "outFile",
                    "help": "Location where results should be saved. If not specified, STDOUT will be used.",
                },
            ),
            "gd": (
                2,
                "--gd",
                {"dest": "gd", "help": "Temporary file to store gene symbol mapping.",},
            ),
        }
        return None

    def process(
        self,
        *,
        xena_hub: str = par_examples.xena_hub,
        gex_prefix: str = par_examples.gextag,
        phenotype_prefix: str = par_examples.phenotypetag,
        cohort: str = par_examples.cohort,
        gene: str = par_examples.target_gene,
        genedict: Union[None, str] = None,
        geneslice: int = 500,
        survival_table: Union[None, str] = None,
        probemap: str = par_examples.probemap,
        gex_basis: str = "gene",
    ) -> Tuple[gex_tools.pd.DataFrame, str]:

        """
        Check the impact of other genes on survival benefit or disadvantage of a gene.

        Parameters
        ----------
        xena_hub
            Url of the data repository hub.
        gex_prefix
            Constant part of the gene expression dataset name.
        phenotype_prefix
            Constant part of the clinical phenotype dataset name.
        cohort
            The TCGA cohort to check.
        gene
            The gene to be queried.
        genedict
            A table mapping gene names to gene symbols.
        geneslice
            The size of gene batches to be querried.
        survival_table
            Manual curated survival data outside the scope of UCSC Xena.
        probemap
            A probemap file for testing direct table download.
        gex_basis
            If gene average or probe values should be checked Type `probes` for probes.
        
        Returns
        -------
        Interaction weights on a given gene for every gene in the genome.
        
        """

        # TODO: This module should introduce the second factor (p53 mutations, or CA20 score). The question is if this factor should be just a list of IDs or a dataframe, or a Series of actual values
        # The list of sample IDs should be read from a file

        dataset = cohort + phenotype_prefix
        gex_dataset = cohort + gex_prefix
        gene = gene.replace('"', "")
        ch = ""
        if len(cohort.split("-")) > 1:
            ch = cohort.split("-")[1]
        if genedict is not None:
            symbol = gex_tools.map_genenames([gene], genedict)[0]
            with open(genedict, "r") as f:
                gd = f.read()
        else:
            symbol = gene[:]
            gd = "NaN\tNaN"

        ### Retrieve clinical information
        if survival_table in ["None", None]:
            clinicals = xena_tools.download_gdc_clinicals(xena_hub, dataset)
            clinicals = xena_tools.fix_phenotype_factorlevels(
                clinicals,
                xena_hub,
                dataset,
                "sample_type.samples",
                leveldict={"NaN": "NaN"},
            )
            clinicals = clinicals.loc[
                clinicals["sample_type.samples"].isin(gdc_features.gdc_any_tumor), :
            ]
            clinicals = survival_tools.fix_gdc_survival_categories(
                clinicals, xena_hub, dataset
            )
        else:
            clinicals = survival_tools.pd.read_csv(survival_table, sep="\t")
            clinicals = clinicals.loc[clinicals["type"].isin([cohort, ch]), :]
            clinicals = clinicals.set_index("sample")
        if gex_basis == "gene":
            clinicals = gex_tools.add_gene_expression_by_genes(
                [symbol], clinicals, xena_hub, gex_dataset
            )
        else:
            clinicals = gex_tools.add_gene_expression_by_probes(
                [symbol], clinicals, xena_hub, gex_dataset
            )

        ### Add data on the expression of the focus gene
        cg = clinicals.loc[clinicals["gex_" + symbol] != "NaN", :]
        cg = gex_tools.split_by_gex_median(cg, symbol)
        mask = cg["cat_" + symbol] == "low"
        try:
            stat = survival_tools.logRankSurvival(cg["time"], cg["event"], mask)
            basestat = -1 * np.log10(stat.p_value)
        except:
            basestat = 0.0
        records = []

        ### Retrieve a list of all genes
        allgenes = np.array(
            xena_tools.xena.dataset_field_examples(xena_hub, gex_dataset, None)
        )
        N_genes = allgenes.shape[0]
        gsn = int(allgenes.shape[0] / geneslice)
        rest = allgenes[geneslice * gsn :]
        allgenes = allgenes[: geneslice * gsn].reshape(-1, geneslice).tolist()
        allgenes.append(rest.tolist())
        # allgenes = allgenes[:2]  ### For testing only!!!

        ### Create a mapping for ENS gene codes
        probes = xena_tools.read_xena_table(probemap, hubPrefix=xena_hub)
        probedict = gex_tools.parse_gene_mapping(probes, probecol="gene", genecol="id")

        ### Loop through genes in chunks
        for chunk in allgenes:
            etdf = gex_tools.add_gene_expression_by_probes(
                chunk, cg, xena_hub, gex_dataset
            )
            for i, interactor in enumerate(chunk):
                tdf = gex_tools.split_by_gex_median(etdf, interactor)
                bmask = tdf["cat_" + symbol] == "low"
                imask = tdf["cat_" + interactor] == "low"
                try:
                    stat = survival_tools.logRankSurvival(
                        tdf["time"], tdf["event"], imask
                    )
                    ibasestat = -1 * np.log10(stat.p_value)
                except:
                    ibasestat = 0.0
                try:
                    focus_enables = (
                        -1
                        * np.log10(
                            survival_tools.logRankSurvival(
                                tdf["time"],
                                tdf["event"],
                                (imask & ~bmask),
                                alternative_mask=(~imask & ~bmask),
                            ).p_value
                        )
                        - ibasestat
                    )
                except:
                    focus_enables = 0.0
                try:
                    focus_inhibits = (
                        -1
                        * np.log10(
                            survival_tools.logRankSurvival(
                                tdf["time"],
                                tdf["event"],
                                (imask & bmask),
                                alternative_mask=(~imask & bmask),
                            ).p_value
                        )
                        - ibasestat
                    )
                except:
                    focus_inhibits = 0.0
                try:
                    interactor_enables = (
                        -1
                        * np.log10(
                            survival_tools.logRankSurvival(
                                tdf["time"],
                                tdf["event"],
                                (~imask & bmask),
                                alternative_mask=(~imask & ~bmask),
                            ).p_value
                        )
                        - basestat
                    )
                except:
                    interactor_enables = 0.0
                try:
                    interactor_inhibits = (
                        -1
                        * np.log10(
                            survival_tools.logRankSurvival(
                                tdf["time"],
                                tdf["event"],
                                (imask & bmask),
                                alternative_mask=(imask & ~bmask),
                            ).p_value
                        )
                        - basestat
                    )
                except:
                    interactor_inhibits = 0.0
                records.append(
                    (
                        cohort,
                        symbol,
                        interactor,
                        focus_enables,
                        focus_inhibits,
                        interactor_enables,
                        interactor_inhibits,
                    )
                )
        records = gex_tools.pd.DataFrame(
            records,
            columns=[
                "cohort",
                "symbol",
                "interactor",
                "focus_enables",
                "focus_inhibits",
                "interactor_enables",
                "interactor_inhibits",
            ],
        )
        records["interactor"] = records["interactor"].map(probedict)

        return records, gd


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
