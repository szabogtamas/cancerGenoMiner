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
from typing import Union, Tuple, Sequence, Callable
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
    "mutlabel": "p53",
    "report_title": '"Title of report"',
    "author_name": '"Author"',
    "lab_name": '"Laboratory or PI"',
    "bibliography": "biblio.bib",
    "slidesfile": "figures.pptx",
    "comments": os.getcwd() + "/pipeline/comments.tex",
    "mpl_backend": "pgf",
    "fn": "survival_table.tsv",
}


def enlist_process_nodes(
    nodes: Union[None, list],
    replacement_nodes: Union[None, dict],
    container_paths: Union[None, dict],
    conda: Union[None, str],
) -> list:
    """
    Helper function returning a list of initialized process objects. Node list, container
    and conda locations gets passed to this function by the pipeline creator function.
    """

    default_nodes = [
        ana_surgex_single.fetchClinicalFile(
            inchannels=["survtab"], outchannels=["local_survtab"], conda=conda
        ),
        ana_surgex_single.getSurvival(
            inchannels=["cohorts", "genes", "genedict", "nicer_survtab"],
            outchannels=["survivals", "gene_nd"],
            conda=conda,
        ),
        plotSurvival(
            inchannels=[
                "survivals",
                "conditiontab",
                "plotgenes",
                "symdict",
                "numplot",
                "mpl_backend",
                "colors",
                "labels",
            ],
            outchannels=[
                "plotnames",
                "quadnames",
                "hcorrnames",
                "hdistnames",
                "mdistnames",
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
        # ana_surgex_single.pdfFromLatex(),
    ]

    return introSpect.flowNodes.checkNodeReplacements(
        nodes, default_nodes, replacement_nodes
    )


def create_pipeline(**kwargs):
    """
    This just wraps main parameters and nodes defined here into the pipeline creator
    function of `ana_surgex_single`.
    """

    kws = {
        "default_main_kws": default_main_kws,
        "node_initializer": enlist_process_nodes,
    }
    kws.update(kwargs)
    return ana_surgex_single.create_pipeline(**kws)


create_pipeline.__doc__ = ana_surgex_single.create_pipeline.__doc__
recipe = create_pipeline


class plotSurvival(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import scipy",
                "import numpy as np",
                "from typing import Union, Tuple, Sequence",
                "from cancerGenoMiner import par_examples, gdc_features, plotting_tools, survival_tools, xena_tools, gex_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder],
        }

    def directives(self):
        return {
            "echo": "true",
            "publishDir": "'../notebooks', mode: 'copy'" + ', pattern: "*.md"',
        }
        
    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from(params.colors)",
                "toList()",
                "map{it.join(',')}",
                "map{it.replaceAll(" + '"#", "\\\\#"' +")}",
                "set{colors}",
            ],
            [
                "Channel",
                "from(params.labels)",
                "toList()",
                "map{it.join(',')}",
                "set{labels}",
            ],
        ]

    def channel_specifications(self):
        return {
            "symdict": ("each", "*symdict", "genedict", None, False),
            "conditiontab": ("each", "conditiontab", "conditiontab", None, True),
            "mpl_backend": ("env", "MPLBACK", None, None, True),
            "plotgenes": ("val", "plotgenes", "genes", None, False),
            "colors": ("val", "colors", "colors", None, False),
            "labels": ("val", "labels", "labels", None, False),
            "survivals": (
                "tuple",
                (
                    "note_title",
                    "note_name",
                    "plotcohort",
                    "gexcohort",
                    '"${plotcohort}_data.tsv"',
                ),
                (None, None, "cohort", None, "clinicals"),
                None,
                False,
            ),
            "stats": ("file", '"${plotcohort}_stats.csv"', "lrt", None, False),
            "plotnames": ("val", '"${plotcohort}_km"', "kmp", None, False),
            "quadnames": ("val", '"${plotcohort}_quad"', "kmq", None, False),
            "hcorrnames": ("val", '"${plotcohort}_hcorr"', "hazardcorr", None, False),
            "hdistnames": ("val", '"${plotcohort}_hdist"', "hazarddist", None, False),
            "mdistnames": ("val", '"${plotcohort}_mdist"', "mutdist", None, False),
            "images": ("file", '"*.png"', None, None, False),
            "plots": (
                "tuple",
                (
                    '"${plotcohort}_km.pgf"',
                    '"${plotcohort}_quad.pgf"',
                    '"${plotcohort}_hcorr.pgf"',
                    '"${plotcohort}_hdist.pgf"',
                    '"${plotcohort}_mdist.pgf"',
                ),
                (None, None, None, None, None),
                None,
                False,
            ),
            "titles": ("file", '"${plotcohort}_title.txt"', "titles", None, False),
            "notebooks": ("file", "note_name", None, None, False),
        }

    def customize_features(self):
        self.modified_kws = {
            p[0]: (i + 1, "--" + p[0], {"dest": p[0], "help": p[1],})
            for i, p in enumerate(
                [
                    (
                        "kmp",
                        "Kaplan-Meier plots of 4 groups, by 2 dichotomous conditions",
                    ),
                    ("kmq", "Kaplan-Meier plot series with every pairs of conditions"),
                    (
                        "hazardcorr",
                        "Correlation of survived hazard with gene expression",
                    ),
                    ("hazarddist", "Distribution of gene expression in risk groups"),
                    ("mutdist", "Distribution of gene expression in WT and mutant"),
                    ("lrt", "Temporary file to store results of log-rank test"),
                    ("titles", "Temporary file to store results of page titles"),
                ]
            )
        }
        self.capturepars = ["note_title", "note_name"]
        return None

    def process(
        self,
        *,
        clinicals: str = "",
        conditiontab: str = "",
        cohort: str = par_examples.cohort,
        genes: list = par_examples.target_genes,
        genedict: Union[None, str] = None,
        mutlabel: Union[None, str] = None,
        labels: list = par_examples.quadKMlabels,
        colors: list = par_examples.quadKMcolors,
        plotrow: int = 5,
        plotcol: int = 4,
    ) -> Tuple[
        plotting_tools.plt.Axes, list, gex_tools.pd.DataFrame, plotting_tools.plt.Axes
    ]:

        """
        A process that plots survival split up by two conditions: typically low vs. high
        expression of a gene and the presence of mutations in another one.

        Parameters
        ----------
        clinicals
            Data frame with information on both survival and gene expression.
        conditiontab
            A list of samples with the seconary condition.
        cohort
            The TCGA cohort to check.
        genes
            The genes to be queried.
        genedict
            A table mapping gene names to gene symbols.
        labels
            Legend labels for the two factors and four conditions.
        colors
        Line colors to be used.
        plotrow
            Number of subplots in a row.
        plotcol
            Number of subplots in a column.
        
        Returns
        -------
        Kaplan–Meier plot and distributions of gene expression.
        """

        ### Retrieve gene symbols, if available
        if genedict is not None:
            symbols = gex_tools.map_genenames(genes, genedict)
        else:
            symbols = genes[:]

        ### Use diagnosis instead of cohort codes
        stats = [cohort]
        if cohort in gdc_features.tcga_cohort_diagnosis:
            titles = (
                cohort
                + "\tImpact of gene expression on survival in "
                + gdc_features.tcga_cohort_diagnosis[cohort]
                + " ("
                + cohort
                + ")"
            )
        else:
            titles = cohort + "\t" + cohort

        ### Read the list of sample affected by the condition (mutation)
        with open(conditiontab, "r",) as f:
            mutants = f.read().split("\n")
        if mutlabel is not None:
            labels[2] = mutlabel + " " + labels[4]
            labels[3] = mutlabel + " " + labels[5]
        else:
            mutlabel = None

        ### Read the prefetched data table, add patient info
        clinicals = gex_tools.pd.read_csv(clinicals, sep="\t")
        cmask = clinicals["sample"].isin(mutants)
        clinicals["mutation"] = clinicals["sample"].apply(
            lambda x: mutlabel if x in mutants else "WT"
        )
        clinicals = survival_tools.calcSurvHazardCat(clinicals, hazardcol="hazard")
        smallclinicals = clinicals.head()

        ### Reorder table to group gene expression by mutation status and hazard
        df = ["gex_" + symbol for symbol in symbols]
        commoncols = ["sample", "mutation", "hazard"]
        df = clinicals.loc[:, commoncols + df]
        df.columns = commoncols + symbols
        df = df.melt(id_vars=commoncols)
        df.columns = commoncols + ["gene", "gex"]

        ### Plot survival for every gene
        plt = plotting_tools.plt
        plotting_tools.set_figure_rc()
        fig, axs = plt.subplots(plotrow, plotcol)
        axs = axs.flatten()

        gN = len(symbols)
        for i in range(gN):
            gene = genes[i]
            gene = gene.replace('"', "")
            symbol = symbols[i]
            ax = axs[i]
            cg = clinicals.loc[clinicals["gex_" + symbol] != "NaN", :]
            cg = gex_tools.split_by_gex_median(cg, symbol)
            mask = cg["cat_" + symbol] == "low"
            if symbol == gene:
                symbol = ""
            else:
                symbol = " (" + symbol + ")"
            try:
                survival_tools.plotKMquad(
                    cg,
                    mask,
                    cmask,
                    title=gene + symbol,
                    ax=ax,
                    make_legend=False,
                    colors=colors,
                )
                statc = []
                statc.append(
                    survival_tools.logRankSurvival(
                        cg["time"], cg["event"], (mask & ~cmask)
                    ).p_value
                )
                statc.append(
                    survival_tools.logRankSurvival(
                        cg["time"], cg["event"], (~mask & ~cmask)
                    ).p_value
                )
                statc.append(
                    survival_tools.logRankSurvival(
                        cg["time"], cg["event"], (mask & cmask)
                    ).p_value
                )
                statc.append(
                    survival_tools.logRankSurvival(
                        cg["time"], cg["event"], (~mask & cmask)
                    ).p_value
                )
                stat = np.min(statc)
                stats.append(-1 * np.log10(stat))
            except:
                ax.text(0.2, 0.5, "Not enough data")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                stats.append(0.0)
        ax = plotting_tools.legend_only(ax=axs[-1], labels=labels[:4], colors=colors)

        ### Plot survival by factors
        fgs = []
        if gN % (plotrow - 1) == 0:
            fN = gN / (plotrow - 1)
        else:
            fN = gN / (plotrow - 1) + 1
        for fi in range(int(fN)):
            fig, axs = plt.subplots(plotrow, 5)
            axs = axs.flatten()
            naxs = []
            has_legendrow = False
            for i in range(plotrow - 1):
                sax = axs[5 * i : 5 * i + 5]
                pN = fi * (plotrow - 1) + i
                if pN < len(genes):
                    gene = genes[pN]
                    gene = gene.replace('"', "")
                    symbol = symbols[pN]
                    cg = clinicals.loc[clinicals["gex_" + symbol] != "NaN", :]
                    cg = gex_tools.split_by_gex_median(cg, symbol)
                    mask = cg["cat_" + symbol] == "low"
                    if symbol == gene:
                        symbol = ""
                    else:
                        symbol = " (" + symbol + ")"
                    sax = survival_tools.plotKMquads(
                        cg, mask, cmask, title=gene + symbol, make_legend=False, axs=sax
                    )
                    naxs.extend(sax)
                else:
                    if has_legendrow:
                        for nonax in sax:
                            nonax.axis("off")
                    else:
                        sax = plotting_tools.make_kmquad_legendrow(sax, labels, colors)
                        has_legendrow = True
                naxs.extend(sax)
            i += 1
            sax = axs[5 * i : 5 * i + 5]
            if has_legendrow:
                for nonax in sax:
                    nonax.axis("off")
            else:
                sax = plotting_tools.make_kmquad_legendrow(sax, labels, colors)
            fgs.append(naxs)

        ### Plot gene expression as a function of survived hazard
        plt = plotting_tools.plt
        plotting_tools.set_figure_rc()
        fig, haxs = plt.subplots(plotrow, plotcol)
        haxs = haxs.flatten()

        gN = len(symbols)
        for i in range(gN):
            gene = genes[i]
            gene = gene.replace('"', "")
            symbol = symbols[i]
            hax = haxs[i]
            cg = clinicals.loc[
                clinicals["gex_" + symbol] != "NaN",
                ["time", "event", "mutation", "gex_" + symbol, "hazard"],
            ]
            cg.columns = ["time", "event", "mutation", "gex", "hazard"]
            if symbol == gene:
                symbol = ""
            else:
                symbol = " (" + symbol + ")"
            hax = survival_tools.plotSurvHazardCat(
                cg,
                hazardcol="hazard",
                featcol="gex",
                catcol="mutation",
                colordict={"WT": colors[4], mutlabel: colors[5]},
                title=gene + symbol,
                make_legend=False,
                ax=hax,
            )
        hax = plotting_tools.legend_only(
            ax=haxs[-1], labels=["WT", mutlabel], colors=colors[4:6]
        )

        ### Label patients according to hazard: low survived hazard is high risk
        df["risk"] = df["hazard"].apply(lambda x: "high" if x <= 0.25 else "")
        df["risk"] = df["risk"].astype(str) + df["hazard"].apply(
            lambda x: "low" if x >= 0.75 else ""
        )
        df["risk"] = df["risk"].apply(lambda x: "mild" if x == "" else x)

        ### Plot distribution of gene expression in hazard groups
        fig, pax = plt.subplots(figsize=(7.2, 3.6))
        pax = plotting_tools.sns.violinplot(
            x="gene",
            y="gex",
            hue="risk",
            hue_order=["high", "mild", "low"],
            data=df,
            linewidth=0.2,
            ax=pax,
        )
        pg1 = pax.scatter(0, 0, s=1, label="High risk")
        pg2 = pax.scatter(0, 0, s=1, label="Mild prognosis")
        pg3 = pax.scatter(0, 0, s=1, label="Low risk")
        pax.scatter(0, 0, color="white", s=1)
        pax.legend(handles=[pg1, pg2, pg3], loc="lower right")
        pax.set_xticklabels(
            [item.get_text() for item in pax.get_xticklabels()], rotation=30, ha="right"
        )
        pax.set_xlabel("")
        pax.set_ylabel("Gene expression (FPKM-UQ)", fontsize=9)
        pax.set_title(
            "Gene expression subset by risk\n(survived hazard retrospectively)",
            fontsize=9,
        )

        ### Calculate statistics and add stars if significant
        bottom, top = pax.get_ylim()
        top = 0.9 * top
        for i, symbol in enumerate(symbols):
            t, p = scipy.stats.ttest_ind(
                df.loc[(df["risk"] == "low") & (df["gene"] == symbol), "gex"].tolist(),
                df.loc[(df["risk"] == "high") & (df["gene"] == symbol), "gex"].tolist(),
                equal_var=False,
            )
            s = "p={:1.5f}".format(p)
            if p < 0.05:
                s = "*"
                if p < 0.01:
                    s += "*"
                    if p < 0.001:
                        s += "*"
            else:
                s = ""
            pax.text(i, top, s)

        ### Plot distribution of gene expression
        fig, gex = plt.subplots(figsize=(7.2, 3.6))
        gex = plotting_tools.sns.violinplot(
            x="gene",
            y="gex",
            hue="mutation",
            hue_order=["WT", mutlabel],
            data=df,
            linewidth=0.2,
            ax=gex,
        )
        pg1 = gex.scatter(0, 0, s=1, label="WT")
        pg2 = gex.scatter(0, 0, s=1, label=mutlabel)
        gex.scatter(0, 0, color="white", s=1)
        gex.legend(handles=[pg1, pg2], loc="lower right")
        gex.set_xticklabels(
            [item.get_text() for item in gex.get_xticklabels()], rotation=30, ha="right"
        )
        gex.set_xlabel("")
        gex.set_ylabel("Gene expression (FPKM-UQ)", fontsize=9)
        gex.set_title(
            "Gene expression subset by mutation status\n(" + mutlabel + ")", fontsize=9,
        )

        bottom, top = gex.get_ylim()
        top = 0.9 * top
        for i, symbol in enumerate(symbols):
            t, p = scipy.stats.ttest_ind(
                df.loc[
                    (df["mutation"] == "WT") & (df["gene"] == symbol), "gex"
                ].tolist(),
                df.loc[
                    (df["mutation"] == mutlabel) & (df["gene"] == symbol), "gex"
                ].tolist(),
                equal_var=False,
            )
            s = "p={:1.5f}".format(p)
            if p < 0.05:
                s = "*"
                if p < 0.01:
                    s += "*"
                    if p < 0.001:
                        s += "*"
            else:
                s = ""
            gex.text(i, top, s)

        return ax, fgs[0], hax, pax, gex, [["cohort"] + genes, stats], titles


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
