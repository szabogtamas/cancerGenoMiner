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
from typing import Union, Tuple, Callable
from cancerGenoMiner import (
    environment_definiton,
    par_examples,
    gdc_features,
    plotting_tools,
    survival_tools,
    xena_tools,
    gex_tools,
    reporting_tools,
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

    Parameters
    ----------
    nodes
        Objects that define processes as nodes linked by Nextflow. If None is supplied,
        the pipeline will consist of the nodes defined here.
    replacement_nodes
        Offers a way of replacing only some processes with custom objects by supplying a
        dictionary with [name_in_default_list]:[custom node] key: value pairs.
    container_paths
        Paths to containers where some processes need to run.
    conda
        Path to a yaml file to be used during environment creation or the conda dir.
    
    Returns
    -------
    List of initialized process objects
    """

    default_nodes = [
        fetchClinicalFile(
            inchannels=["survtab"], outchannels=["local_survtab"], conda=conda
        ),
        getSurvival(
            inchannels=["cohorts", "genes", "genedict", "nicer_survtab"],
            outchannels=["survivals", "gene_nd"],
            conda=conda,
        ),
        plotSurvival(
            inchannels=["survivals", "plotgenes", "symdict", "numplot", "mpl_backend",],
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
        makeHeatmap(
            inchannels=["stat", "heatmap_f"],
            outchannels=["heatmap", "heatimage", "cohort_order"],
            conda=conda,
        ),
        pptFromFigures(conda=conda),
        getRefs(
            inchannels=["bibliography", "drive_key"],
            outchannels=["bibtex"],
            conda=conda,
        ),
        compileReport(
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
        pdfFromLatex(),
    ]

    return introSpect.flowNodes.checkNodeReplacements(
        nodes, default_nodes, replacement_nodes
    )


def create_pipeline(
    *,
    location: str = os.getcwd(),
    nodes: Union[None, list] = None,
    node_initializer: Callable = enlist_process_nodes,
    replacement_nodes: Union[None, dict] = None,
    main_kws: Union[None, dict] = None,
    default_main_kws: dict = default_main_kws,
    comment_on_methods: Union[None, str] = None,
    comment_location: Union[None, str] = None,
    conda: Union[None, str] = None,
    general_configs: Union[None, dict] = None,
    container_paths: Union[None, dict] = None,
    verbose: bool = True,
) -> str:
    """
    Create a Nextflow pipeline that take a list of genes and a list of cohorts, compares
    survival in high and low expression groups, plots a heatmap of survival impact and 
    also Kaplan-Meier plots for each gene in every cohort.

    Parameters
    ----------
    location
        Path where the pipeline directory should be created.
    nodes
        Objects that define processes as nodes linked by Nextflow. If set to None, it
        checks for globally defined node list.
    node_initializer
        A function that takes 4 arguments (nodes, replacement nodes, container paths and
        conda) and returns a list of initialized node objects.
    replacement_nodes
        Offers a way of replacing only some processes with custom objects by supplying a
        dictionary with [name_in_default_list]:[custom node] key: value pairs.
    main_kws
        Initial pipeline parameters.
    default_main_kws
        Default values for the initial pipeline parameters.
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

    ### Add process nodes and compile the pipeline into a temporary folder
    nodes = node_initializer(nodes, replacement_nodes, container_paths, conda)
    introSpect.flowNodes.channelNodes(
        *nodes,
        main_kws=main_kws,
        location=location + "/pipeline",
        generalSettings=general_configs,
        containerPaths=container_paths,
        verbose=verbose,
    )
    hint(verbose, "Pipeline compiled to folder", location + "/pipeline")

    ### Small adjustments to get comment locations right
    if comment_location is None:
        if "comments" in main_kws:
            comment_location = main_kws["comments"]
        else:
            comment_location = os.getcwd() + "/pipeline/comments.tex"
    else:
        main_kws["comments"] = comment_location

    ### Save comments to a file
    if comment_on_methods is None:
        comment_on_methods = ""
    with open(comment_location, "w") as f:
        f.write(comment_on_methods)

    return os.path.dirname(location + "/pipeline/main.nf")


recipe = create_pipeline


class fetchClinicalFile(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import os",
                "import connectDrive",
                "from typing import Union",
                "from cancerGenoMiner import par_examples, survival_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder, drive_folder],
        }

    def channel_specifications(self):
        return {
            "survtab": ("val", "survtab", "survival_table", None, True),
            "local_survtab": ("stdout", "survival_table", None, None, False),
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
        }
        return None

    def process(
        self,
        *,
        fn: str = "survival_table.tsv",
        survival_table: Union[None, str] = None,
        drive_key: str = par_examples.gdrive_secret,
    ) -> str:

        """
        Download survival table from Google drive.

        Parameters
        ----------
        survival_table
            Manual curated survival data outside the scope of UCSC Xena.
        drive_key
            Location of the credential file for Google Drive access.
        
        Returns
        -------
        Manually curated clinical endpoints.
        
        """

        if survival_table in ["None", None]:
            return "None"
        else:
            drive = connectDrive.io.login_to_drive(drive_key)
            gd_projectFolderID, datasetDir = connectDrive.io.explore_drive_content(
                drive
            )
            clinicals = survival_tools.curatedSurvival(
                survival_table, par_examples.pancan_sampletypes, drive, datasetDir
            )
            clinicals.to_csv(fn, sep="\t")
        return os.path.realpath(fn)


class getSurvival(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import connectDrive",
                "from typing import Union, Tuple",
                "from cancerGenoMiner import par_examples, gdc_features, survival_tools, xena_tools, gex_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder, drive_folder],
        }

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "echo": "true"}

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from(params.cohorts)",
                "map{['"
                + '"'
                + "'+"
                + 'params.report_title.replaceAll("^\\"|\\"\\$", "") + '
                + "' (' +"
                + "it +')"
                + '"'
                + "'"
                + ", it+'.md', it, it+'_expression']}",
                "set{cohorts}",
            ],
            [
                "Channel",
                "value(params.genes)",
                "map{it.join('" + '","' + "')}",
                "map{'" + '"' + "'+ it +'" + '"' + "'}",
                "into{genes; plotgenes}",
            ],
            [
                "Channel",
                "from(params.genedict)",
                "map{it.join('\\t')}",
                "collectFile(name: 'genedict.tsv', newLine: true)",
                "into{genedict; symdict}",
            ],
            ["local_survtab", "map{it.trim()}", "set{nicer_survtab}",],
        ]

    def channel_specifications(self):
        return {
            "genedict": ("each", "*genedict", "genedict", None, False),
            "genes": ("val", "genes", "genes", None, False),
            "cohorts": (
                "tuple",
                ("note_title", "note_name", "cohort", "xplot"),
                (None, None, "cohort", None),
                None,
                False,
            ),
            "nicer_survtab": ("each", "survtab", "survival_table", None, False),
            "survivals": (
                "tuple",
                ("note_title", "note_name", "cohort", "xplot", '"${cohort}_data.tsv"'),
                (None, None, None, None, "outFile"),
                None,
                False,
            ),
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
        genes: list = par_examples.target_genes,
        genedict: Union[None, str] = None,
        survival_table: Union[None, str] = None,
        gex_basis: str = "gene",
    ) -> Tuple[gex_tools.pd.DataFrame, str]:

        """
        Check survival of patients in a cohort and link gene expression to this data.

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
        genes
            The genes to be queried.
        genedict
            A table mapping gene names to gene symbols.
        survival_table
            Manual curated survival data outside the scope of UCSC Xena.
        gex_basis
            If gene average or probe values should be checked Type `probes` for probes.
        
        Returns
        -------
        The Kaplan–Meier plot.
        
        """

        dataset = cohort + phenotype_prefix
        gex_dataset = cohort + gex_prefix
        ch = ""
        if len(cohort.split("-")) > 1:
            ch = cohort.split("-")[1]
        if genedict is not None:
            symbols = gex_tools.map_genenames(genes, genedict)
            with open(genedict, "r") as f:
                gd = f.read()
        else:
            symbols = genes[:]
            gd = "NaN\tNaN"
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
                symbols, clinicals, xena_hub, gex_dataset
            )
        else:
            clinicals = gex_tools.add_gene_expression_by_probes(
                symbols, clinicals, xena_hub, gex_dataset
            )
        return clinicals, gd


class plotSurvival(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import scipy",
                "import numpy as np",
                "from typing import Union, Tuple",
                "from cancerGenoMiner import par_examples, gdc_features, plotting_tools, survival_tools, xena_tools, gex_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder],
        }

    def directives(self):
        return {"publishDir": "'../notebooks', mode: 'copy'" + ', pattern: "*.md"'}

    def channel_specifications(self):
        return {
            "symdict": ("each", "*symdict", "genedict", None, False),
            "mpl_backend": ("env", "MPLBACK", None, None, True),
            "plotgenes": ("val", "plotgenes", "genes", None, False),
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
            "plotnames": ("val", "plotcohort", "outFile", None, False),
            "gexnames": ("val", "gexcohort", "gex", None, False),
            "images": ("file", '"*.png"', None, None, False),
            "plots": (
                "tuple",
                ('"${plotcohort}.pgf"', '"${gexcohort}.pgf"'),
                (None, None),
                None,
                False,
            ),
            "titles": ("file", '"${plotcohort}_title.txt"', "titles", None, False),
            "notebooks": ("file", "note_name", None, None, False),
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
            "lrt": (
                2,
                "--lrt",
                {
                    "dest": "lrt",
                    "help": "Temporary file to store results of log-rank test.",
                },
            ),
            "titles": (
                3,
                "--titles",
                {
                    "dest": "titles",
                    "help": "Temporary file to store results of page titles.",
                },
            ),
            "gex": (
                4,
                "--gex",
                {"dest": "gex", "help": "Distribution of gene expressions.",},
            ),
        }
        self.capturepars = ["note_title", "note_name"]
        return None

    def process(
        self,
        *,
        clinicals: str = "",
        cohort: str = par_examples.cohort,
        genes: list = par_examples.target_genes,
        genedict: Union[None, str] = None,
        plotrow: int = 5,
        plotcol: int = 4,
    ) -> Tuple[
        plotting_tools.plt.Axes, list, gex_tools.pd.DataFrame, plotting_tools.plt.Axes
    ]:

        """
        Convenience functions to plot and compare Kaplan–Meier curves of individuals
        split up in two groups.

        Parameters
        ----------
        clinicals
            Data frame with information on both survival and gene expression.
        cohort
            The TCGA cohort to check.
        genes
            The genes to be queried.
        genedict
            A table mapping gene names to gene symbols.
        plotrow
            Number of subplots in a row.
        plotcol
            Number of subplots in a column.
        
        Returns
        -------
        The Kaplan–Meier plot and gene expression distribution.
        
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

        ### Read the prefetched data table
        clinicals = gex_tools.pd.read_csv(clinicals, sep="\t")
        smallclinicals = clinicals.head()

        ### Plot survival for every gene
        plt = plotting_tools.plt
        plotting_tools.set_figure_rc()
        fig, axs = plt.subplots(plotrow, plotcol)
        axs = axs.flatten()

        for i in range(len(symbols)):
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
                survival_tools.plotKMpair(
                    cg, mask, title=gene + symbol, ax=ax, make_legend=False
                )
                stat = survival_tools.logRankSurvival(cg["time"], cg["event"], mask)
                stats.append(-1 * np.log10(stat.p_value))
            except:
                ax.text(0.2, 0.5, "Not enough data")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                stats.append(0.0)
        ax = plotting_tools.legend_only(ax=axs[-1])

        ### Create prognosis categories based on survival quartiles
        lowquart, highquart = clinicals.time.quantile([0.25, 0.75]).tolist()
        df = ["gex_" + symbol for symbol in symbols]
        df = clinicals.loc[:, ["time"] + df]
        df.columns = ["time"] + symbols
        df = df.melt(id_vars=["time"])
        df.columns = ["time", "gene", "gex"]
        df["prognosis"] = df["time"].apply(lambda x: "poor" if x <= lowquart else "")
        df["prognosis"] = df["prognosis"].astype(str) + df["time"].apply(
            lambda x: "good" if x >= highquart else ""
        )
        df["prognosis"] = df["prognosis"].apply(lambda x: "mild" if x == "" else x)

        ### Plot distribution of gene expression
        fig, gex = plt.subplots(figsize=(7.2, 3.6))
        gex = plotting_tools.sns.violinplot(
            x="gene",
            y="gex",
            hue="prognosis",
            hue_order=["poor", "mild", "good"],
            data=df,
            linewidth=0.2,
            ax=gex,
        )
        # gex = plotting_tools.sns.boxplot( x="gene", y="gex", hue="prognosis", hue_order=["poor", "mild", "good"], dodge=True, data=df, whis=np.inf, color="white", linewidth=0.4, ax=gex,)
        # gex = plotting_tools.sns.stripplot(x="gene", y="gex", hue="prognosis", hue_order=["poor", "mild", "good"], dodge=True, data=df, ax=gex, size=1, jitter=0.3,)
        pg1 = gex.scatter(0, 0, s=1, label="Poor prognosis")
        pg2 = gex.scatter(0, 0, s=1, label="Mild prognosis")
        pg3 = gex.scatter(0, 0, s=1, label="Good prognosis")
        gex.scatter(0, 0, color="white", s=1)
        gex.legend(handles=[pg1, pg2, pg3], loc="lower right")
        gex.set_xticklabels(
            [item.get_text() for item in gex.get_xticklabels()], rotation=30, ha="right"
        )
        gex.set_xlabel("")
        gex.set_ylabel("Gene expression (FPKM-UQ)", fontsize=9)
        gex.set_title(
            "Gene expression subset by survival quartiles\n(mild prognosis is 2nd and 3rd quartile)",
            fontsize=9,
        )

        bottom, top = gex.get_ylim()
        top = 0.9 * top
        for i, symbol in enumerate(symbols):
            t, p = scipy.stats.ttest_ind(
                df.loc[
                    (df["prognosis"] == "poor") & (df["gene"] == symbol), "gex"
                ].tolist(),
                df.loc[
                    (df["prognosis"] == "good") & (df["gene"] == symbol), "gex"
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

        return ax, [["cohort"] + genes, stats], titles, gex


class makeHeatmap(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import numpy as np",
                "from typing import Tuple",
                "from cancerGenoMiner import survival_tools, reporting_tools",
            ],
            "inhouse_packages": [cgm_folder, intro_folder],
        }

    def channel_pretreat(self):
        return [
            [
                "stats",
                "collectFile(name: 'stats.txt', newLine: true, keepHeader: true)",
                "set{stat}",
            ],
            [
                "Channel",
                "from('heatmap.tex', 'heatmap.png', 'tree1.pgf', 'tree2.pgf')",
                "toList()",
                "set{heatmap_f}",
            ],
        ]

    def channel_specifications(self):
        return {
            "stat": ("file", "stat", "fn", None, False),
            "heatmap_f": (
                "tuple",
                ("ltx", "img", "left", "bottom"),
                (None, None, None, None),
                None,
                False,
            ),
            "heatmap": (
                "tuple",
                ("'heatmap.tex'", "'tree1.pgf'", "'tree2.pgf'"),
                ("ltx", "left", "bottom"),
                None,
                False,
            ),
            "heatimage": ("file", "'heatmap.png'", "img", None, False),
            "cohort_order": ("file", "'cohort_order.txt'", "cohort_order", None, False),
        }

    def customize_features(self):
        self.modified_kws = {
            p[0]: (i + 1, "--" + p[0], {"dest": p[0], "help": p[1],})
            for i, p in enumerate(
                [
                    ("ltx", "LaTeX for a heatmap"),
                    ("img", "Heatmap image"),
                    ("left", "Left dendrogram"),
                    ("bottom", "Bottom dendrogram"),
                    ("cohort_order", "Order of cohort based on dendrogram"),
                ]
            )
        }
        return None

    def process(
        self, fn: str,
    ) -> Tuple[
        str,
        reporting_tools.sns.matrix.ClusterGrid,
        survival_tools.plt.Axes,
        survival_tools.plt.Axes,
    ]:
        """
        Make a table and heatmaps of statistics.

        Parameters
        ----------
        fn
            Table with p-values.
        
        Returns
        -------
        A latex format table, an image of the heatmap, plus dendrograms and cohort order.
        """

        ptable = reporting_tools.get_ptable(fn)
        row_linkage, col_linkage = reporting_tools.tree_linkages(ptable)
        ptable, left_tree, bottom_tree = reporting_tools.linkorder_table(
            ptable, row_linkage, col_linkage
        )
        latex_heatmap = reporting_tools.latex_heatmap(ptable)
        heatmap_image = reporting_tools.heatmap_image(ptable)
        return (
            latex_heatmap,
            heatmap_image,
            left_tree,
            bottom_tree,
            ",".join(ptable.index.values.tolist()),
        )


class getRefs(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import os",
                "from typing import Union",
                "from connectDrive import io",
            ],
            "inhouse_packages": [intro_folder, drive_folder],
        }

    def channel_specifications(self):
        return {
            "bibliography": ("val", "bibliography", "bib", None, True),
            "drive_key": ("val", "drive_key", "drive_key", None, True),
            "bibtex": ("file", "'biblio.bib'", "outFile", None, False),
        }

    def process(
        self, *, bib: Union[None, str] = None, drive_key: str = "~/.key",
    ) -> str:
        """
        Save references to a bibtex file either form a path or from Google drive.

        Parameters
        ----------
        bib
            Bibliography in bobtex format.
        drive_key
            Credentials to access Google Drive.
        
        Returns
        -------
        Bibtex with references.
        """

        if bib is None:
            return ""
        else:
            if os.path.isfile(bib):
                with open(bib, "r") as f:
                    bibtex = f.read()
            else:
                bibtex = io.read_file_by_id(bib, drive_key)

        return bibtex


class pptFromFigures(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def directives(self):
        return {
            "publishDir": "'../', mode: 'copy'",
        }

    def channel_pretreat(self):
        return [
            [
                "heatimage",
                "map{'"
                + '"'
                + "' + sdf.format(date) + '\\n\\n---\\n\\n![Heatmap of p-values](' + it +')\\n"
                + '"'
                + "'}",
                "set{slide_head}",
            ],
            [
                "images",
                "flatten()",
                "map{[it.getName().replaceAll('.png', ''), it]}",
                "map{it.join('](')}",
                "toList()",
                "map{it.join(')\\n\\n---\\n\\n![')}",
                "map{'" + '"' + "' + '\\n\\n---\\n\\n![' + it +')\\n" + '"' + "'}",
                "set{slide_body}",
            ],
            ["Channel", "fromPath(params.slidesfile)", "set{slide_show}",],
        ]

    def customize_features(self):
        self.manualDoc = "Collect output figures and save them in a ppt\n"
        self.inputs = [
            "val report_title from params.report_title",
            "val author_name from params.author_name",
            "val slide_head",
            "val slide_body",
            "file slide_show",
        ]
        self.outputs = ["file slide_show into slides"]
        self.command = (
            'echo "% " $report_title "\\n% " $author_name "\\n% " $slide_head $slide_body'
            " > figures.md" + "\n            "
            "pandoc figures.md -o $slide_show"
        )
        return None

    def compile_command(self):
        return self.command


class compileReport(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def dependencies(self):
        return {
            "imports": [
                "import numpy as np",
                "import os",
                "from typing import Union",
                "from cancerGenoMiner import reportex_templates",
            ],
            "inhouse_packages": [cgm_folder, intro_folder],
        }

    def channel_pretreat(self):
        return [
            [
                "titles",
                "collectFile(name: 'titles.txt', newLine: true)",
                "set{page_titles}",
            ],
            ["plots", "map{it.join(':')}", "toList()", "set{plotsr}",],
            ["cohort_order", "map{it.text}", "set{ordered_cohorts}",],
        ]

    def channel_specifications(self):
        return {
            "heatmap": (
                "tuple",
                ("'heatmap.tex'", "'tree1.pgf'", "'tree2.pgf'"),
                ("textable", "left", "bottom"),
                None,
                False,
            ),
            "plotsr": ("val", "plots", "plotlist", None, False),
            "page_titles": ("val", "pages", "pagetitles", None, False),
            "bibtex": ("file", "bibtex", "bibliography", None, False),
            "comments": ("val", "*comments", "method_comments", None, True),
            "ordered_cohorts": ("val", "cohort_order", "cohort_order", None, False),
            "reportex": ("file", "'report.tex'", "outFile", None, False),
        }

    def process(
        self,
        textable: str,
        left: str,
        bottom: str,
        plotlist: list,
        *,
        report_title: str = "Report",
        author_name: str = "Author",
        lab_name: str = "Lab",
        pagetitles: Union[None, str] = None,
        cohort_order: Union[None, list] = None,
        method_comments: Union[None, str] = None,
        bibliography: Union[None, str] = None,
    ) -> str:
        """
        Compile a latex file that sums up all the results.

        Parameters
        ----------
        textable
            Latex table of p-values.
        left
            Left denrogram.
        bottom
            Bottom dendrogram.
        plotlist
            List of plots to include as separate pages.
        report_title
            Title of the report, describing the aim of the analysis.
        author_name
            Name of the author.
        lab_name
            Name of the lab.
        page_title
            A table providing alternative (to cohort name) titles for the plot pages.
        cohort_order
            Order of cohorts after hierarchical clustering (as shown on the heatmap).
        method_comments
            Description of methods.
        bibliography
            Path to the bibtex file with references.
        
        Returns
        -------
        A latex file for the report.
        """

        report_title = report_title.strip("'")
        report_title = report_title.strip('"')
        author_name = author_name.strip("'")
        author_name = author_name.strip('"')
        lab_name = lab_name.strip("'")
        lab_name = lab_name.strip('"')

        plotlist = [x.replace("[", "") for x in plotlist]
        plotlist = [x.replace("]", "") for x in plotlist]
        plotlist = [x.split(":") for x in plotlist]

        with open(textable, "r") as f:
            textable = f.read()
        with open(method_comments, "r") as f:
            method_comments = f.read()
        bibliography = os.path.realpath(bibliography)

        tex = reportex_templates.general_report(
            report_title=report_title,
            author_name=author_name,
            lab_name=lab_name,
            heat_table=textable,
            left_tree=left,
            bottom_tree=bottom,
            plot_list=plotlist,
            pagetitles=pagetitles,
            cohort_order=cohort_order,
            method_comments=method_comments,
            bibliography=bibliography,
        )
        return tex


class pdfFromLatex(nextflowProcess):
    """
    Nextflow process to execute the function below.
    """

    def directives(self):
        return {
            "publishDir": "'../', mode: 'copy'",
        }

    def customize_features(self):
        self.manualDoc = "Convert the LaTeX format report into pdf.\n"
        self.inputs = ["file reportex"]
        self.outputs = ['file("*.pdf")']
        self.command = "pdflatex -interaction nonstopmode -halt-on-error -file-line-error $reportex\n            "
        self.command += "biber ${reportex.baseName}\n            "
        self.command += "pdflatex -interaction nonstopmode -halt-on-error -file-line-error $reportex\n"

        self.container = environment_definiton.latex_container

        return None

    def compile_command(self):
        return self.command


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
