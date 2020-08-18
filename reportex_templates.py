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
from typing import Union


def general_report(
    *,
    report_title: str = "Report",
    author_name: str = "Author",
    lab_name: str = "Lab",
    heat_table: Union[None, str] = None,
    left_tree: str = "",
    bottom_tree: str = "",
    plot_list: list = [],
    pagetitles: Union[None, str] = None,
    cohort_order: Union[None, list] = None,
    method_comments: str = "",
    bibliography: Union[None, str] = None,
    verbose: bool = True,
) -> str:

    """
    LaTeX template for reports that sum up the most significant results.

    Parameters
    ----------
    textable
        Latex table of p-values.
    img
        Heatmap in image format.
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
    verbose
        Flag to turn on messages displaying intermediate results.
    
    Returns
    -------
    A .tex file.
    
    """

    if heat_table is None:
        heat_table = """
        \\begin{tabular}{c | c | c}
        & \\rotatebox{90}{TP53} & \\rotatebox{90}{MKI67}\\\\ 
        TCGA-BRCA & \\color[rgb]{0.03137254901960784,0.18823529411764706,0.4196078431372549} \\cellcolor[rgb]{1.0,0.9607843137254902,0.9411764705882353} \\hyperlink{page.2}{0.49} & \\color[rgb]{0.03137254901960784,0.18823529411764706,0.4196078431372549} \\cellcolor[rgb]{1.0,0.9607843137254902,0.9411764705882353} \\hyperlink{page.2}{0.13}\\\\ 
        TCGA-LIHC & \\color[rgb]{0.03137254901960784,0.18823529411764706,0.4196078431372549} \\cellcolor[rgb]{1.0,0.9607843137254902,0.9411764705882353} \\hyperlink{page.3}{0.44} & \\color[rgb]{0.03137254901960784,0.18823529411764706,0.4196078431372549} \\cellcolor[rgb]{1.0,0.9607843137254902,0.9411764705882353} \\hyperlink{page.3}{3.38}\\\\ 
        \\end{tabular}
        """

    if bibliography is None:
        bibliography = ("", "")
    else:
        bibliography = (
            "\\addbibresource{" + bibliography + "}",
            "\n    \\flushleft\n    \\leftskip=60pt\n    \\rightskip=60pt\n    \\printbibliography[title={Citations for data and tools}]",
        )

    if pagetitles is None:
        pagetitles = dict()
    else:
        pagetitles = pd.read_csv(pagetitles, sep="\t", names=["cohort", "title"])
        pagetitles = pagetitles.set_index("cohort")
        pagetitles = pagetitles.to_dict()
        pagetitles = pagetitles["title"]

    page_list = []
    chi = 0
    if cohort_order is None:
        for plots in plot_list:
            if type(plots) is list:
                path = plots[0]
            else:
                path = plots
            e = os.path.basename(path)
            e = e[:-4]
            e = e.split("_")[0]
            if e in pagetitles:
                e = pagetitles[e]
            page_list.append((plots, e))
    else:
        d = dict()
        for plots in plot_list:
            if type(plots) is list:
                path = plots[0]
            else:
                path = plots
            e = os.path.basename(path)
            e = e[:-4]
            e = e.split("_")[0]
            if type(plots) is list:
                a = []
                for b in plots: # filter for nested fig lists here
                    if b[-4:] == ".txt":
                        with open(b, "r") as f:
                            for line in f:
                                print([line])
                                a.append(os.path.realpath(line))
                    a.append(os.path.realpath(b))
                d[e] = a
            else:
                d[e] = plots
        for cohort in cohort_order:
            if cohort in pagetitles:
                page_list.append((d[cohort], pagetitles[cohort]))
            else:
                page_list.append((d[cohort], cohort))

    bottom_tree = os.path.realpath(bottom_tree)
    left_tree = os.path.realpath(left_tree)

    tex = (
        """
    \\documentclass[10pt]{article}
    \\usepackage[a4paper, total={7.5in, 10in}]{geometry}
    \\usepackage[utf8]{inputenc}
    \\usepackage{pgfplots}
    \\usepackage{hyperref}
    \\usepackage{colortbl}
    \\usepackage{subcaption}
    \\usepackage{fancyhdr}
    \\usepackage{authoraftertitle}
    \\usepackage[style=nature]{biblatex}
    """
        + bibliography[0]
        + """

    \\defbibenvironment{bibliography}
    {\\list
        {\\printtext[labelnumberwidth]{%
        \\printfield{labelprefix}%
        \\printfield{labelnumber}}}
        {\\setlength{\\labelwidth}{\\labelnumberwidth}%
        \\setlength{\\leftmargin}{\\labelwidth}%
        \\setlength{\\labelsep}{\\biblabelsep}%
        \\addtolength{\\leftmargin}{\\labelsep}%
        \\addtolength{\\leftmargin}{\\leftskip}%
        \\addtolength{\\rightmargin}{\\rightskip}%
        \\setlength{\\itemsep}{\\bibitemsep}%
        \\setlength{\\parsep}{\\bibparsep}}%
        \\renewcommand*{\\makelabel}[1]{\\hss##1}}
    {\\endlist}
    {\\item}

    \\pagestyle{fancy}
    \\fancyhf{}
    \\rhead{"""
        + lab_name
        + """ \\hspace{60pt}}
    \\lhead{\hspace{60pt} \\MyTitle}
    \\rfoot{Generated by \\MyAuthor \\hspace{0pt} on \\MyDate \\hspace{60pt}}
    \\lfoot{\hspace{60pt} Page \\thepage}

    \\setlength{\\tabcolsep}{1pt}

    \\makeatletter
    \\def\@maketitle{%
    \\newpage
    \\null
    \\vskip -1em%
    \\begin{center}%
    \\let \\footnote \\thanks
        {\\LARGE \@title \\par}%
        \\vskip 2em%
        {\\large © """
        + lab_name
        + """, as of \\@date}%
    \\end{center}%
    \\par
    \\vskip 2em
    }
    \\makeatother

    \\title{"""
        + report_title
        + """}
    \\author{"""
        + author_name
        + """}
    \\date{\\today}

    \\begin{document}

    \\maketitle
    \\thispagestyle{fancy}

    \\vspace*{-90pt}
    \\begin{table}[!ht]
    \\begin{subtable}{0.1\\textwidth}
        \\begin{flushright}
            \\vspace*{50pt}
            \\hspace*{20pt}
            \\input{"""
        + left_tree
        + """}
        \\end{flushright}
    \\end{subtable}
    \\begin{subtable}{0.9\\textwidth}\\scriptsize
    \\renewcommand{\\arraystretch}{1.8}

    """
        + heat_table
        + """
    
    \\end{subtable}
    \\begin{subtable}{\\textwidth}
        \\vspace*{-110pt}
        \\hspace*{45pt}\\input{"""
        + bottom_tree
        + """}
    \\end{subtable}
    \\end{table}

    \\vspace*{-30pt}
    \\begin{center}
    \\large
    Color represents -log10(p value) of the log-rank test. 
    \\end{center}



    """
        + "\n".join(
            [
                "\\newpage\n\\begin{center}\n\\section{"
                + pagetitle
                + "}"
                + "\\end{center}\\vspace*{-30pt}\n"
                + "\n".join(["\\input{" + plot + "}\n" for plot in plots])
                + "\\par\n\\vspace*{30pt}\n\\centering\\hyperlink{page.1}{Back to overview table}\n"
                for plots, pagetitle in page_list
            ]
        )
        + """

    \\newpage

    \\flushleft
    \\leftskip=60pt
    \\rightskip=60pt

    """
        + method_comments
        + """
    
    """
        + bibliography[1]
        + """

    \end{document}
    """
    )
    return tex


recipe = general_report


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
