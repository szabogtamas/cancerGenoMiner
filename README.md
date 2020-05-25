# cancerGenoMiner
Functions facilitating download and analysis of cancer genomics data from databases like TCGA or CCLE.
Project is under constant development and not yet production-ready!

## Interaction with data repositories
Database connectivity relies heavily on UCSC Xena and the Python API [xenaPython](https://github.com/ucscXena/xenaPython). Manually curated data tables, complied from distinct papers (like progression-free survival or CA20 score) can also be imported from Google Drive after setting up access with [PyDrive](https://pythonhosted.org/PyDrive/).

## Data analysis
Besides classicel dataanalysis libriaries, like numpy, scipy, pandas, seaborn and matplotlib, a fair bit of functions rely on [lifelines](https://github.com/CamDavidsonPilon/lifelines) for survival analysis.
