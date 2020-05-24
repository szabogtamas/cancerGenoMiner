import os

__doc__ = """
A python package facilitating large-scale access to cancer genomic datasets (exploits
xenaPython) and find genomic alterations associated with clinical phenotype, focusing on
survival (lifeLines).
"""

homefolder = os.path.abspath(os.path.dirname(__file__))

for module in os.listdir(homefolder):
    if module == "__init__.py" or module[-3:] != ".py":
        continue
    __import__(module[:-3], globals(), locals(), level=1)
del module
