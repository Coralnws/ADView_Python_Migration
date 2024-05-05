import codecs
import os
from setuptools import setup, find_packages


with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.1'
DESCRIPTION = ('pip package to visually compare phylogenetic trees, utilizing Aggregated Dendrogram for phylogenetic tree visualization. ')
LONG_DESCRIPTION = 'pip package to visually compare phylogenetic trees, utilizing Aggregated Dendrogram for phylogenetic tree visualization. '


setup(
    name="ADViewpy",
    version=VERSION,
    author="Ng Weng Shan",
    author_email="",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[
        'dendropy',
        'ipycanvas',
        'ipywidgets',
        'scikit-learn',
        'numpy'
    ],
    keywords=['python', 'phylogenetic tree', 'aggregrated dendrogram','tree comparison'],
)
