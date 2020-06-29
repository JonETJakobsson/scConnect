.. image:: https://readthedocs.org/projects/scconnect/badge/?version=latest
    :target: https://scconnect.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

=========
scConnect
=========

The end goal of scConnect is to provide information about putative interaction between cell populations 
in single cell RNA-seq data. 
To do so gene expression levels are used to asses precence or absence of ligands and receptors. 
This information is then integrated with ligand receptor interaction data from `Guide to pharmacology`__ 
to detect putative connections.
The resulting NetworkX graph can be analysed and browsed using a dash web app.

__ https://www.guidetopharmacology.org/

?????????????
Documentation
?????????????
For helpful tutorials and information about the API please consult the `documentation`__.

__ https://scconnect.readthedocs.io/en/latest/

????????????
Installation
????????????

To install scConnect, for now, please clone this repository and install using
..code:: bash
    `cd scConnect/`
    `pip install .`

When the initial version is released we will make the packaged available at PyPI.