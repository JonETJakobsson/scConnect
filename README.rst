.. image:: https://readthedocs.org/projects/scconnect/badge/?version=latest
    :target: https://scconnect.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
 
.. image:: https://travis-ci.com/JonETJakobsson/scConnect.svg?branch=master
    :target: https://travis-ci.com/JonETJakobsson/scConnect
    
.. image:: https://badge.fury.io/py/scConnect.svg
    :target: https://badge.fury.io/py/scConnect

.. image:: scConnect/assets/logo.png
  :width: 100px
  :align: center
  :height: 100px
 
===========================================
scConnect: Like Tinder but for single cells
===========================================

**What can I do with scConnect?**
You can investigate expression of ligands and receptors among the cell types in your scRNA-seq dataset. You can detect putative interactions between cell types which can be good starting points for further investigations *in vivo*. You can detect specific interaction between cell types, which can be good drug targets as the effect would be limited to those cell types.

**How does it work?**
scConnect estimate expression of ligands and receptors for cell types in scRNA-seq datasets. scConnect also estimate expression of molecular ligands that are synthezised by many enzymes, hence integrating gene expression related to synthesis, transport, reuptake etc. Using interaction information from `Guide to pharmacology <https://www.guidetopharmacology.org/>`__ putative cell-cell interactions can be identified. Using random permutation of cell type lables, scConnect estimate the specificity of each ligand and receptor for each cell type, and use this information to estimate the specificity of each interaction. All interactions are stored in a multi-directional graph structure and scConnect provide multitude of tools to analyse this data, including an interactive web application and several plotting functions. scConnect integrate well with Scanpy and  can be appended to  any scanpy analysis pipeline.

=========
Usecases:
=========

* Identify putative cell-cell communication in a tissue
* Infer neuronal networks based on ligand receptor compatability
* Study connectivity changes following treatment


??????????????????????????
Documentation and tutorial
??????????????????????????
For a short introduction to how scConnect can be used please check out this `tutorial`_ and for information about the API please consult the `documentation`_.

The quickest and easiest way to try out the tutorial is to run through this binder:

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/JonETJakobsson/scConnect/master?filepath=tutorial%2FConnecting%20brain%20regions.ipynb   
 
.. _tutorial: https://github.com/JonETJakobsson/scConnect/blob/master/tutorial/Connecting%20brain%20regions.ipynb
.. _documentation: https://scconnect.readthedocs.io/en/latest/




????????????
Installation
????????????

To install the latest stable version install using PYPI

.. code-block:: bash

    pip install scConnect
    
This will install all dependencies needed to run the tutorial, and utilize all features of scConnect.

To install the master branch, please clone this repository and install using

.. code-block:: bash

  cd scConnect/
  pip install .
  
  
or
 
.. code-block:: bash
 
   pip install git+https://github.com/JonETJakobsson/scConnect
    



  

===========
Change log:
===========

?????
1.0.3
?????


Major changes:

* Added a method that calulates z-scores and p-values for each ligand and receptor score (:code:`connect.specificity()`):
   * Utilize permutation to assess the random distribution of ligand and receptor score for each cell type.
   * Calculates a Z-score for each ligand and receptor score given this random distribution.
   * Calculates emperical p-values from the permutated random distribution. 
   * Calculates multiple test corrected p-values using Benjamini/Hochberg (false discovery rate) correction.
   * Estimate interaction specificity by wieghting both ligand and receptor p-values .

* Updates to the web app:
   * Summmize and filter edges based on specificity in the network graph.
   * Added possibility to download current network graph view as a svg file.
   * Filter based on specificity in sankey graph.
   * Added a scatter plot for interaction of selected edge, where x axis is log(interaction score), y axis is specificity and color is interaction importance.
   * Selection of interactions in the graph also filters the interaction table.
   * Added a scatter plot for ligands and receptors where the x axis is log(score) and y axis -log(p-value)
   * Selected ligands or receptors filters the table under the graph.

* retrieving data from graph:
   * Retrieve interaction data using :code:`graph.edge_list()` and plot a dotplot using :code:`graph.dotplot()`
   * Retrieve information about ligands and receptors using :code:`graph.get_ligand_df()` and :code:`graph.get_receptor_df()`

* Save progress
    * Save calculated specificity using :code:`connect.save_specificity()` and :code:`connect.load_specificity()`.

Minor Changes:

* Updated GTP database to 2020-5 from 2019-5.


?????
1.0.2
?????

* Fixed documentation issues (added .readthedocs.yml)
* removed requirement.txt, build is constructed entirely from setup.py

?????
1.0.1
?????

Bugfixes:

* Fixed a bug in connect.py which cased a crash when connecting ligands and receptors.


?????
1.0.0
?????

Initial release.
