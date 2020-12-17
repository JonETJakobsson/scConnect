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
 
scConnect integrate gene expression profiles in scRNA-seq datasets with ligand and receptor interaction information from `Guide to pharmacology <https://www.guidetopharmacology.org/>`__ to construct a graph containing all putative interaction between cell types in the dataset. scConnect integrate well with Scanpy and  can be appended to  any scanpy analysis pipeline.

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

* Added a method that calulates z-scores and p-values for each ligand and receptor score (:code:`connect.significance()`):
   * Utilize bootstapping to assess the random distribution of ligand and recepto score for each cell type.
   * Calculates a Z-score for each ligand and receptor score given this random distribution.
   * Calculates multiple test corrected p-values using Benjamini/Hochberg (false discovery rate) correction.
   * Estimate interaction significance by wieghting both ligand and receptor p-values 
   * :math:`I_{LR_{significance}} = -\log_{10} \frac{R_{p-value} + L_{p-value}}{2}`


   * Specify specific interactions where corrected p-value for both ligand and receptor are under 0.05.

* Updates to the web app:
   * Filter for significant interactions (where both ligand and receptor p-value are under 0.05) in the network graph.
   * Network graph filtering is propagated to the sankey graph.
   * Added a scatter plot for interaction of selected edge, where x axis is ligand z-score, y axis is receptor z-score size is log(interaction score) and color is interaction significance
   * Selection of interactions in the graph also filters the interaction table.
   * Added a scatter plot for ligands and receptors where the x axis is log(score) and y axis -log(p-value)
   * Selected ligands or receptors filters the table under the graph.



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
