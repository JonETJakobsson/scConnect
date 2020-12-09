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

* Added scConnect.connect.significance(). Calculates z-score, p-values and adjusted p-values for ligand and receptor scores based on
  the distribution of these scores when randomly assigning group lables. Must be run in order to filter the graph based on these metrices.
* Added z-score and p-value (adjusted using fales discovery rate benjamini/hochberg) to interaction list in interactive webapp.
* Integrated significance information in the web app:

   * Filter the full graph on ligands and receptors that are significanly upregulated in specific cell types

   * Filtering is propagated to the sankey graph
   
   * New ligand and receptor graph visulize score vs p-values. This makes it easier to find wither specific or highly expressed ligands and receptors.

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
