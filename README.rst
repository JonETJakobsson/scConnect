.. image:: https://readthedocs.org/projects/scconnect/badge/?version=latest
    :target: https://scconnect.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: scConnect/assets/logo.png
  :width: 100px
  :align: center
  :height: 100px
 
The end goal of scConnect is to provide information about putative interaction between cell populations 
in single cell RNA-seq data. 
To do so gene expression levels are used to asses precence or absence of ligands and receptors. 
This information is then integrated with ligand receptor interaction data from `Guide to pharmacology`__ 
to detect putative connections.
The resulting NetworkX graph can be analysed and browsed using a dash web app.

__ https://www.guidetopharmacology.org/

??????????????????????????
Documentation and tutorial
??????????????????????????
For a short introduction to how scConnect can be used please check out this `tutorial`_ and for information about the API please consult the `documentation`_.

The quickest and easiest way to try out the tutorial is to run through this binder:

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/JonETJakobsson/scConnect/master?filepath=tutorial%2FConnecting%20brain%20regions.ipynb   
 
.. _tutorial: https://github.com/JonETJakobsson/scConnect/blob/master/tutorial/Connecting%20brain%20regions.ipynb
.. _documentation: https://scconnect.readthedocs.io/en/latest/


 **Note:** Unfortunatly, to run the web application, you will have to install scConnect and run on your local machine.


????????????
Installation
????????????

To install scConnect, for now, please clone this repository and install using

.. code-block:: bash

  cd scConnect/
  pip install .

When the initial version is released we will make the package available at PyPI.

