=============
API
=============

++++++++
Database
++++++++
.. note:: 
    The database files are constructed from csv files downloaded from  http://www.guidetopharmacology.org/download.jsp.
    Current version is py:Values:'scConnect.database.version'

.. automodule:: scConnect.database
    :members:


++++++++++
Gene calls
++++++++++
Aggregates the gene expression on the cell clusters or cell populations. 
Values returned are between 0 and 1.

.. automodule:: scConnect.genecall
   :members:

+++++++
Connect
+++++++
Find the links between gene expression and ligands and receptors. 
Also find connections between ligands and receptors between clusters.

.. automodule:: scConnect.connect
   :members:

+++++
Graph
+++++
Build a graph with edge list sand node lists. Graph can be analysed and plotted in different ways.

.. automodule:: scConnect.graph
   :members:

