from setuptools import setup, find_packages
import glob

files = glob.glob("scConnect/**/*.*", recursive=True) # add all files under scConnect
files = [path.replace("\\", "/") for path in files] # Change from windows style to requred style for setuptools
files = [path.replace("scConnect/", "") for path in files]

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
   name='scConnect',
   version='1.0.3',
   description='scConnect integrate gene expression profiles in scRNA-seq datasets with ligand and receptor interaction information from Guide to pharmacology to construct a graph containing all putative interaction between cell types in the dataset.',
   long_description=long_description,
   author='Jon E.T. Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   Project_URL={"GitHub": "https://github.com/JonETJakobsson/scConnect",
                "Documentation": "https://scconnect.readthedocs.io/en/latest/"},
   license='License :: OSI Approved :: MIT License',
   key_words=['scRNA-seq, connectome, ligands, receptors, interactions, sequencing'],
   packages=["scConnect"],
   package_data={
      "scConnect": files},
   
   install_requires=[#external packages as dependencies
      'scanpy>=1.4',
      'pandas>=1.0',
      'numpy',
      'networkx',
      'dash',
      'dash-cytoscape',
      'jupyter-dash',
      'jupyter-server-proxy',
      'ipywidgets',
      'holoviews',
      'notebook',
      'xlrd'],
   Requires_Python='>=3.6', #utilize ordered dictionaries
   
)
