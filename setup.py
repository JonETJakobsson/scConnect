from setuptools import setup find_packages

setup(
   name='scConnect',
   version='0.1',
   description='scConnect: like tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=find_packages(),
   package_data={"scConnect": ["data/Gene_annotation/*/*.csv", #  interaction.csv
                               "data/Gene_annotation/*/*/*.csv", #  ligands and receptors csv's
                               "data/GTP_tables/*.txt", #  licence for GTP
                               "data/GTP_tables/*/*.csv", # GTP tables for specific version
                               "data/GTP_tables/*/*.csv", # cleaned GTP tables
                               "assets/*.*" # all assets for web app
                              ]},
   
   install_requires=[#external packages as dependencies
      'scanpy',
      'pandas',
      'numpy',
      'networkx',
      'dash',
      'dash-cytoscape',
      'holoviews',
      'notebook'], 
)
