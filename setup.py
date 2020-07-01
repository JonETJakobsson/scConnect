from setuptools import setup, find_packages
import glob

csv_files = glob.glob("scConnect/data/**/*.csv", recursive=True)

setup(
   name='scConnect',
   version='0.1',
   description='scConnect: Tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=["scConnect"],
   package_data={"scConnect": csv_files+["assets/*.*"]},
   
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
