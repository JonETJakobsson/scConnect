from setuptools import setup, find_packages
import glob

csv_files = glob.glob("scConnect/data/**/*.csv", recursive=True)
csv_files = [path.replace("\\", "/") for path in csv_files]
csv_files = [path.replace("scConnect/data", "data") for path in csv_files]


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
