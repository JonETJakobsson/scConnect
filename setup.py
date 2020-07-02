from setuptools import setup, find_packages
import glob

files = glob.glob("scConnect/**/*.*", recursive=True) # add all files under scConnect
files = [path.replace("\\", "/") for path in files] # Change from windows style to requred style for setuptools
files = [path.replace("scConnect/", "") for path in files]


setup(
   name='scConnect',
   version='0.1',
   description='scConnect: Tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=["scConnect"],
   package_data={
      "scConnect": files},
   
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
