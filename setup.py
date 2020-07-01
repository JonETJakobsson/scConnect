from setuptools import setup, find_packages

setup(
   name='scConnect',
   version='0.1',
   description='scConnect: Tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=["scConnect"],
   package_data={"scConnect": ["data/*.*", #  add all files under data/
                               "data/*/*.*", #  add all files under data/
                               "data/*/*/*.*", #  add all files under data/
                               "data/*/*/*/*.*", #  add all files under data/
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
