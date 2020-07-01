from setuptools import setup

setup(
   name='scConnect',
   version='0.1',
   description='scConnect: like tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=['scConnect'],  #same as name
   package_dir={"scConnect": "scConnect"},
   package_data={"scConnect": ["data/*", "assets/*"]},
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
