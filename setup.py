from setuptools import setup

setup(
   name='scConnect',
   version='0.1',
   description='scConnect: like tinder for single cells',
   author='Jon Jakobsson',
   author_email='jon.jakobsson@neuro.uu.se',
   packages=['scConnect'],  #same as name
   install_requires=[#external packages as dependencies
      'scanpy',
      'pandas',
      'numpy',
      'networkx'], 
)