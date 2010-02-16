#!/usr/bin/env python
#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

from distribute_setup import use_setuptools
use_setuptools()

from setuptools import setup,find_packages

setup(name='Astropysics',
      version='0.1a0',
      description='Astrophysics libraries for Python',
      
      packages=find_packages(),
      package_data={'astropysics':'data/*'},
      scripts=['scripts/spylot'],
      requires=['numpy','scipy'],
      provides=['astropysics'],
      
      author='Erik Tollerud',
      author_email='etolleru@uci.edu',
      licence = 'Apache License 2.0',
      url='http://www.physics.uci.edu/~etolleru/software.html#astropysics',
      
      
      long_description="""
      ``astropysics`` contains a variety of utilities and algorithms for 
      reducing, analyzing, and visualizing astronomical data.
      
      while ``astropysics`` requres only ``numpy`` and ``scipy``, other 
      packages are necessary for some of the functionality.  These include: 
      ``matplotlib``,``Traits``, ``TraitsGUI``, ``chaco``.
      """,
     )
