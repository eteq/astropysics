#!/usr/bin/env python
#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

from  sys import argv
try:
    argv.pop(argv.index('ez_setup'))
    from ez_setup import use_setuptools
    use_setuptools()
except ValueError:
    pass

try:
    from setuptools import setup
    stls=True
except ImportError:
    from distutils.core import setup
    stls=False

setup(name='Astropysics',
      version='0.1',
      description='Astrophysics libraries for Python',
      author='Erik Tollerud',
      author_email='etolleru@uci.edu',
      licence = 'Apache License 2.0',
      url='http://www.physics.uci.edu/~etolleru/software.html#astropysics',
      requires=['numpy','scipy'],
      provides=['astropysics'],
      long_description="""
      ``astropysics`` contains a variety of utilities and algorithms for 
      reducing, analyzing, and visualizing astronomical data.
      
      while ``astropysics`` requres only ``numpy`` and ``scipy``, other 
      packages are necessary for some of the functionality.  These include: 
      ``Traits``, ``TraitsGUI``, ``PyEphem``, and ``pymc``
      """,
      packages=['astropysics','astropysics.gui'],
      package_data={'astropysics':'data/*'},
      scripts=['scripts/spylot']
     )
