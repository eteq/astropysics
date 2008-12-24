#!/usr/bin/env python

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
      url='http://www.physics.uci.edu/~etolleru/software.html',
      packages=['astropysics'],
      requires=['numpy','scipy'],
      provides=['astropysics']
     )
