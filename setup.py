#!/usr/bin/env python

from distutils.core import setup

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
