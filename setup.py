#!/usr/bin/env python
#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 



from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup,find_packages
from distutils.command.build_py import build_py as du_build_py

from astropysics.version import version as versionstr

#custom build_py overwrites version.py with a version overwriting the revno-generating version.py
class apy_build_py(du_build_py):
    def run(self):
        from os import path
        res = du_build_py.run(self)
        
        versfile = path.join(self.build_lib,'astropysics','version.py')
        print 'freezing version number to',versfile
        with open(versfile,'w') as f: #this overwrites the actual version.py
            f.write(self.get_version_py())
        
        return res
        
    def get_version_py(self):
        import datetime
        from astropysics.version import _frozen_version_py_template
        from astropysics.version import version,major,minor,bugfix,dev
        
        
        timestamp = str(datetime.datetime.now())
        t = (timestamp,version,major,minor,bugfix,dev)
        return _frozen_version_py_template%t

setup(name='Astropysics',
      version=versionstr,
      description='Astrophysics libraries for Python',
      
      packages=find_packages(),
      package_data={'astropysics':['data/*']},
      scripts=['scripts/spylot','scripts/fitsinfo'],
      install_requires=['numpy','scipy'],
      provides=['astropysics'],
      extras_require={'plots':'matplotlib',
                      'guis':['traits','traitsGUI','chaco'],
                      'gui3d':'mayavi',
                      'fits':'pyfits',
                      'all':['matplotlib','traits','traitsGUI','chaco','pyfits','mayavi']},
      
      author='Erik Tollerud',
      author_email='etolleru@uci.edu',
      license = 'Apache License 2.0',
      url='http://www.physics.uci.edu/~etolleru/software.html#astropysics',
      long_description="""
      ``astropysics`` contains a variety of utilities and algorithms for 
      reducing, analyzing, and visualizing astronomical data.
      
      while ``astropysics`` requres only ``numpy`` and ``scipy``, other 
      packages are necessary for some of the functionality.  These include: 
      ``matplotlib``,``Traits``, ``TraitsGUI``, ``chaco``.
      """,
      cmdclass = {'build_py' : apy_build_py}
     )
