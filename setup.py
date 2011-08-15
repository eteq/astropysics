#!/usr/bin/env python
#Copyright (c) 2008 Erik Tollerud (erik.tollerud@gmail.com) 
from __future__ import division,with_statement

from glob import glob
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup,find_packages
from distutils.command.build_py import build_py as du_build_py
from distutils.core import Command

from astropysics.version import version as versionstr
from astropysics.config import _recpkgs,_guipkgs

    
descrip = """
`astropysics` contains a variety of utilities and algorithms for reducing, analyzing, and visualizing astronomical data.
      
See http://packages.python.org/Astropysics/ for detailed documentation.
"""

apyspkgs = find_packages(exclude=['tests'])
scripts = glob('scripts/*')

#recommended/gui packages are stored in config module - used in extras
recpkgs = [pkg.name for pkg in _recpkgs]
guipkgs = [pkg.name for pkg in _guipkgs]


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
        
        
#custom sphinx builder just makes the directory to build if it hasn't already been made
try:
    from sphinx.setup_command import BuildDoc
    
    class apy_build_sphinx(BuildDoc):
        def finalize_options(self):
            from os.path import isfile    
            from distutils.cmd import DistutilsOptionError
            
            if self.build_dir is not None:
                if isfile(self.build_dir):
                    raise DistutilsOptionError('Attempted to build_sphinx into a file '+self.build_dir)
                self.mkpath(self.build_dir)
            return BuildDoc.finalize_options(self)
            
except ImportError: #sphinx not present
    apy_build_sphinx = None
    
    
#command to count the number of lines of code (mostly for curiosity's sake) in the main dirs
class CountLines(Command):
    # Brief (40-50 characters) description of the command
    description = "Print the number of lines in the major directories to the terminal."

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('includeempty', 'e',
                     "Include empty lines in the count"),
                   ]
    
    def initialize_options (self):
        self.includeempty = False

    def finalize_options (self):
        pass
    
    def visit_files(self,lists,dirname,fnames):
        lcountlist,fcountlist = lists
        from os import path
        
        #prefilter for valid extentions
        if dirname != 'scripts':
            fnames = [fn for fn in fnames if (fn.endswith('.py') or fn.endswith('.pyx')) ]
        
        cnt = 0
        for fn in fnames:
            fn = path.join(dirname,fn)
            with open(fn) as f:
                if self.includeempty:
                    for l in f:
                        cnt += 1
                else:
                    for l in f:
                        if l.strip()!='':
                            cnt += 1
        
        lcountlist.append(cnt)
        fcountlist.append(len(fnames))

    def run(self):
        from os import path
        
        dir,name = path.split(__file__)
        
        apydir = path.join(dir,'astropysics')
        apyllst,apyflst = [],[]
        path.walk(apydir,self.visit_files,(apyllst,apyflst))
        self.apylinecount = sum(apyllst)
        self.apyfilecount = sum(apyflst)
        
        scrdir = path.join(dir,'scripts')
        scrllst,scrflst = [],[]
        path.walk(scrdir,self.visit_files,(scrllst,scrflst))
        self.scrlinecount = sum(scrllst)
        self.scrfilecount = sum(scrflst)
        
        tstdir = path.join(dir,'tests')
        tstllst,tstflst = [],[]
        path.walk(tstdir,self.visit_files,(tstllst,tstflst))
        self.tstlinecount = sum(tstllst)
        self.tstfilecount = sum(tstflst)
        
        self.linecount = self.apylinecount + self.scrlinecount + self.tstlinecount
        self.filecount =  self.apyfilecount + self.scrfilecount + self.tstfilecount
        
        print 'Astropysics source directory has %i lines in %i files'%(self.apylinecount,self.apyfilecount)
        print 'Scripts directory has %i lines in %i files'%(self.scrlinecount,self.scrfilecount)
        print 'Tests directory has %i lines in %i files'%(self.tstlinecount,self.tstfilecount)
        print 'Total %i lines in %i files'%(self.linecount,self.filecount)
         

cmdclassd = {'build_py' : apy_build_py,'count_lines':CountLines}
if apy_build_sphinx is not None:
    cmdclassd['build_sphinx'] = apy_build_sphinx
    
setup(name='Astropysics',
      version=versionstr,
      description='Astrophysics libraries for Python',
      
      packages=apyspkgs,
      package_data={'astropysics':['data/*']},
      scripts=scripts,
      requires=['numpy','scipy'],
      install_requires=['numpy'],
      provides=['astropysics'],
      extras_require={'all':recpkgs+guipkgs,
                      'nogui':recpkgs},  
      author='Erik Tollerud',
      author_email='erik.tolleru@gmail.com',
      license = 'Apache License 2.0',
      url='http://packages.python.org/Astropysics/',
      long_description=descrip,
      cmdclass = cmdclassd
     )
