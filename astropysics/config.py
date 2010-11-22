#Copyright 2010 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""

=================================
config -- configuration and setup 
=================================

The :mod:`config` module contains functions and prompts to manage and access the
persistent astropysics configuration, as well as setup and installation actions.
It also includes utilities to set up the suggested ipython environment.

.. add this back in if any classes are added
    Classes and Inheritance Structure
    ---------------------------------

    .. inheritance-diagram:: astropysics.config
       :parts: 1

Module API
----------

"""
from __future__ import division,with_statement

#recommended packages and associated import names
_recpkgs = {'matplotlib':'matplotlib',
            'pyfits':'pyfits',
            'ipython':'IPython',
            'networkx':'networkx',
            'pygraphviz':'pygraphviz'}
_guipkgs = {'traits':'enthought.traits',
            'traitsGUI':'enthought.traits.ui.api',
            'chaco':'enthought.chaco',
            'mayavi':'enthought.mayavi'}
            
class VersionError(Exception): pass 
class BuildError(Exception): pass
class InstallError(Exception): pass 
class DownloadError(Exception): pass 

def _check_if_installed(pkgs):
    importable = {}
    for name,mod in pkgs.iteritems():
        try:
            __import__(mod)
            importable[name] = True
        except ImportError:
            importable[name] = False
    return importable

def install_package(pkgname,dldir=None,overwrite=False):
    """
    Attempt to install the package with the provided name.  
    
    :param str pkgname: 
        The name of the package to install. May include version requirement
        after name (e.g. 'numpy>1.0').
    :param bool overwrite: 
        If True, downloaded package archives will be overwritten instead of
        being re-used.
    
    :returns: True if installation was sucessful, False if not
    
    """
    import os
    
    if dldir is None:
        dldir = os.path.join(get_config_dir(),'install_pkgs')
        if not os.path.isdir(dldir):
            if os.path.exists(dldir):
                raise IOError(dldit+" is a file - can't make directory!")
            os.mkdir(dldir)
            
    try:
        fn = _download_package(pkgname,dldir,overwrite)
        _do_install(fn,dldir)
    except IOError,e:
        if 'CRC check failed' in e.args[0]:
            print 'Problem with downloaded package file',fn,'- re-downloading.'
            install_package(pkgname,dldir,True)
        else:
            raise
    
def _download_package(pkgname,dldir,overwrite=False):
    import os,urllib,xmlrpclib
    from pkg_resources import parse_version
    
    #parse out the version if it's in the package name     
    reqvers = reqstr = None
    for s in ('>=','<=','=','<','>'):
        if s in pkgname:
            pkgname = pkgname.split(s)
            reqvers = pkgname[-1]
            pkgname = pkgname[0]
            reqstr = s
                
    client = xmlrpclib.ServerProxy('http://pypi.python.org/pypi')
    vers = client.package_releases(pkgname,True)
    
    #now identify the correct version
    if reqvers is None:
        ver = vers[0]
    else:
        preqvers = parse_version(reqvers)
        if reqstr=='>=':
            if parse_version(vers[0])<preqvers:
                raise VersionError('cannot get requested version %s of package %s'%(reqvers,pkgname))
            ver = vers[0]
        elif reqstr=='<=':
            for v in vers:
                if parse_version(v)<=preqvers:
                    ver = v
                    break
            else:
                raise VersionError('cannot get requested version %s of package %s'%(reqvers,pkgname))
        elif reqstr=='=':
            for v in vers:
                if parse_version(v)<=preqvers:
                    ver = v
                    break
            else:
                raise VersionError('cannot get requested version %s of package %s'%(reqvers,pkgname))
        elif reqstr=='>':
            if parse_version(vers[0])<=preqvers:
                raise VersionError('cannot get requested version %s of package %s'%(reqvers,pkgname))
            ver = vers[0]
        elif reqstr=='<':
            for v in vers:
                if parse_version(v)<preqvers:
                    ver = v
                    break
            else:
                raise VersionError('cannot get requested version %s of package %s'%(reqvers,pkgname))
    
    
    for durl in client.release_urls(pkgname, ver):
        if durl['packagetype']=='sdist':
            url = durl['url']
            fn = durl['filename']
            break
    else:
        raise DownloadError('Could not find a source distribution for %s %s'%(pkgname,ver))
    
    fn = os.path.join(dldir,fn)
    if not overwrite and os.path.exists(fn):
        print fn,'already exists'
    else:
        print 'Downloading',url,'to',fn
        urllib.urlretrieve(url,fn)
    
    return fn

def _do_install(tgzfn,dldir):
    import tarfile,subprocess,os,sys,shutil
    from contextlib import closing
    
    print 'Untarring',tgzfn,'to',dldir
    with closing(tarfile.open(tgzfn)) as f:
        m0 = f.getmembers()[0]
        if not m0.isdir():
            idir = os.path.join(dldir,tgzfn.replace('.tgz','').replace('.tar.gz'))
            os.mkdir(idir)
            f.extractall(idir)
        else:
            idir = os.path.join(dldir,m0.name)
            f.extractall(dldir)
            
    print 'Building in',idir
    pb = subprocess.Popen(sys.executable+' setup.py build',shell=True,cwd=idir)
    bretcode = pb.wait()
    if bretcode != 0:
        raise BuildError('build of %s failed'%pkgname)
    
    print 'Installing in',idir
    pi = subprocess.Popen(sys.executable+' setup.py install',shell=True,cwd=idir)
    iretcode = pi.wait()
    if iretcode == 0:
        print '\nInstall successful, deleting install directory.',idir,'\n'
        shutil.rmtree(idir)
    else:
        raise InstallError('install of %s failed'%pkgname)
       
        

def run_install_tool():
    """
    Starts the console-based interactive installation tool.
    """
    try:
        import numpy #should be impossible to even install, but better to check
    except ImportError:
        print 'Numpy not installed - get it (http://numpy.scipy.org/) before continuing.'
    try:
        import scipy
    except ImportError:
        print 'Scipy not installed - get it (http://www.scipy.org/) before continuing.'
    
    quit = False
    while not quit:
        recs = _check_if_installed(_recpkgs)
        guis = _check_if_installed(_guipkgs)
        pkgs = recs.keys()
        insts = recs.values()
        pkgs.extend(guis.keys())
        insts.extend(guis.values())
        
        print 'Recommended packages:'
        for n,i in recs.iteritems():
            print '%i-%s: %s'%(pkgs.index(n)+1,n,'Installed' if i else 'NOT installed')
        print '\nGUI packages:'
        for n,i in guis.iteritems():
            print '%i-%s: %s'%(pkgs.index(n)+1,n,'Installed' if i else 'NOT installed')
            
        if all(insts):
            print '\nAll packages Installed - nothing left to do.\n'
            quit = True
        else:
            inpt = None
            while inpt is None:
                inpt = raw_input("\nSelect individual package to install (#),  'a' to install everything not yet installed, or 'q' to quit:")
                if inpt.strip()=='q':
                    print ''
                    quit = True
                elif inpt.strip()=='a':
                    for n,i in recs.iteritems():
                        if not i:
                            try:
                                install_package(n)
                            except Exception,e:
                                print 'Installation of',n,'Failed:',e,'Skipping...'
                else:
                    try:
                        inpt = int(inpt)-1
                    except ValueError:
                        print 'Invalid entry.'
                        inpt=None
                    try:
                        install_package(pkgs[inpt])
                    except Exception,e:
                        print 'Installation of',pkgs[inpt],'Failed:',e
    
def run_ipython_setup():
    """
    Starts the console-based ipython setup tool.
    """
    try:
        import IPython
    except ImportError:
        print 'IPython not installed - install it before running ipython setup.'
        return
    
    raise NotImplementedError

def get_config_dir(create=True):
    """
    Returns the astropysics configuration directory name.
    
    :param bool create: 
        If True, the directory will be created if it doesn't exist.
    :returns: The absolute path to the config directory as a string
    """
    import os,sys
    from os import environ as env

    #First find the home directory - this is inspired by the scheme ipython
    #uses to identify "home"
    if os.name == 'posix':
        # Linux, Unix, AIX, OS X
        if 'HOME' in env:
            homedir = env['HOME'].decode(sys.getfilesystemencoding())
        else:
            raise OSError('Could not find unix home directory to search for astropysics config dir')
    elif os.name == 'nt': # This is for all modern Windows (NT or after)
        #Try for a network home first
        if 'HOMESHARE' in env: 
            homedir = env['HOMESHARE'].decode(sys.getfilesystemencoding())
        #See if there's a local home
        elif 'HOMEDRIVE' in env and 'HOMEPATH' in env: 
            homedir = os.path.join(env['HOMEDRIVE'],env['HOMEPATH'])
            homedir = homedir.decode(sys.getfilesystemencoding())
        #Maybe a user profile?
        elif 'USERPROFILE' in env:
            homedir = os.path.join(env['USERPROFILE']).decode(sys.getfilesystemencoding())
        else:            
            try:
                import _winreg as wreg
                key = wreg.OpenKey(
                    wreg.HKEY_CURRENT_USER,
                    "Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders"
                )
                homedir = wreg.QueryValueEx(key,'Personal')[0]
                homedir = homedir.decode(sys.getfilesystemencoding())
                key.Close()
            except:
                #As a final possible resort, see if HOME is present for some reason
                if 'HOME' in env:
                    homedir = env['HOME'].decode(sys.getfilesystemencoding())
                else:
                    raise OSError('Could not find windows home directory to search for astropysics config dir')
    else:
        #for other platforms, try HOME, although it probably isn't there
        if 'HOME' in env:
            homedir = env['HOME'].decode(sys.getfilesystemencoding())
        else:
            raise OSError('Could not find a home directory to search for astropysics config dir - are you on an unspported platform?')
        
    configdir = os.path.realpath(os.path.join(homedir,'.astropysics'))
        
    if create and not os.path.isdir(configdir):
        #try to create it
        os.mkdir(configdir)
    return configdir

def get_data_dir(create=True):
    """
    Returns the directory name for data that astropysics downloads. See
    :func`astropysics.io.get_data` to work with data in this directory.
    
    :param bool create: 
        If True, the directory will be created if it doesn't exist.
    :returns: The absolute path to the data directory as a string.
    """
    import os
    datadir = os.path.join(get_config_dir(create),'data')
    if create and not os.path.isdir(datadir):
        os.mkdir(datadir)
    return datadir

