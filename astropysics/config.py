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

The :mod:`config` module contains functions to manage and access the persistent
astropysics configuration. It also includes utilities to install recommended
packages and set up the ipython environment.

.. add this back in if any class structure is added
    Classes and Inheritance Structure
    ---------------------------------

    .. inheritance-diagram:: astropysics.config
       :parts: 1

Module API
----------


"""
from __future__ import division,with_statement

class VersionError(Exception): pass
class InstallError(Exception): pass 
class DownloadError(Exception): pass 

#packages avaiable for install and the associated objects
class PackageInstaller(object):
    """
    Represents a package to be downloaded and installed.
    """
    def __init__(self,name,importmod=None):
        """
        :param name: The name of the package.
        :param importmod: 
            The module name to import to test if the pacakge is installed. If
            None, will be assumed to match `name`
        
        """
        self.name = name
        if importmod is None:
            self.importmod = name
        else:
            self.importmod = importmod
        
    def getUrl(self):
        """
        Override this to get a URL from somewhere other than PyPI.
        
        :returns: 
            (url,fn) where `url` is the URL to the source distribution, and `fn`
            is the filename to use locally if None, the end of the URL will be
            used)
            
        """
        return self._getUrlfromPyPI()
    
    def _getUrlfromPyPI(self):
        import os,xmlrpclib
        from pkg_resources import parse_version
        
        #TODO: version string in object
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
        
        return url,fn
    
    def isInstalled(self):
        """
        Test if the package is installed.
        
        :returns: True if the module is installed, otherwise False.
        """
        try:
            __import__(self.importmod)
            return True
        except ImportError:
            return False
        
    def download(self,dldir=None,overwrite=False):
        """
        Attempt to install the package with the provided name.  
        
        :param str dldir: 
            The directory to download to.  If None, the standard configuration
            directory will be used.
        :param bool overwrite: 
            If True, downloaded package archives will be overwritten instead of
            being re-used.
            
        :returns: The full path to the downloaded file.
        
        :raises DownloadError: If the pacakge could not be located
        """        
        import os,urlparse,urllib
        
        pkgname = self.name
        if dldir is None:
            dldir = os.path.join(get_config_dir(),'install_pkgs')
            if not os.path.isdir(dldir):
                os.mkdir(dldir)
        
        url,fn = self.getUrl()
        if fn is None:
            fn = os.path.split(urlparse.urlparse(url).path)[-1]
        
        dlfn = os.path.join(dldir,fn)
        if not overwrite and os.path.exists(dlfn):
            print dlfn,'already exists - not overwriting.'
        else:
            print 'Downloading',url,'to',dlfn
            urllib.urlretrieve(url,dlfn)
        
        return dlfn
        
        
        
    def install(self,dldir=None,overwrite=False):
        
        """
        Download the package, if necessary, and install it.
        
        param str dldir: 
            The directory to download to.  If None, the standard configuration
            directory will be used.
        :param bool overwrite: 
            If True, downloaded package archives will be overwritten instead of
            being re-used.
        
        :raises InstallError: If the install fails
        """
        import tarfile,subprocess,os,sys,shutil
        from contextlib import closing
                
        fn = self.download(dldir,overwrite)
        dldir,tgzfn = os.path.split(fn)[0]
        
        try:
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
        except IOError,e:
            if 'CRC check failed' in e.args[0]:
                print 'Problem with downloaded package file',fn,'- re-downloading.'
                install_package(pkgname,dldir,True)
            else:
                raise
            
class _MatplotlibInstaller(PackageInstaller):
    def __init__(self):
        PackageInstaller.__init__(self,'matplotlib')
        
class _PyfitsInstaller(PackageInstaller):
    def __init__(self):
        PackageInstaller.__init__(self,'pyfits')
    

_recpkgs = [PackageInstaller('ipython','IPython'),
            _MatplotlibInstaller(),
            _PyfitsInstaller(),
            PackageInstaller('networkx'),
            PackageInstaller('pygraphviz')]
_guipkgs = [PackageInstaller('traits','enthought.traits'),
            PackageInstaller('traitsUI','enthought.traits.ui.api'),
            PackageInstaller('chaco','enthought.chaco'),
            PackageInstaller('mayavi','enthought.mayavi')]
            
            
def run_install_tool():
    """
    Starts the console-based interactive installation tool.
    """
    #TODO: remork guts for PackageInstaller class
    
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
                    print '\n'
                else:
                    try:
                        inpt = int(inpt)-1
                    except ValueError:
                        print 'Invalid entry.'
                        inpt=None
                    try:
                        install_package(pkgs[inpt])
                    except Exception,e:
                        print 'Installation of',pkgs[inpt],'Failed:',e,'\n\n'
    
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

