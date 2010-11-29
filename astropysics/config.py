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
from HTMLParser import HTMLParser as _HTMLParser

class VersionError(Exception): pass
class InstallError(Exception): pass 
class DownloadError(Exception): pass 

#packages avaiable for install and the associated objects
class PackageInstaller(object):
    """
    Represents a python package to be downloaded and installed.
    """
    def __init__(self,name,importmod=None,version=None,buildargs='',instargs='',
                      extrainfo=None,verbose=True):
        """
        :param name: The name of the package.
        :param importmod: 
            The module name to import to test if the pacakge is installed. If
            None, will be assumed to match `name`
        :param str version: 
            A version request for finding the package on PyPI, such as
            '0.2','>0.1' (greater than 0.1), or '<0.3'. Can also be None to get
            the most recent version. If the PyPI entry only has a download link,
            this is ignored.
        :param buildargs: 
            Arguments to be given to the "python setup.py build [args]" stage.
            Can be either a string or a sequence of strings.
        :param instargs:
            :param buildargs: 
            Arguments to be given to the "python setup.py install [args]" stage.
            Can be either a string or a sequence of strings.
        :param str extrainfo:
            A string with additional information about the pacakge (to be shown
            if the user requests it in the install tool). Can be None to
            indicate no extra info.
        :param bool verbose:
            If True, information will be printed to standard out about steps in
            the install process.
            
        *Subclassing*
        
        If a package needs some additional install steps, the
        :meth:`postInstall` and :meth:`preInstall` methods can be overridden
        (default does nothing). If the package isn't in PyPI, the :meth:`getURL`
        method should be overridden to return the necessary URL.
        
        """
        self.name = name
        if importmod is None:
            self.importmod = name
        else:
            self.importmod = importmod
        self.version = version
        self.buildargs = buildargs
        self.instargs = instargs
        self.extrainfo = extrainfo
        self.verbose = verbose
        
    def getUrl(self):
        """
        Returns the URL to download to get this package. Override this to get a
        URL from somewhere other than PyPI.
        
        :returns: 
            (url,fn) where `url` is the URL to the source distribution, and `fn`
            is the filename to use locally if None, the end of the URL will be
            used)
            
        """
        return self._getUrlfromPyPI()
    
    def _getUrlfromPyPI(self):
        import os,xmlrpclib,urllib2
        from pkg_resources import parse_version
        from urlparse import urlparse
        from contextlib import closing
        
        pkgname = self.name
        reqvers = self.version
        
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
            rd = client.release_data(pkgname,ver)
            if 'download_url' in rd:
                finder = _DownloadURLFinder()
                with closing(urllib2.urlopen(rd['download_url'])) as uf:
                    finder.feed(uf.read())
                finder.close()
                
                #find *first* plausible download link and assume that's it.
                for url in finder.urls:
                    if any([ext in url for ext in ('.tar.gz','.tgz','.zip')]):
                        break
                else:
                    raise DownloadError('download URL %s is not a source distribution'%url)
                
                #find the element in the path with the correct extension.  Assume
                #that's the correct download name.
                for pathpiece in urlparse(url).path.split('/'):
                    if any([pathpiece.endswith(ext) for ext in ('.tar.gz','.tgz','.zip')]):
                        fn = pathpiece
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
        Attempt to download this package  
        
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
            if self.verbose:
                print dlfn,'already exists - not overwriting.'
        else:
            if self.verbose:
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
        
        :raises InstallError: 
            If the install fails (:meth:`postInstall` will be called immediately
            before).
        """
        import tarfile,zipfile,subprocess,os,sys,shutil
        from contextlib import closing
                
        fn = self.download(dldir,overwrite)
        dldir,cfn = os.path.split(fn)
        
        try:
            if cfn.endswith('.tar.gz') or cfn.endswith('.tgz'):
                if self.verbose:
                    print 'Untarring',cfn,'to',dldir
                cfntype = tarfile
                cfnbase = cfn.replace('.tgz','').replace('.tar.gz','')
            elif cfn.endswith('.zip'):
                if self.verbose:
                    print 'Unzipping',cfn,'to',dldir
                cfntype = zipfile
                cfnbase = cfn.replace('.zip','')
            else:
                raise InstallError('Downloaded file %s is not a zip or tar source archive'%cfn)
                
            with closing(cfntype.open(fn)) as f:
                m0 = f.getmembers()[0]
                if not m0.isdir():
                    idir = os.path.join(dldir,cfnbase)
                    os.mkdir(idir)
                    f.extractall(idir)
                else:
                    idir = os.path.join(dldir,m0.name)
                    f.extractall(dldir)
                    
            self.preInstall(idir)
                    
            if self.verbose:
                print 'Building in',idir
            if not isinstance(self.buildargs,basestring):
                buildargs = ' '.join(self.buildargs)
            else:
                buildargs = ' ' + self.buildargs
            pb = subprocess.Popen(sys.executable+' setup.py build'+buildargs,shell=True,cwd=idir)
            bretcode = pb.wait()
            if bretcode != 0:
                raise BuildError('build of %s failed'%pkgname)
            
            if self.verbose:
                print 'Installing in',idir
            if not isinstance(self.instargs,basestring):
                instargs = ' '.join(self.instargs)
            else:
                instargs = ' ' + self.instargs
            pi = subprocess.Popen(sys.executable+' setup.py install'+instargs,shell=True,cwd=idir)
            iretcode = pi.wait()
            
            if iretcode == 0:
                self.postInstall(idir,True)
                if self.verbose:
                    print '\nInstall successful, deleting install directory.',idir,'\n'
                shutil.rmtree(idir)
            else:
                self.postInstall(idir,False)
                raise InstallError('install of %s failed'%pkgname)
            
        except IOError,e:
            if 'CRC check failed' in e.args[1]:
                if self.verbose:
                    print 'Problem with downloaded package file',fn,'- re-downloading.'
                self.install(dldir,True)
            else:
                raise
            
    def preInstall(self,idir):
        """
        Subclasses can override this method to do something before building and
        installing occurs.
        
        :meth str idir: The path to the directory in which the package is built.
        
        """
        pass
    
    def postInstall(self,idir,success):
        """
        Subclasses can override this method to do something after building and
        installing occurs. Only called if install succeeds.
        
        :meth str idir: The path to the directory in which the package is built.
        :meth bool success: 
            If True, the install was sucessful. Otherwise, it failed.
        
        """
        pass
    
class _DownloadURLFinder(_HTMLParser):
    """
    Used in _getUrlfromPyPI to find source download URLs.
    """
    def __init__(self):
        _HTMLParser.__init__(self)
        self.urls = []
    
    def handle_starttag(self,tag, attrs):
        if tag.lower()=='a':
            for name,val in attrs:
                if name.lower()=='href':
                    self.urls.append(val)
        
class _PyfitsInstaller(PackageInstaller,_HTMLParser):
    def __init__(self):
        extrainfo = 'Requires a C-compiler to install.'
        PackageInstaller.__init__(self,'pyfits',extrainfo=extrainfo)
        _HTMLParser.__init__(self)
        self.intr = False
        self.lastlink = None
        self.dlurl = None
        
    def getUrl(self):
        import urllib2,urlparse,os
        dlurl = 'http://www.stsci.edu/resources/software_hardware/pyfits/Download'
        uf = urllib2.urlopen(dlurl)
        try:
            self.feed(uf.read())
        finally:
            uf.close()
            
        url = self.dlurl
        fn= os.path.split(urlparse.urlparse(url).path)[-1]
            
        return url,fn
    
    #HTMLParser methods
    def handle_starttag(self,tag,attrs):
        if tag.lower()=='tr':
            self.intr = True
        elif tag.lower()=='a':
            for n,v in attrs:
                if n.lower()=='href':
                    self.lastlink = v
    def handle_endtag(self,tag):
        if tag=='tr':
            self.intr = False
    def handle_data(self,data):
        if self.intr and 'Current stable release' in data:
            self.dlurl = self.lastlink
        
  
#<-------------------Recommended Packages-------------------------------------->

_recpkgs = [PackageInstaller('ipython','IPython'),
            PackageInstaller('matplotlib',extrainfo='Requires a C-compiler to install.'),
            _PyfitsInstaller(),
            PackageInstaller('networkx'),
            PackageInstaller('pygraphviz')]
_guipkgs = [PackageInstaller('traits','enthought.traits'),
            PackageInstaller('traitsUI','enthought.traits.ui.api',extrainfo='Requires WxWidgets or Qt to be installed'),
            PackageInstaller('chaco','enthought.chaco'),
            PackageInstaller('mayavi','enthought.mayavi',extrainfo='Requires VTK to be installed.')]
            
            
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
        pkgs = []
        insts = []
        
        print 'Recommended packages:'
        i = 1
        for pkg in _recpkgs:
            pkgs.append(pkg)
            n = pkg.name
            infostr = '' if pkg.extrainfo is None else ' (i)'
            insts.append(pkg.isInstalled())
            print '%i-%s%s: %s'%(i,n,infostr,'Installed' if insts[-1] else 'NOT installed')
            i += 1
        print '\nGUI packages:'
        for pkg in _guipkgs:
            pkgs.append(pkg)
            n = pkg.name
            infostr = '' if pkg.extrainfo is None else ' (i)'
            insts.append(pkg.isInstalled())
            print '%i-%s%s: %s'%(i,n,infostr,'Installed' if insts[-1] else 'NOT installed')
            i += 1
            
        if all(insts):
            print '\nAll packages Installed - nothing left to do.\n'
            quit = True
        else:
            inpt = None
            while inpt is None:
                inpt = raw_input("\nSelect individual package to install (#),  'a' to install everything not yet installed, 'i#' for information about a package, or 'q' to quit installer:")
                linpt = inpt.strip().lower()
                if linpt=='q':
                    print ''
                    quit = True
                elif linpt=='a':
                    for pkg in pkgs:
                        if not pkg.isInstalled():
                            try:
                                pkg.install()
                            except Exception,e:
                                print 'Installation of',pkg.name,'Failed:',e,'Skipping...'
                    print '\n'
                elif linpt.startswith('i'):
                    try:
                        i = int(linpt[1:])-1
                    except ValueError:
                        print 'Invalid package number.\n'
                        continue
                    
                    pkg = pkgs[i]
                    if pkg.extrainfo is None:
                        print '\nNo extra information for package',pkg.name,'\n'
                    else:
                        print '\nExtra info for package',pkg.name+':\n',pkg.extrainfo,'\n'
                        raw_input('(Press enter to continue)')
                else:
                    try:
                        i = int(linpt)-1
                    except ValueError:
                        print 'Invalid package number.\n'
                        continue
                    
                    try:
                        pkgs[i].install()
                    except Exception,e:
                        print 'Installation of',pkgs[i].name,'Failed:',e.__class__.__name__,e,'\n\n'
                        
    
def run_ipython_setup():
    """
    Runs the console-based ipython setup tool.
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
