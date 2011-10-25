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

from __future__ import division,with_statement

#Warning: this file is overwritten in the build_py stage of setup.py to freeze
#the version number - this file is only used in dev mode or during build.

#these components can be changed for new versions
major = 0
minor = 1
bugfix = None
dev = True


#everything below here is derived
release = not dev
_warntxt = ', so no revision number available for dev build of astropysics - this may prevent this install from overwriting the previous version'


def _get_bzr_devstr():
    """
    Returns the bzr devstring.
    
    :returns: The revision devstring (to be appended to the version number)
    :rtype: str
    
    :except ValueError: If this is a release version.
    """
    import os
    from os import path

    if release:
        raise ValueError('revsion devstring not valid for a release version')
    
    try:
        bzrdir = path.join(path.split(__file__)[0],'..','.bzr')
        if not path.exists(bzrdir):
            raise IOError('.bzr directory does not exist, cannot get revision')
        
        lrevfn = path.join(bzrdir,'branch','last-revision')
        if not path.exists(lrevfn):
            raise IOError('last-revision file does not exist, cannot get revision')
            
        with open(lrevfn) as f:
            s = f.read()
            
        rnum = int(s.split()[0])
        return 'dev-r'+str(rnum)
    except IOError:
        from warnings import warn
        warn('Could not find bzr'+_warntxt)
        return 'dev'

def _get_hg_devstr(idtype='rev'):
    """
    Returns the mercurial devstring.
    
    :param idtype: Can be 'rev','short',or 'hex'
    :type idtype: str
    
    :returns: The revision devstring (to be appended to the version number)
    :rtype: str
    
    :except ValueError: If this is a release version
    :except TypeError: If idtype is invalid
    """
    from os import path
    
    if release:
        raise ValueError('revsion devstring not valid for a release version')
    
    try:
        from mercurial import hg,ui
        from mercurial.node import hex,short
        from mercurial.error import RepoError
        
        try:
            basedir = path.join(path.split(__file__)[0],'..')
            repo = hg.repository(ui.ui(),basedir)
                        
            if idtype == 'rev':
                rstr = 'r'+str(repo['tip'].rev())
            elif idtype == 'short':
                rstr = 'hg-'+short(repo.lookup('tip'))
            elif idtype == 'hex':
                rstr = 'hg-'+hex(repo.lookup('tip'))
            else:
                raise TypeError('invalid idtype')
            
            return 'dev-'+rstr
            
        except RepoError:
            from warnings import warn
            warn('No mercurial repository present'+_warntxt)
            return 'dev'
        
    except ImportError:
        from warnings import warn
        warn('Could not find mercurial'+_warntxt)
        return 'dev'
    

def _get_git_devstr(sha=False):
    """
    Returns the git devstring.
    
    :param bool sha: 
        If True, the full SHA1 hash will be at the end of the devstr, otherwise,
        a revision number.
    :returns: The revision devstring (to be appended to the version number)
    :rtype: str
    
    :except ValueError: If this is a release version
    """
    from os import path
    from subprocess import Popen,PIPE
    from warnings import warn
    
    if release:
        raise ValueError('revsion devstring not valid for a release version')

    currdir = path.abspath(path.split(__file__)[0])
#    gitdir = path.join(currdir,'.git')

#    while not path.exists(gitdir):
#        currdir = path.split(currdir)[0]
#        gitdir = path.join(currdir,'.git')
#        if currdir=='/' and not path.exists(gitdir):
#            warn('No git repository present'+_warntxt)
#            return 'dev'
    try:
        p = Popen(['git','rev-list','HEAD'],cwd=currdir,
                  stdout=PIPE,stderr=PIPE,stdin=PIPE)
        stdout,stderr = p.communicate()
    except OSError:
        stdout = stderr = None
      
    if stdout is None:
        warn("git couldn't be run to determine git version number")
        return 'dev'
    elif p.returncode == 128:
        warn('No git repository present'+_warntxt)
        return 'dev'
    elif p.returncode != 0:
        warn('Git failed: '+stderr)
        return 'dev'
    
    if sha:
        return 'dev-git-'+stdout[:40]
    else:
        nrev = stdout.count('\n')
        return  'dev-r%i'%nrev
    
_get_devstr = _get_git_devstr
    
version = str(major)+'.' +\
          str(minor) +\
          (('.'+str(bugfix)) if bugfix else '') +\
          ('.'+_get_devstr() if dev else '')
          
_frozen_version_py_template = """#Autogenerated in astropysics setup.py on %s
#these components can be changed for new versions
version = '%s'

major = %s
minor = %s
bugfix = %s
dev = %s

release = not dev
"""
