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

============================================================
publication -- tools for preparing astronomy papers in LaTeX 
============================================================

The :mod:`publication` module contains classes and funtions to assist in
preparing papers or other publications with an emphasis on common astronomy
journals. The tools here are for use with `LaTeX
<http://www.latex-project.org/>`_ as the actual authoring tools.

.. todo:: examples

Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.publication
   :parts: 1

Module API
----------

"""

from __future__ import division,with_statement
import re as _re

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
    from collections import Sequence,MutableSequence
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    class MutableSequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
    class Sequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary

class TeXNode(object):
    """
    An element in the TeX parsing tree.  The main shared characteristic is
    that calling the node will return a string with the combined text.
    
    *Subclassing*
    
    Subclasses must implement :meth:`getSelfText` (see docstring for details)
    
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self,parent):
        self.parent = parent
        self.children = []
    
    def getChildrenText(self):
        """
        Returns the text for the child nodes (or empty string if leaf)
        """
        txt = [c() for c in self.children]
        return ''.join(txt)
    
    @abstractmethod
    def getSelfText(self):
        """
        Subclass implementations must return a 2-tuple of strings such that
        the child text goes in between the tuple elements.
        """
        raise NotImplementedError
    
    def __call__(self):
        """
        Returns the text associated with this node, composed by combining 
        the self text and children text
        """
        st1,st2 = self.getSelfText()
        return st1 + self.getChildrenText() + st2
    
    def isLeaf(self):
        """
        Returns True if this node is a leaf (has no children)
        """
        return self.children is None or len(self.children) == 0
    
    def isRoot(self):
        """
        Returns True if this node is a root (has no parent)
        """
        return self.parent is None
    
class TeXFile(TeXNode):
    """
    TODO:DOC
    """
    def __init__(self,fn=None):
        self.parent = None
        self.children = []
        if fn is not None:
            with open(fn) as f:
                s = f.read()
            self._parse(str)
    
    def _parse(self,f):
        """
        Parse the file - `f` can be a file object or a string w/newlines
        """
        if isinstance(f,basestring):
            f = f.split('\n')
        for l in f:
            
            raise NotImplementedError
    
    def save(fn):
        with open(fn,'w') as f:
            f.write(self())
    
class TeXt(TeXNode):
    """
    A node that stores generic text. This is always a leaf.
    """
    def __init__(self,parent,text):
        self.parent = parent
        self.text = text
        
    @property
    def children(self):
        return tuple()
    
    def getSelfText(self):
        return (self.text,'')
    
class Newline(TeXt):
    """
    A node that stores just a new line. This is always a leaf.
    """
    text = '\n'
    def __init__(self,parent):
        self.parent = parent
    
class Preamble(TeXNode):
    """
    TODO:DOC
    """
    def __init__(self,parent,content):
        self.parent = parent
        self.children = []
        
        raise NotImplementedError
    
    
class Environment(TeXNode):
    """
    A LaTex environment.  
    
    *Subclassing*
    Subclasses should implement the :meth:`parse` method - see the method for 
    syntax.  They must also be registered with the :meth:`registerEnvironment` 
    static method to have them be parsed with the default :class:`TeXFile` 
    parser.
    
    """
    
    #this finds anything that begins with \begin{<name>} and ends with \end{<name}
    #group 1 is the whole env, 2 is the name, 3 is the content
    _envre = _re.compile(r'(\\begin{(.*?)}(.*?)\\end{\2})',_re.DOTALL)
    #This finds anthing like \command{...}{...}... w/ group 1 the whole command
    _commandre = _re.compile(r'(\\.*?(?:{.*})*)\W') 
    
    def __init__(self,parent,content,envname=None):
        """
        If envname is None, it will be taken from the class-level envname object
        """
        self.parent = parent
        self.children = c = []
        
        if envname is not None:
            self.envname = envname
        elif not hasattr(self,'name'):
            raise ValueError('Environment must have a name')
        
        contentl = []
        #now split out nested environments and commands and build them 
        splitenv = self._envre.split(content) 
        envstrs = splitenv[1::4]
        for i,txts in eumerate(slitenv[::4]):
            #for each text chunk, split out the command nodes
            splitcomm = _commandre.split(txts)
            for i in range(1,len(splitcomm),2):
                splitcomm[i] = command_factory(self,splitcomm[i])
            contentl.extend(splitcomm)
            if i < len(envstrs):
                contentl.append(environment_factory(self,envstrs[i]))
        
        parsed = self.parse(contentl)
        
        #now add all the nodes to the child list, 
        #transforming text into TeXt and maybe Newlines if present
        for p in parsed:
            if isinstance(p,basestring):
                c.append(Newline(self))
                for txt in p.split('\n'):
                    c.append(TeXt(self,txt))
                    c.append(Newline(self))
                
            elif isinstance(p,TeXNode):
                c.append(p)
            else:
                raise TypeError('invalid item returned from parsing %s: %s'+(self.__class__,p))
        
    def parse(self,contentlist):
        """
        This method should be overridden in subclasses to add particular 
        children or functionality.  
        
        :param contentlist: 
            A list of :class:`Command` nodes, :class:`Environment` nodes, and/or
            strings that contain the in-order content of this environment.
        :returns: 
            A list of :class:`TeXNode` objects possibly interspersed with strings. Strings will be
            automatically converted to :class:`TeXt` nodes.
        
        """
        
        return contentlist
    
    def getSelfText(self):
        b = '\\begin{'+self.envname+'}'
        e = '\\end{'+self.envname+'}'
        return (b,e)
    
    #used for static functions on the registry
    _registry = {}
    
    @staticmethod
    def registerEnvironment(envclass):
        """
        Registers the provided `envclass` in the environment registry.  Also
        returns the class to allow use as a decorator.
        
        :except TypeError: 
            If the provided class is not a :class:`Environment` subclass. 
        :except ValueError: 
            If the :attr:`envname` attribute matches one already in the registry.
        """
        if not issubclass(envclass,Environment):
            raise TypeError('envclass must be an Environment subclass')
        for e in Environment._registry:
            if envclass.envname in Environment._registry:
                raise ValueError('envname %s already present in class %s'%(envclass.envname,e))
        Environment._registry[envclass.envname] = envclass
        return envclass
    
    @staticmethod
    def unregisterEnvironment(envclass):
        """
        Removes the `envclass` :class:`Environment` object from the registered 
        environment list
        
        :param envclass: 
            The :class:`Environment` object to be removed, or its associated
            envname.
        """
        if isinstance(envclass,basestring):
            del Environment._registry[envclass]
        else:
            regclass = Environment._registry.pop(envclass.envname)
            if regclass is not envclass:
                Environment._registry[envclass.envname] = regclass
                raise KeyError("envname %s found, but doesn't match class %s"%(envclass.envname,envclass))
        
    @staticmethod
    def getEnvironments():
        """
        Returns a tuple of the registered environments.
        """
        return tuple(Environment._registry.values())
    
def environment_factory(parent,texstr):
    """
    This function takes a string from a TeX document starting with '\begin' and 
    ending in '\end{...}' and uses it to construct the appropriate Environment
    object.
    """
    enstart = texstr.index('{')+1
    enend = texstr.index('}')
    envname = texstr[enstart:enend]
    content = texstr[enend+1:textstr.rindex('\\end')]
    if envname in Environment._registry:
        envcls = Environment._registry[envname]
        return envcls(parent,content)
    else:
        return Environment(parent,content,envname)
    
@Environment.registerEnvironment
class Document(Environment):
    envname = 'document'

@Environment.registerEnvironment
class Figure(Environment):
    envname = 'figure'

class Command(TeXNode):
    """
    TODO:DOC
    """
    def __init__(self,parent,content):
        raise NotImplementedError 

def command_factory(parent,texstr):
    raise NotImplementedError 
    
    
#issue warning if abstract too long
#strip comments
#copy over bbl if necessary
#.tar.gz up with appropriate name and date            
def prep_for_arxiv_pub(texfile):
    raise NotImplementedError

#Tasks: redo figures into f##[l].eps and move the files
#set class to aastex
#strip comments
#fix any deluxetable* -> deluxetable (and rotate)
#issue warning if abstract too long
#copy over bib and bbl if necessary
#.tar.gz up with appropriate name and date
def prep_for_apj_pub(texfile,authorlast):
    raise NotImplementedError