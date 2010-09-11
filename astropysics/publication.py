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

#TODO:prep functions

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
    
    @abstractmethod
    def getSelfText(self):
        """
        Subclass implementations must return a 2-tuple of strings such that the
        child text goes in between the tuple elements. Alternatively, it can
        return a 3-tuple (before,between,after), and the resulting text will be
        "<beforetext><childtext1><between><childtext2>...<after>". It can also
        be None, in which case just the strings from the children will be
        returned. Otherwise, it can return a string, which will be returned as
        the full text.
        """
        raise NotImplementedError
    
    def __call__(self):
        """
        Returns the text associated with this node, composed by combining 
        the self text and children text (or just the self text if 
        :meth:`getSelfText` returns a string)
        """
        st = self.getSelfText()
        if st is None:
            return ''.join([c() for c in self.children])
        elif isinstance(st,basestring):
            return st
        elif len(st) == 2:
            st1,st2 = st
            return st1 + ''.join([c() for c in self.children]) + st2
        elif len(st) == 3:
            st1,btwn,st2 = st
            return st1 + btwn.join([c() for c in self.children]) + st2
        else:
            raise TypeError('Self text for node '+str(self)+' is invalid length')
        
        
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
    A TeX Document loaded from a file.
    """
    def __init__(self,fn=None):
        self.parent = None
        self.children = []
        if fn is not None:
            with open(fn) as f:
                s = f.read()
            self._parse(s)
    
    def _parse(self,f):
        """
        Parse the file - `f` can be a file object or a string w/newlines
        """
        if isinstance(f,basestring):
            f = f.split('\n')
        
        preamblecontent = []
        content = None
        for l in f:
            if content is None:
                if r'\begin{document}' in l:
                    content = [l]
                else:
                    preamblecontent.append(l)
            else:
                content.append(l)
        preamblestr = '\n'.join(preamblecontent)+'\n'
        contentstr = '\n'.join(content)+'\n'
        
        self.preamble = Preamble(self,preamblestr)
        self.children = _text_to_nodes(self,contentstr)
        self.children.insert(0,self.preamble)
    def getSelfText(self):
        return None
        
    def save(self,fn):
        """
        Save the content of this object to a new file.
        
        :param str fn: The name of the file to save.
        """
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
    
#<-----------------Useful regexes for parsing documents------------------------>
#this finds anything that begins with \begin{<name>} and ends with \end{<name}
#group 1 is the whole env, 2 is the name, 3 is the content
_envre = _re.compile(r'(\\begin{(.*?)}(.*?)\\end{\2})',_re.DOTALL)
#This finds anthing like \command{...}{...}... w/ group 1  is the command name,
# 2 the arguments
_commandre = _re.compile(r'\\(.*?)((?:(?:(?:{.*?})|(?:\[.*?\]))+)|\s)')
#_commandre = _re.compile(r'\\(.*?)(?:(?:((?:\[.*?\])*)((?:{.*?})+))|\s)') # group 1  is the command name, 2 the optional arguments, and 3 the required arguments (including all braces)
#_commandre = _re.compile(r'\\(.*?)((?:{.*})|\W)')
#_commandre = _re.compile(r'\\(.*?)((?:{.*})|(?:\[.*})|\s)')
class Environment(TeXNode):
    """
    A LaTex environment.  
    
    *Subclassing*
    Subclasses must implement the :meth:`parse` method - see the method for 
    syntax.  They should also be registered with the :meth:`registerEnvironment` 
    static method to have them be parsed with the default :class:`TeXFile` 
    parser. Generally, they should also have a class attribute named `envname`
    that gives the name of the environment (this name will be automatically used
    to determine which environments the subclass represents)
    """
    
    def __init__(self,parent,content,envname=None):
        """
        If envname is None, it will be taken from the class-level envname object
        """
        self.parent = parent
        self.children = c = []
        
        if envname is not None:
            self.name = envname
        elif not hasattr(self,'name'):
            raise ValueError('Environment must have a name')
        
        self.children = self.postParse(_text_to_nodes(self,content))
        
    def postParse(self,nodes):
        return nodes
    
    def getSelfText(self):
        b = '\\begin{'+self.name+'}'
        e = '\\end{'+self.name+'}'
        return (b,e)
    
    #used for static functions on the registry
    _registry = {}
    
    @staticmethod
    def registerEnvironment(envclass):
        """
        Registers the provided `envclass` in the environment registry.  Also
        returns the class to allow use as a decorator.
        
        :param envclass: 
            The :class:`Environment` object to be registered.
        
        :except TypeError: 
            If the provided class is not a :class:`Environment` subclass. 
        :except ValueError: 
            If the :attr:`name` attribute of `envclass` matches one already in 
            the registry.
        """
        if not issubclass(envclass,Environment):
            raise TypeError('envclass must be an Environment subclass')
        for e in Environment._registry:
            if envclass.name in Environment._registry:
                raise ValueError('Environment name %s already present in class %s'%(envclass.name,e))
        Environment._registry[envclass.name] = envclass
        return envclass
    
    @staticmethod
    def unregisterEnvironment(envclass):
        """
        Removes the `envclass` :class:`Environment` object from the registered 
        environment list
        
        :param envclass: 
            The :class:`Environment` object to be removed, or its associated
            name.
        """
        if isinstance(envclass,basestring):
            del Environment._registry[envclass]
        else:
            regclass = Environment._registry.pop(envclass.name)
            if regclass is not envclass:
                Environment._registry[envclass.name] = regclass
                raise KeyError("Environment name %s found, but doesn't match class %s"%(envclass.name,envclass))
        
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
    content = texstr[enend+1:texstr.rindex('\\end')]
    if envname in Environment._registry:
        envcls = Environment._registry[envname]
        return envcls(parent,content)
    else:
        return Environment(parent,content,envname)
    
@Environment.registerEnvironment
class Document(Environment):
    name = 'document'

@Environment.registerEnvironment
class Figure(Environment):
    name = 'figure'

class RequiredArgument(TeXNode):
    """
    An argument to a macro that is required (i.e. enclosed in curly braces)
    """
    children = None #arguments are always a leaf
    def __init__(self,parent,text):
        self.parent = parent
        self.text = text
        
    def getSelfText(self):
        return self.text
    
class OptionalArgument(TeXNode):
    """
    An argument to a macro that is required (i.e. enclosed in square brackets)
    """
    children = None #arguments are always a leaf
    def __init__(self,parent,text):
        self.parent = parent
        self.text = text
        
    def getSelfText(self):
        return self.text

_cmdargsepre = _re.compile(r'}|\]')
class Command(TeXNode):
    """
    A LaTeX command (anything with leading backslash and possible arguments
    that isn't an environment)
    """
    def __init__(self,parent,content):
        """
        Content can either be a string with the command text, or a (name,args)
        tuple where args is a sequence of strings or a string 
        '{arg1}[oarg1]{arg2}...'
        """
        if isinstance(content,basestring):
            if not content.startswith('\\'):
                raise ValueError("Tried to make a Command that doesn't start with a backslash")
            curlyind = content.find('{')
            squareind = content.find('[')
            if curlyind==-1:
                if squareind==-1:
                    cmdname = content[1:]
                    args = ''
                else:
                    cmdname = content[1:squareind]
                    args = content[squareind:]
            elif squareind==-1:
                cmdname = content[1:curlyind]
                args = content[curlyind:]
            else:
                ind = min(curlyind,squareind)
                cmdname = content[1:ind]
                args = content[ind:]
        elif len(content)==2:
            cmdname,args = content
        else:
            raise ValueError('Invalid content passed to Command')
        
        self.parent = parent 
        self.name = cmdname
        
        if args.strip() == '':
            self.children = []
        else:
            argspl = _cmdargsepre.split(args)
            argnodes = []
            #last one should always be empty string
            assert argspl[-1] == ''
            argspl = argspl[:-1]
            for a in argspl:
                if a[0]=='{':
                    argnodes.append(RequiredArgument(self,a[1:]))
                elif a[0]=='[':
                    argnodes.append(OptionalArgument(self,a[1:]))
                else:
                    raise ValueError('invalid set of arguments '+args)
            
            self.children = argnodes
        
    @property
    def reqargs(self):
        """
        A list of strings with the text of the required arguments (arguments 
        enclosed in curly braces).
        """
        return [c.text for c in self.children if isinstance(c,RequiredArgument)]
        
    @property
    def optargs(self):
        """
        A list of strings with the text of the optional arguments (arguments 
        enclosed in square brackets)
        """
        return [c.text for c in self.children if isinstance(c,OptionalArgument)]
        
    def getSelfText(self):
        args = ['\\',self.name]
        for c in self.children:
            if isinstance(c,RequiredArgument):
                args.append('{'+c.text+'}')
            elif isinstance(c,OptionalArgument):
                args.append('['+c.text+']')
            else:
                raise TypeError('Command children should be RequiredArguments or OptionalArguments')
        return ''.join(args)
    
    def __str__(self):
        return 'Command@%i:%s'%(id(self),self.name)
    
class Preamble(TeXNode):
    """
    The preamble of a TeX File (i.e. everything before \begin{document} )
    """
    def __init__(self,parent,content):
        self.parent = parent
        self.children = _text_to_nodes(self,content)
        
        self.docclass = None
        self.packages = []
        for c in self.children:
            if isinstance(c,Command):
                if c.name == 'documentclass':
                    self.docclass = c.reqargs[0]
                if c.name == 'usepackage':
                    self.packages.append(c.reqargs[0])
                
    def getSelfText(self):
        return None

def _text_to_nodes(parent,txt):
    """
    Converts a string into a list of corresponding TeXNodes, splitting out
    commands, environments, and Newlines.
    """
    txtl = []
    #now split out nested environments and commands and build them 
    splitenv = _envre.split(txt) 
    envstrs = splitenv[1::4]
    for i,txts in enumerate(splitenv[::4]):
        #for each text chunk, split out the command nodes and then add
        #the chunks back to the contentlist
        splitcomm = _commandre.split(txts)
        for j in reversed(range(1,len(splitcomm),3)):
            splitcomm[j] = Command(parent,(splitcomm[j],splitcomm[j+1]))
            del splitcomm[j+1]
        txtl.extend(splitcomm)
        if i < len(envstrs):
            txtl.append(environment_factory(parent,envstrs[i]))
    
    #now add all the nodes to the child list, 
    #transforming text into TeXt and maybe Newlines if present
    contentl = []
    for t in txtl:
        if isinstance(t,basestring):
            for txt in t.split('\n'):
                contentl.append(TeXt(parent,txt))
                contentl.append(Newline(parent))
            del contentl[-1] #extra newline
        elif isinstance(t,TeXNode):
            contentl.append(t)
        else:
            raise TypeError('invalid item %s returned from parsing text to nodes'%t)
    
    return contentl
        
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