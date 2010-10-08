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
        
#<----------Useful regular expressions, mostly used in text_to_nodes----------->
#this finds anything that begins with \begin{<name>} and ends with \end{<name}
#group 1 is the whole env, 2 is the name, 3 is the content
_envre = _re.compile(r'(\\begin{(.*?)}(.*?)\\end{\2})|(\${1,2}.*?\${1,2})',_re.DOTALL)

#This finds anthing like \command{...}{...}... w/ group 1  as the command name,
# 2 the arguments
_commandstr = r'\\(\w.*?)((?:(?:(?:{.*?})|(?:\[.*?\]))+)|(?=\W))'
_commandre = _re.compile(_commandstr)

#this matches either '}' or ']' if it is followed by { or [ or the end of the string
_cmdargsepre = _re.compile(r'(?:}|])(?=\s*?(?:$|{|\[))')


class TeXNode(object):
    """
    An element in the TeX parsing tree.  The main shared characteristic is
    that calling the node will return a string with the combined text.
    
    *Subclassing*
    
    Subclasses must implement :meth:`getSelfText` (see docstring for details)
    
    """
    
    __metaclass__ = ABCMeta
    
    #:The parent of this node in the node tree, or None if this is a root.
    parent = None
    #:A list of child nodes of this node.
    children = tuple()
    
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
    
    #: The :class:`Preamble` object for this file
    preamble = None
    #: The first :class:`Document` environment in this file or None if there isn't one
    document = None
    
    
    def __init__(self,fn=None):
        self.parent = None
        self.children = []
        if fn is not None:
            with open(fn) as f:
                s = f.read()
            self._parse(s)
    
    def _parse(self,f):
        """
        Parse the file - `f` can be a file object or a string w/newlines.
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
        contentstr = '\n'.join(content)
        
        self.preamble = Preamble(self,preamblestr)
        self.children = text_to_nodes(self,contentstr)
        self.children.insert(0,self.preamble)
        
        for c in self.children:
            if isinstance(c,Document):
                self.document = c
                break
        
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
    
    #: The text in this object
    text = ''
    
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
        

class Environment(TeXNode):
    """
    A LaTex environment.  
    
    *Subclassing*
    Subclasses can implement the :meth:`postParse` method - see the method for
    syntax. They should also be registered with the :meth:`registerEnvironment`
    static method to have them be parsed with the default :class:`TeXFile`
    parser. Generally, they should also have a class attribute named `name` that
    gives the name of the environment (this name will be automatically used to
    determine which environments the subclass represents)
    """
    
    #: The name of this environment.
    name = ''
    
    def __init__(self,parent,content,envname=None):
        """
        If `envname` is None, it will be taken from the class-level :attr:`name`
        """
        self.parent = parent
        self.children = c = []
        
        if envname is not None:
            self.name = envname
        elif not hasattr(self,'name'):
            raise ValueError('Environment must have a name')
        
        self.children = self.postParse(text_to_nodes(self,content))
        
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
   
class MathMode(TeXNode):
    """
    A math environment surrounded by $ symbols or $$ (for display mode)
    """
    name = ''
    
    #: determines if the MathMode is in display mode ($$) or not ($)
    displaymode = False
    
    def __init__(self,parent,content):
        """
        `content` should be the full string (including $) or (displaymode,content)
        """
        if isinstance(content,basestring):
            if content.startswith('$$'):
                displaymode = True
                content = content[2:-2]
            elif content.startswith('$'):
                displaymode = False
                content = content[1:-1]
            else:
                raise ValueError('Tried to make MathMode object with no $')
        else:
            displaymode,content = content
            
        self.parent = parent
        self.displaymode = displaymode
        self.children = text_to_nodes(self,content)
    
    def getSelfText(self):
        msymb = '$$' if self.displaymode else '$'
        return (msymb,msymb)

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
        self.children = argnodes = []
        
        if args:
            argspl = _cmdargsepre.split(args)
            if argspl[-1] != '':
                raise ValueError('invalid argument string in Command:"'+args+'"')
            argspl = argspl[:-1]
            for a in argspl:
                if a[0]=='{':
                    argnodes.append(RequiredArgument(self,a[1:]))
                elif a[0]=='[':
                    argnodes.append(OptionalArgument(self,a[1:]))
                else:
                    raise ValueError('invalid set of arguments '+args)
            
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
    
class TrailingCharCommand(Command):
    r"""
    A special command that allows a single trailing character of any type - 
    used for '\left{' '\right]' and similar.
    """
    
    children = tuple() #technically a leaf, but should be an iterable for compatibility w/:class:`Command`
    
    def __init__(self,parent,content):
        """
        `content` can be a (name,char) tuple, or a command string
        """
        if isinstance(content,basestring):
            if not content.startswith('\\'):
                raise ValueError("Tried to make a Command that doesn't start with a backslash")
            name = content[1:-1]
            tchar = content[-1]
        else:
            name,tchar = content
            
        self.parent = parent
        self.name = name
        self.trailingchar = tchar
        
    def getSelfText(self):
        return '\\'+self.name+self.trailingchar
        
class RequiredArgument(TeXNode):
    """
    An argument to a macro that is required (i.e. enclosed in curly braces)
    """
    children = None #arguments are always a leaf
    #: The text of this argument object
    text = ''
    
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
    
class EnclosedDeclaration(TeXNode):
    r"""
    A TeX construct of the form {\name{op}[op] content}. Note that declarations
    terminated by the \end command will not be treated as this kind of object.
    """
    
    #: Whitespace string between the opening brace and the command
    padstr = ''
    #: A :class:`Command` object with the command in the declaration
    cmd = None
    
    def __init__(self,parent,content):
        """
        Content can either be a (padstr,commandnode,innercontent) tuple where
        `padstr` is the string before the command, `commandnode` is a
        :class:`Command` object with the command portion of the declaration, and
        `innercontent` is the content after the command either as a string or a
        list of nodes (possibly mixed with strings). Alternatively, it can be
        the full string including the outermost braces.
        """
        if isinstance(content,basestring):
            if not content.startswith('{\\'):
                raise ValueError("Tried to make enclosed declaration that does't start with {\\")
            if not content.endswith('}'):
                raise ValueError("Tried to make enclosed declaration that does't end with }")
            
            content = content[1:-1] #remove braces
            cmdspl = _commandre.search(content)
            if cmdspl is None:
                raise ValueError('Tried to make an enclosed declaration with no enclosed command')
            padstr = content[:cmdspl.start()]
            
            innercontent = content[cmdspl.end():]
            cmd = Command(self,cmdspl.group(1),'')

        else:
            if len(content) != 3:
                raise ValueError('EnclosedDeclaration content tuple is not length-3')
            padstr,cmd,innercontent = content
            cmd.parent = self
            
        self.parent = parent
        self.padstr = padstr
        self.cmd = cmd
        #content after the command
        if isinstance(innercontent,basestring):
            self.children = text_to_nodes(self,innercontent)
        else:
            self.children = children = []
            for e in innercontent:
                if isinstance(e,basestring):
                    children.extend(text_to_nodes(self,e))
                else:
                    e.parent = self
                    children.append(e)
        
    def getSelfText(self):
        return ('{'+self.padstr+self.cmd(),'}')
    
class Preamble(TeXNode):
    """
    The preamble of a TeX File (i.e. everything before \begin{document} )
    """
    def __init__(self,parent,content):
        self.parent = parent
        self.children = text_to_nodes(self,content)
        
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



def text_to_nodes(parent,txt):
    """
    Converts a string into a list of corresponding TeXNodes.
    """
    txtnodel = []
    #now split out nested environments and commands and build them 
    splitenv = _envre.split(txt) 
    envstrs = splitenv[1::5]
    for i,txts in enumerate(splitenv[::5]):
        #for each text chunk, split out the command  nodes and then add
        #the chunks back to the contentlist
        splitedcms = _commandre.split(txts)
        for j in range(0,len(splitedcms)-1,3):
            txtnodel.append(splitedcms[j])
            txtnodel.append(Command(parent,(splitedcms[j+1],splitedcms[j+2])))
            txtnodel[-1].origtext = splitedcms[j+2]
        #add in last one excluded from above
        if len(splitedcms)>0:
            txtnodel.append(splitedcms[-1])
            
        if i < len(envstrs):
            if envstrs[i] is not None:
                txtnodel.append(environment_factory(parent,envstrs[i]))
            else:
                txtnodel.append(MathMode(parent,splitenv[5*i+4]))
            
    #Fix up any Command nodes that are \left or \right
    for i in reversed(range(len(txtnodel))):
        if isinstance(txtnodel[i],Command):
            if txtnodel[i].name == 'left' and len(txtnodel[i].children)>0:
                #this is for '\left{' and '\left[' because they match to the 
                #command regular expression incorrectly
                cmd = txtnodel[i]
                txt = cmd.origtext
                tchar = txt[0]
                
                tcmd = TrailingCharCommand(parent,(cmd.name,tchar))
                toins = text_to_nodes(parent,txt[1:])
                
                txtnodel[i] = tcmd
                insi = i+1
                for j,e in enumerate(toins):
                    txtnodel.insert(i+j+1,e)
                for c in txtnodel[i].children:
                    print c()
            elif txtnodel[i].name in ('left','right') and \
                    (i<len(txtnodel)-1) and \
                    isinstance(txtnodel[i+1],basestring) and \
                    not txtnodel[i+1][0].isspace():
                     
                nm = txtnodel[i].name 
                tchar = txtnodel[i+1][0]
                txtnodel[i] = TrailingCharCommand(parent,(nm,tchar))
                newtxt = txtnodel[i+1][1:]
                if newtxt == '':
                    del txtnodel[i+1]
                else:
                    txtnodel[i+1] = newtxt
            
            
    #now find Command nodes that are actually supposed to be EnclosedDeclarations
    todel = []
    for i in range(1,len(txtnodel)-1):
        if isinstance(txtnodel[i],Command) and isinstance(txtnodel[i-1],basestring) and isinstance(txtnodel[i+1],basestring):
            #an enclosed dec needs to be "{<maybewhitespace>\command<whitespace>"
            prebracket = txtnodel[i-1].rstrip().endswith('{')
            postwhitespace = len(txtnodel[i+1])>0 and txtnodel[i+1][0].isspace()
            if prebracket and postwhitespace:
                cmdnode = txtnodel[i]
                
                #search for the matched closing brace
                nbraces = 1
                encdeccontent = []
                for j in range(1,len(txtnodel)-i):
                    txtnode = txtnodel[i+j]
                    if isinstance(txtnode,basestring):
                        bracescount = txtnode.count('{') - txtnode.count('}')
                        if (nbraces+bracescount)>1:
                            encdeccontent.append(txtnode)
                            nbraces += bracescount
                            todel.append(i+j)
                        else: #need to find exactly where it finishes
                            for k,char in enumerate(txtnode): #find character that is '}' closer
                                if char=='{':
                                    nbraces+=1
                                elif char=='}':
                                    nbraces-=1
                                if nbraces==0:
                                    break
                            if nbraces!=0:
                                from warnings import warn
                                warn('apparent EnclosedDeclaration has unbalanced braces, skipping')
                                encdeccontent = None
                                break
                                #raise ValueError('apparent EnclosedDeclaration has unbalanced braces')
                            
                            encdeccontent.append(txtnode[:k])
                            #replace the current node with the left-over text
                            txtnodel[i+j] = txtnode[k+1:]
                            break
                    else: #everything else should be balanced, and hence should just get added to the encdec
                        encdeccontent.append(txtnode)
                        
                if encdeccontent is not None:
                    #remove starting '{' and determine the whitespace padding before the command
                    braceind = txtnodel[i-1].rfind('{')
                    padstr = txtnodel[i-1][braceind+1:]
                    txtnodel[i-1] = txtnodel[i-1][:braceind]
                    #add declaration
                    txtnodel[i] = EnclosedDeclaration(parent,(padstr,cmdnode,encdeccontent))
                
    for i in reversed(todel):
        del txtnode[i]
        
    #remove origtxt from any nodes that have it
    for n in txtnodel:
        if hasattr(n,'origtext'):
            del n.origtext 
    
    #now add all the nodes to the child list, 
    #transforming remaining text into TeXt and maybe Newlines if present
    nodel = []
    for t in txtnodel:
        if isinstance(t,basestring):
            for txt in t.split('\n'):
                nodel.append(TeXt(parent,txt))
                nodel.append(Newline(parent))
            del nodel[-1] #extra newline
        elif isinstance(t,TeXNode):
            nodel.append(t)
        else:
            raise TypeError('invalid item %s returned from parsing text to nodes'%t)
    
    return nodel
        
#issue warning if abstract too long
#strip comments
#copy over bbl if necessary
#.tar.gz up with appropriate name and date            
def prep_for_arxiv_pub(texfile):
    raise NotImplementedError

#Tasks: redo figures into f##[l].eps and move the files
#set class to aastex
#strip comments
#fix any deluxetable* -> deluxetable (and add rotate where necessary, also remove [t!] and the like?)
#issue warning if abstract too long
#copy over bib or bbl if necessary
#.tar.gz up with appropriate name and date
#?test compilation?
#remove \comment?
#tablenotemark{$...$} warning?
def prep_for_apj_pub(texfile,authorlast):
    raise NotImplementedError

