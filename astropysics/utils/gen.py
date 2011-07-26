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
The :mod:`gen` module contains classes and functions of general utility used in
various places throughout `astropysics`. These are python-specific utilities and
hacks - general data-processing or numerical operations are in the
:mod:`~astropysics.utils.alg` module.

.. todo:: examples?


Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.utils.gen
   :parts: 1

Module API
^^^^^^^^^^

"""

#TODO: implement generalized grid-inverter, possibly tied to modcore.modelgrid1d

from __future__ import division,with_statement
import numpy as np
import re as _re

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
    from collections import MutableMapping
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    MutableMapping = type

def check_type(types,val,acceptnone=True):
    """
    Call this function to check if the value matches the provided types.
    
    :param types:  
        A single type, a sequence of types, or a callable that accepts one
        argument and will be called with the value - if it returns True, the
        value is the correct type.
    :param val: The value to type-check.
    :param acceptnone: Always accept the value None regardless of types
    :type acceptnone: bool
    
    :except TypeError: if the type-checking fails
    :except ValueError: if the `types` are invalid
    
    .. note::
        If any of the types are a :class:`numpy.dtype`, and the value is an 
        :class:`numpy.array`, the data type will also be checked.
    
    """
    if val is None:
        if acceptnone:
            return
        else:
            raise TypeError('None is not a valid value')
        
    if types is not None:
        from operator import isSequenceType
        if isSequenceType(types):
            err = 'Type checking problem'
            for ty in types:
                if isinstance(ty,np.dtype):
                    if not isinstance(val,np.ndarray):
                        err = 'Value %s not a numpy array'%val
                        continue
                    if val.dtype != ty:
                        err = 'Array %s does not match dtype %s'%(val,ty)
                        continue
                elif not isinstance(val,ty):
                    err = 'Value %s is not of type %s'%(val,ty)
                    continue
                return
            raise TypeError(err)
        elif isinstance(types,type) or isinstance(types,np.dtype):
            check_type((types,),val)
        elif callable(types):
            if not types(val):
                raise TypeError('custom function type-checking failed')
        else:
            raise ValueError('invalid type to check')
        



class SymmetricMapping(MutableMapping):
    """
    A dict-like object that maps in two directions - the keys for
    one direction are the value for the other, and vice versa. 
    
    Note that this means all values most be unique and non-mutable.
    
    .. warning::
        This class is probably not at all thread safe.
        
    """
    def __init__(self,*args):
        """
        initialized same as a dict
        """
        self._right = _LinkedDict(None)
        self._left = _LinkedDict(self._right)
        self._right.link = self._left
        if len(args)>0:
            self._left.update(*args) #hackish - better way?
        
    @property
    def forward(self):
        return self._left
    
    @property
    def backward(self):
        return self._right
    
    def __getitem__(self,key):
        return self._left[key]
    def __setitem__(self,key,val):
        self._left[key] = val
    def __delitem__(self,key):
        del self._left[key]
    def __contains__(self,key):
        return key in self._left
    def __len__(self):
        return len(self._left)
    def __iter__(self):
        return iter(self._left)
    def __str__(self):
        return str(self._left)+'/Symmetric'
        
class _LinkedDict(dict):
    """
    A dictionary that applies the reverse of it's operation to a linked
    mapping
    
    .. warning::
        This class is probably not at all thread safe.

    """
    def __init__(self,link,*args,**kwargs):
        self.link = link
        self.setting = False
        self.deling = False
        super(_LinkedDict,self).__init__(*args,**kwargs)
        
    def __getitem__(self,key):
        return super(_LinkedDict,self).__getitem__(key)
    def __setitem__(self,key,val):
        if not self.setting:
            oldval = self[key] if key in self else None
            if isinstance(self.link,_LinkedDict) and self.link.setting and (oldval is not None):
                        raise KeyError('Tried to set a symmetric dict with value %s already present'%key)
            
            super(_LinkedDict,self).__setitem__(key,val)
            if self.link is not None:
                try:
                    self.setting = True
                    self.link[val] = key
                except:
                    if oldval is None:
                        super(_LinkedDict,self).__delitem__(key)
                    else:
                        super(_LinkedDict,self).__setitem__(key,oldval)
                    raise
                finally:
                    self.setting = False
    def __delitem__(self,key):
        if not self.deling:
            val = self[key]
            super(_LinkedDict,self).__delitem__(key)
            if self.link is not None:
                try:
                    self.deling = True
                    del self.link[val]
                finally:
                    self.deling = False
    def update(self,val):
        for k,v in dict(val).iteritems():
            self[k] = v
        
    def pop(self,*args):
        if len(args) > 2:
            raise TypeError('pop expected at most 2 arguments, got '+str(len(args)))
        elif len(args) == 2:
            if key in self:
                val = self[args[0]]
                del self[args[0]]
                return val
            else:
                return args[1]
        elif len(args) == 1:
            if key in self:
                val = self[args[0]]
                del self[args[0]]
                return val
            else:
                raise KeyError(args[0])
        else:
            raise TypeError('pop expected at most 2 arguments, got 0')
        if key in self:
            val = self[key]
            
        else:
            return default
        
        
class DataObjectRegistry(dict):
    """
    A class to register data sets used throughout a module and enable easy 
    access using string names.
    """
    def __init__(self,dataname='data',datatype=None):
        dict.__init__(self)
        self._groupdict = {}
        self.dataname = dataname
        self.datatype = datatype
        
    def __getitem__(self,val):
        if val in self._groupdict:
            data = self._groupdict[val]
            return dict([(d,self[d]) for d in data])
        else:
            
            return dict.__getitem__(self,val)
    def __delitem__(self,key):
        dict.__delitem__(key)
        for k,v in self._groupdict.items():
            if key in v:
                del self.groupdict[k]
                
    def __getattr__(self,name):
        if name in self.keys():
            return self[name]
        else:
            raise AttributeError('No %s or attribute %s in %s'%(self.dataname,name,self.__class__.__name__))
    
    def register(self,objects,groupname=None):
        """
        Register a set of objects, possibly in a group
        
        :param objects: map of the datum name to the associated objects to store
        :type objects: 
            dictionary mapping strings to objects with type matching :attr:`datatype` 
        :param groupname: A common name that applies to all of the objects
        :type groupname: string or None
        
        """
        from operator import isMappingType,isSequenceType
        
        if not isMappingType(objects):
            raise ValueError('input must be a map of bands')
        for k,v in objects.iteritems():
            if self.datatype is not None and not isinstance(v,self.datatype):
                raise ValueError('an object in the %s set is not a %s'%(self.dataname,self.datatype.__name__))
            self[str(k)]=v
            
        if groupname:
            if type(groupname) is str:
                self._groupdict[groupname] = objects.keys()
            elif isSequenceType(groupname):
                for s in groupname:
                    self._groupdict[s] = objects.keys()
            else:
                raise ValueError('unrecognized group name type')
    
    @property        
    def groupnames(self):
        return self._groupdict.keys()
    
    def addToGroup(self,key,group):
        """
        Adds the registered object with the provided name to a group with the 
        given name.
        
        :param key: The object to register
        :type key: string
        :param group: The group name to apply
        :type group: string
        """
        obj = self[key] 
        group = self._groupdict.setdefault(group,[])
        if key not in group:
            group.append(key)
            
    def removeGroup(group):
        """
        Removes the provided group
        
        :param group: name of group to remove
        :type group: string
        """
        del self._groupdict[group]
        
    def getGroupData(self,groupname):
        return self._groupdict[groupname][:]
    
    def getObjects(self,objectstrs,addmissing=False):
        """
        Translates a string, sequence of strings, or mixed sequence of data 
        objects and strings into a sequence of data objects.
        
        :param objectstrs: 
            The requested object names, mixed with objects if desired.
        :type: 
            String, sequence of strings, or mixed sequence of strings and data objects.
        :param addmissing: 
            If True, data objects in `objectstrs` will be added to the registry.
            If the object has a `name` attribute, it will be used as the name,
            otherwise a name will be generated of the form "dataname##".
        :type addmissing: bool
        
        :returns: A sequence of data objects.
        
        """
        from operator import isMappingType,isSequenceType
        
        if self.datatype is not None and isinstance(objectstrs,self.datatype):
            #forces a lone object to be treated as a sequence by itself
            objectstrs = [objectstrs] 
        
        if isinstance(objectstrs,basestring):
            if objectstrs.lower() == 'all':
                objs = self.values()
            elif ',' in objectstrs:
                objs = [self[s.strip()] for s in objectstrs.split(',')]
            elif objectstrs in self:
                obj = self[objectstrs]
                if isMappingType(obj):
                    objs = obj.values()
                else:
                    objs = (obj,)
            elif objectstrs in self.groupnames:
                objs = [self[s] for s in self.getGroupData(objectstrs)]
            else: #assume each character is a data object name
                objs = [self[s] for s in objectstrs.strip()]
        elif isSequenceType(objectstrs):
            notstrs = [not isinstance(s,basestring)for s in objectstrs]
            objs = [self[s] if isinstance(s,basestring) else s for s in objectstrs]
            if self.datatype is not None:
                for i,o in enumerate(objs):
                    if notstrs[i]:
                        if  not isinstance(o,self.datatype):
                            raise ValueError('input object %s is not valid %s'%(o,self.dataname))
                        if addmissing and o not in self.itervalues():
                            if hasattr(o,'name'):
                                if o.name in self:
                                    raise KeyError('%s with name %s already present in this registry'%(self.dataname,o.name))
                                name = o.name
                            else:
                                j = 1
                                while self.dataname+str(j) not in self:
                                    j+=1
                                name = self.dataname+str(j)
                            self[name] = o
                                
        else:
            raise KeyError('unrecognized band(s)')
        
        return objs
    
    def getObjectNames(self,retdict=False):
        """
        Generates and returns  list of names of the objects in this registry, 
        extracted from the `name` attribute of the object if present, otherwise 
        from the registry key.  The order is the same as that generated by
        keys() or values()
        
        :param retdicy: 
            If True, the return value is a dictionary mapping keys from the
            registry to object names.
            
        :returns: List of names.
        """
        ns = []
        ks = []
        for k,o in self.iteritems():
            if hasattr(o,'name'):
                ns.append(o.name)
                ks.append(k)
            else:
                ns.append(k)
                ks.append(k)
        if retdict:
            return dict(zip(ks,ns))
        else:
            return ns
        

#<-----------------------Tools for manipulating docstrings--------------------->
# A regular expression used to determine the amount of space to
# remove.  It looks for the first sequence of spaces immediately
# following the first newline, or at the beginning of the string.
_find_dedent_regex = _re.compile("(?:(?:\n\r?)|^)( *)\S")
# A cache to hold the regexs that actually remove the indent.
_dedent_regex = {}
def change_indentation(s,ind=0):
    """
    Sets the starting indentation of the provided docstring to the given 
    indentation level.
    
    :param s: The string to dedent
    :type s: string
    :param ind: Amount of indentation, either as a string (that will be
                substituted) or an integer number of spaces.
    :type ind: str or int
    
    :returns: a string like `s` with leading indentation removed.
    """
    # This function is based on the matplotlib.cbook.dedent function
    
    if not isinstance(ind,basestring):
        ind = ' '*ind
    indsub = '\n'+ind
    
    if not s:      # includes case of s is None
        return ''

    match = _find_dedent_regex.match(s)
    if match is None:
        return s

    # This is the number of spaces to remove from the left-hand side.
    nshift = match.end(1) - match.start(1)
    if nshift == 0:
        return s

    # Get a regex that will remove *up to* nshift spaces from the
    # beginning of each line.  If it isn't in the cache, generate it.
    unindent = _dedent_regex.get(nshift, None)
    if unindent is None:
        unindent = _re.compile("\n\r? {0,%d}" % nshift)
        _dedent_regex[nshift] = unindent

    result = unindent.sub(indsub, s).strip()
    return result
        
        
def _add_docs_deco(nmdocs,func):
    """
    Decorator partialized for used in `add_docs` and similar.
    
    Note that if astropysics._ignore_add_docs is True, this actually just 
    removes the section in question.
    """
    from .. import _ignore_add_docs
    
    match = _find_dedent_regex.match(func.__doc__)
    if match is None:
        indent = 0
    else:
        indent = match.end(1) - match.start(1)
    doc = func.__doc__
    if doc is None:
        doc = ''
    for n,d in nmdocs:
        replstr = '{docstr:%s}'%n
        d = change_indentation(d,indent) if not _ignore_add_docs else ''
        if replstr in doc:
            doc = doc.replace(replstr,d)
        else:
            doc += ('\n' + (indent*' ') + d)
    func.__doc__ = doc    
    return func
        
def add_docs(*args):
    """
    This class is a decorator indicating that the decorated object should 
    have part (or all) of its docstring replaced by another object's docstring.
    
    The arguments are objects with a `__name__` and `__doc__` attribute. If a
    string of the form '{docstr:Name}' is present in the decorated object's
    docstring, it will be replaced by the docstring from the passed in object
    with `__name__` of 'Name'. Any arguments without a '{docstr:Whatever}' to
    replace will be appended to the end of the decorated object's docstring.
    
    **Examples**   
    
    >>> def f1(x):
    ...     '''Docstring 1'''
    ...     pass
    >>> def f2(x):
    ...     '''
    ...     Docstring 2
    ...     and more!
    ...     '''
    ...     pass
    >>> @add_docs(f1)
    ... def f3(x):
    ...     '''
    ...     Docstring 3
    ...     '''
    ...     pass
    >>> @add_docs(f2)
    ... def f4(x):
    ...     '''
    ...     Docstring 3
    ...     '''
    ...     pass
    >>> @add_docs(f1)
    ... def f5(x):
    ...     '''
    ...     Docstrong 2 {docstr:f1}
    ...     '''
    ...     pass
    >>> f1.__doc__
    'Docstring 1'
    >>> f2.__doc__
    '\\n    Docstring 2\\n    and more!\\n    '
    >>> f3.__doc__
    '\\n    Docstring 3\\n    \\n    Docstring 1'
    >>> f4.__doc__
    '\\n    Docstring 3\\n    \\n    Docstring 2\\n    and more!'
    >>> f5.__doc__
    '\\n    Docstrong 2 Docstring 1\\n    '

    """
    from functools import partial
    
    nmdocs = [(obj.__name__,obj.__doc__) for obj in args]
    return partial(_add_docs_deco,nmdocs)

def add_docs_and_sig(*args):
    """
    Does the same thing as :func:`replace_docs`, but also adds the function
    signature of the argument function to the replaced (followed by a newline).
    Note that this requires that the argument object be a function (and not
    anything with a `__name__` and `__doc__` attribute). This is typically
    useful for functions that do ``f(*args,**kwargs)`` to wrap some other
    function.
    """
    from functools import partial
    from inspect import getargspec
    
    nmdocs = []
    for obj in args:
        #get the wrapped function for classmethod and staticmethod
        if hasattr(obj,'__get__'):
            obj = obj.__get__(0)
        args, varargs, varkw, defaults = getargspec(obj)
        for i,d in enumerate(reversed(defaults)):
            if isinstance(d,basestring):
                ds = "'"+d+"'"
            else:
                ds = str(d)
            args[-1-i] = args[-1-i]+'='+ds
        if varargs:
            args.insert(0,'*'+varargs)
        if varkw:
            args.append('**'+varkw)
        sigstr = obj.__name__+'('+','.join(args)+')'
        newdoc = sigstr+'\n'+obj.__doc__
        nmdocs.append((obj.__name__,newdoc))
    return partial(_add_docs_deco,nmdocs)





del ABCMeta,abstractmethod,abstractproperty #clean up namespace
