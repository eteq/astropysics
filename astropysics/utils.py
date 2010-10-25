#Copyright 2009 Erik Tollerud
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

======================================
utils -- utility classes and functions
======================================

The :mod:`utils` module contains classes and funtions of general utility used in
multiple places throughout `astropysics`. Some of these are
astropyhysics-specific algorithms while others are more python tricks.

.. todo:: examples/tutorials


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.utils
   :parts: 1

Module API
----------

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
    
#<------------------------Python/generic utilities----------------------------->

def funpickle(fileorname,number=0,usecPickle=True):
    """
    Unpickle a file specified by  and return the pickled contents.
    
    :param fileorname: The file from which to unpickle objects
    :type fileorname: a file name string or a :class:`file` object
    :param number: 
        The number of objects to unpickle - if <1, returns a single object.
    :type number: int
    :param usecPickle: 
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster).
    :type usecPickle: bool 
    
    :returns: A list of length given by `number` or a single object if number<1
    
    """
    if usecPickle:
        import cPickle as pickle
    else:
        import pickle
        
    if isinstance(fileorname,basestring):
        f = open(fileorname,'r')
        close = True
    else:
        f = fileorname
        close = False
        
    try:
        if number > 0:
            res = []
            for i in range(number):
                res.append(pickle.load(f))
        elif number < 0:
            res = []
            eof = False
            while not eof:
                try:
                    res.append(pickle.load(f))
                except EOFError:
                    eof = True
        else: #number==0
            res = pickle.load(f)
    finally:
        if close:
            f.close()
            
    return res

def fpickle(object,fileorname,usecPickle=True,protocol=None,append=False):
    """
    Pickle an object to a specified file.
    
    :param object: the python object to pickle
    :param fileorname: The file from which to unpickle objects
    :type fileorname: a file name string or a :class:`file` object
    :param usecPickle: 
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster).
    :type usecPickle: bool 
    :param protocol: 
        Pickle protocol to use - see the :mod:`pickle` module for details on
        these options. If None, the most recent protocol will be used.
    :type protocol: int or None
    :param append:
        If True, the object is appended to the end of the file, otherwise the
        file will be overwritten (if a file object is given instead of a 
        file name, this has no effect).
    :type append: bool
    
    """
    
    if usecPickle:
        import cPickle as pickle
    else:
        import pickle
        
    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL
    
    if isinstance(fileorname,basestring):
        f = open(fileorname,'a' if append else 'w')
        close = True
    else:
        f = fileorname
        close = False 
        
    try:
        pickle.dump(object,f,protocol=protocol)
    finally:
        if close:
            f.close()

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
        if np.iterable(types):
            err = 'Type checking problem'
            for ty in types:
                if isinstance(ty,np.dtype):
                    if not isinstance(val,np.ndarray):
                        err = 'Value %s not a numpy array'%val
                        continue
                    if self.value.dtype != ty:
                        err = 'Array %s does not match dtype %s'%(val,ty)
                        continue
                elif not isinstance(val,ty):
                    err = 'Value %s is not of type %s'%(val,ty)
                    continue
                return
            raise TypeError(err)
        elif isinstance(types,type):
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
    from . import _ignore_add_docs
    
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
     
    .. testsetup::
    
        from astropysics.utils import add_docs
    
    .. doctest::    
    
        >>> def f1(x):
        >>>     '''Docstring 1'''
        >>>     pass
        >>> def f2(x):
        >>>     '''
        >>>     Docstring 2
        >>>     and more!
        >>>     '''
        >>>     pass
        >>> @add_docs(f1)
        >>> def f3(x):
        >>>     '''
        >>>     Docstring 3
        >>>     '''
        >>>     pass
        >>> @add_docs(f2)
        >>> def f4(x):
        >>>     '''
        >>>     Docstring 3
        >>>     '''
        >>>     pass
        >>> @add_docs(f1)
        >>> def f5(x):
        >>>     '''
        >>>     Docstrong 2 {docstr:f1}
        >>>     '''
        >>>     pass
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

#<--------------------Analysis/simple numerical functions---------------------->
def estimate_background(arr,method='median'):
    """
    Estimates the background of the provided array following a technique
    specified by the `method` keyword:
    
    * 'median' : median of all the data
    * 'mean' : mean of all the data
    * '32' : 3*median - 2*mean
    * '2515' : 2.5*median - 1.5*mean
    * '21' : 2*median - 1*mean
    * a callable that takes a 1D array and outputs a scalar
    * a scalar : always returns that value
    * None : always returns a 0 background
    
    outputs a scalar background estimate
    """
    arr = np.array(arr,copy=False).ravel()
    
    if method is None:
        res = 0
    elif np.isscalar(method) or (hasattr(method,'shape') and method.shape is tuple()):
        res = method
    elif callable(method):
        res = method(arr)
    elif method == 'median':
        res = np.median(arr)
    elif method == 'mean':
        res = np.mean(arr)
    elif method == '32':
        res = (3*np.median(arr) - 2*np.mean(arr))
    elif method == '2515':
        res = (2.5*np.median(arr) - 1.5*np.mean(arr))
    elif method == '21':
        res = (2*np.median(arr) - np.mean(arr))
    else:
        raise ValueError('unrecognized offset type')
    
    return res

def moments(arr,ms,axes=None,bgmethod=None,norm=True,std=False):
    """
    Compute the moments of the provided n-d array.  That is 
    
    :param arr: The input array for which the moments are desired.
    :param ms: 
        The desired order of the moment to be returned.  Can be:
            
            * A scalar
                The same order will be used for each dimension, and the return 
                type will be an array of the moment along each dimension.
            * A sequence 
                The sequence must match the number of dimensions in `arr`, and
                specifies the order along each dimension.  e.g. for a 3D 
                cartesian array, ms=[1,3,2] means the moemnt comupted will be
                :math:`x^1 y^3 z^2`.
    :param axes: 
        If None, the spatial axes are taken to be 0,1,2...,nd-1 for each of the
        dimenstions in the input. Otherwise, `axes` must be a seqence of arrays
        specifying the spatial locations along each axis. This sequence must be
        of length matching the number of dimensions in the input, and Each
        element must be of the same length as the corresponding input array
        dimension.
    :param bgmethod: 
        A background-estimation method (see :func:`estimate_background` for
        options), a scalar to subtract from the array, or None to do no
        subtraction.
    :param bool norm: 
        If True, the moment will be normalized - i.e. ``sum(x^m*arr)/sum(arr)``
        instead of just ``sum(x^m*arr)``
    :param bool std: 
        If True, the output will be standardized (mth moment divided by standard
        deviation to the mth power)
    
    :returns: 
        Either the computed moment if `ms` is a sequence, or a 1D array of
        moments for each dimension if `ms` is a scalar.
    
    """
   
    arr = np.array(arr,copy=False)
    if bgmethod:
        arr = arr - estimate_background(arr,bgmethod)
    shp = arr.shape
    
    #setup/check axes
    if axes is None:
        axes = [np.arange(s) for s in shp]
    elif len(axes) != len(shp):
        raise ValueError('incorrect number of axes provided')
    else:
        axmatch = np.array([len(ax)==s for ax,s in zip(axes,shp)])
        if np.any(~axmatch):
            raise ValueError('axes %s do not match input array'%np.where(~axmatch))
        
    newax = []
    for i,ax in enumerate(axes):
        shi = np.ones(len(shp))
        shi[i] = len(ax)
        newax.append(ax.reshape(shi))
    axes = newax
    
    if np.isscalar(ms):
        res = np.array([np.sum(arr*ax**ms) for ax in axes])
    else:
        if len(ms) != len(shp):
            raise ValueError('moment sequence does not match data')
        #bcast = np.broadcast_arrays(*[ax**m for m,ax in zip(ms,axes)])
        #res = np.sum(arr*np.prod(bcast,axis=0))
        axprod = reduce(np.multiply,[ax**m for m,ax in zip(ms,axes)])
        res = np.sum(arr*axprod)
    
    if std:
        if np.isscalar(ms):
            res = res/np.array([np.std(ax)**ms for ax in axes])
        else:
            res = res/np.prod([np.std(ax)**m for m,ax in zip(ms,axes)])
    
    if norm:
        res/=np.sum(arr)
    return res

def centroid(val,axes=None,offset=None):
    """
    Convinience function calling :func:`moments`  to get just the first
    normalized moment (e.g. the centroid).
    
    :param val: n-d array for which to compute the centroid.
    :type val: array-like
    :param axes: 
        None to use default 0-based location scheme (see :func:`moments`) or an
        array of poisition values for the centriod (must have as many elements
        as the dimension of `val`)
    :type axes: array-like or None
    :param offset: 
        A fixed offset to subtract from the data before centroiding, or a
        `bgmethod` like those accepted by :func:`estimate_background`.
    """
    if axes is not None and np.isscalar(axes[0]):
        axes = (axes,)
    return moments(val,1,axes,offset,True,False)

def sigma_clip(data,sig=3,iters=1,bkgmeth='median',cenfunc=np.var,maout=False):
    """
    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.
    
    :param data: input data (will be flattened to 1D)
    :type data: array-like
    :param sig: 
        The number of standard deviations to use as the clipping limit, or 
        the square root of the variance limit.
    :type sig: scalar
    :param iters: number of iterations to perform clipping
    :type iters: int
    :param cenfunc: 
        The technique to compute the center for the clipping - may be any valid
        input to :func:`estimate_background`
    :param varfunc: 
        The method to compute the variance about the. Should take a 1D array as
        input and output a scalar. This will be compared to the square of the
        data as if it is the variance.
    :type varfunc: a callable
    :param maout: If True, return a masked array (see return value for details).
    :type maout: bool
    
    :returns: 
        A :class:`numpy.ma.Maskedarray` with the rejected points masked, if
        `maout` is True. If maout is False, a tuple (filtereddata,mask) is
        returned where the mask is False for rejected points (and matches the
        shape of the input).
    
    """
    data = np.array(data,copy=False)
    oldshape = data.shape
    data = data.ravel()
    
    mask = np.ones(data.size,bool)
    for i in range(iters):
        do = data-estimate_background(data[mask],cenfunc)
        mask = do*do <= varfunc(data[mask])*sig**2
        
    if maout:
        return np.ma.MaskedArray(data,~mask,copy='maout'=='copy')
    else:
        return data[mask],mask.reshape(oldshape)
    
def nd_grid(*vecs):
    """
    Generates a grid of values given a sequence of 1D arrays. The inputs 
    will be converted to 1-d vectors and the output is an array with
    dimensions (nvecs,n1,n2,n3,...) varying only on the dimension 
    corresponding to its input order
    
    
    **Examples**::
    
        x = linspace(-1,1,10)
        y = x**2+3
        z = randn(13)
        result = nd_grid(x,y,z)
        xg,yg,zg = result
    
    result will be a (3,10,10,13) array, and each of xg,yg, and zg are 
    (10,10,13)
    """
    vecs = [np.array(v,copy=False).ravel() for v in vecs]
    shape = tuple([v.size for v in vecs])
    sz = np.prod(shape)
    
    arrs = []
    for i,(v,n) in enumerate(zip(vecs,shape)):
        newshape = list(shape)
        newshape.insert(0,newshape.pop(i))
        
        arr = np.repeat(v,sz/n).reshape(newshape)
        
        transorder = range(len(arr.shape))[1:]
        transorder.insert(i,0)
        
        arrs.append(arr.transpose(transorder))
    return np.array(arrs)

def lin_to_log_rescale(val,lower=1,upper=3,base=10):
    """
    Linearly rescales input values onto the range [base^lower,base^upper] and
    then applies the requested base of logarithm to the rescaled values.
    
    :param val: input linear values
    :type val: array-like
    :param lower: lower value for output values
    :type lower: scalar
    :param upper: upper value for output values
    :type upper: scalar
    :param base: base of logarithm
    :type base: base of logarithm
    
    :returns: logarithm of rescaled input
    
    """
    """
    maps linear values onto the range [lower,upper] and then applies the
    requested base of logarithm
    """
    if lower > upper:
        raise ValueError('lower must be less than upper')
    
    lower = base**lower
    upper = base**upper
    
    val = np.array(val,copy=False)
    #offset to [0,something]
    val = val - val.min()
    #reacale to [0,range]
    val *= ((upper-lower)/val.max())
    val += lower

    if base is None:
        return np.log(val)
    elif base==10:
        return np.log10(val)
    else:
        return np.log(val)/np.log(base)

def crossmask(x,threshold=0,belowtoabove=True):
    """
    Generates a boolean mask for the point where the input crosses or has
    complted crossing (if between) from below (or above) a threshold.
    
    If `belowtoabove` is True, the returned masks is for where the input
    transitions from below to above.  Otherwise, from above to below.
    
    :param x: input (will be flattened to 1D)
    :type x: array-like
    :param threshold: the transition value for where the crossings occur
    :type threshold: scalar
    :param belowtoabove: 
        If True, the returned mask is for where the input transitions from below
        to above. Otherwise, from above to below.
    :type belowtoabove: bool
    
    :returns: A mask that is True where crossing occurs, False everywhere else.
    :rtype: bool :class:`~numpy.ndarray`
    
    **Example**
    
    .. testsetup::
    
        from astropysics.utils import crossmask
        from numpy import array,where
    
    .. doctest::
        
        >>> xup = [-2,-1,0,1,2]
        >>> xdn = [2,1,0,-1,-2]
        
        >>> print crossmask(xup,0,True)
        [False False  True False False]
        >>> print crossmask(xup,-0.5,True)
        [False False  True False False]
        >>> print crossmask(xup,0.5,True)
        [False False False  True False]
        
        >>> print crossmask(xdn,0,True)
        [False False False False False]
        
        >>> print crossmask(xdn,0,False)
        [False False  True False False]
        >>> print crossmask(xdn,0.5,False)
        [False False  True False False]
        >>> print crossmask(xdn,-0.5,False)
        [False False False  True False]
        
        >>> xupdnup = [-2,-1,0,1,2,1,0,-1,-2,-1,0,1,2]
        >>> where(crossmask(xupdnup,0.5,True))
        (array([ 3, 11]),)
        >>> print array(xupdnup)[crossmask(xupdnup,0,True)]
        [0 0]


        
    """
    x = np.array(x,copy=False).ravel()
    if belowtoabove:
        a = x <  threshold
        b = x >= threshold
    else:
        a = x >  threshold
        b = x <= threshold
        
    mask = np.roll(a,1)&b
    mask[0] = False
    return mask

#<--------------------------Rotations------------------------------------------>
def rotation_matrix(angle,axis='z',degrees=True):
    """
    Generate a 3x3 rotation matrix in cartesian coordinates for rotation about
    the requested axis.
    
    :param axis:
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate about. If
        'x','y', or 'z', the rotation sense is counterclockwise looking down the
        + axis (e.g. positive rotations obey left-hand-rule).
    :type axis: string or 3-sequence
    :param degrees: If True the input angle is degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: A :class:`numpy.matrix` unitary rotation matrix.
    """
    from math import sin,cos,radians,sqrt
    if degrees:
        angle = radians(angle)
        
    
    
    if axis == 'z':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, s, 0),
                          (-s, c, 0),
                          ( 0, 0, 1)))
    elif axis == 'y':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, 0,-s),
                          ( 0, 1, 0),
                          ( s, 0, c)))
    elif axis == 'x':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( 1, 0, 0),
                          ( 0, c, s),
                          ( 0,-s, c)))
    else:
        x,y,z = axis
        w = cos(angle/2)
        
        #normalize
        if w == 1:
            x=y=z=0
        else:
            l = sqrt((x*x + y*y + z*z)/(1-w*w))
            x /= l
            y /= l
            z /= l
        
        wsq = w*w
        xsq = x*x
        ysq = y*y
        zsq = z*z
        return np.matrix((( wsq+xsq-ysq-zsq, 2*x*y-2*w*z, 2*x*z+2*w*y),
                          ( 2*x*y+2*w*z, wsq-xsq+ysq-zsq,2*y*z-2*w*x),
                          ( 2*x*z-2*w*y, 2*y*z+2*w*x, wsq-xsq-ysq+zsq)))
                          
                          
def angle_axis(matrix,degrees=True):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix.
    
    :param matrix: the rotation matrix
    :type matrix: a 3x3 :class:`numpy.ndarray`
    :param degrees: if True, output is in degrees
    :type degrees: boolean
    
    :returns:
        an (angle,axis) tuple where the angle is in degrees if `degrees` is
        True, otherwise in radians
        
    """
    from math import sin,cos,acos,degrees,sqrt
    
    m = np.asmatrix(matrix)
    if m.shape != (3,3):
        raise ValueError('matrix is not 3x3')
    
    
    
    angle = acos((m[0,0] + m[1,1] + m[2,2] - 1)/2)
    denom = sqrt(2*((m[2,1]-m[1,2])+(m[0,2]-m[2,0])+(m[1,0]-m[0,1])))
    axis = np.array((m[2,1]-m[1,2],m[0,2]-m[2,0],m[1,0]-m[0,1]))/denom
    axis /= sqrt(np.sum(axis**2)) 
    
    if degrees:
        return degrees(angle),axis
    else:
        return angle,axis
    
def rotation_matrix_xy(x,y):
    """
    Computes the rotation matrix that moves the +z-axis (pole) to a new vector
    (x,y,z') where x and y are specified and z' is constrained by requiring
    the vector to be length-1.
    
    :param x: The x-value of the new pole in the starting coordinate system
    :type x: float
    :param y: The y-value of the new pole in the starting coordinate system
    :type y: float
    
    :returns: A :class:`numpy.matrix` unitary rotation matrix.
    
    """
    from math import sqrt
    
    xsq = x*x
    ysq = y*y
    xy = x*y
    
    z = sqrt(1-xsq-ysq)
    b = 1/(1+z)
    
    return np.matrix([[1-b*xsq,-b*xy,-x],
                      [-b*xy,1-n*ysq,-y],
                      [x,y,1-b*(xsq+ysq)]])

#<--------------------------Robust statistics---------------------------------->    
def interquartile_range(values,scaletonormal=False):
    """
    Computes the interquartile range for the provided sequence of values, a
    more robust estimator than the variance.
    
    :param values: the values for which to compute the interquartile range
    :type values: array-like, will be treated as 1D
    :param scaletonormal: Rescale so that a normal distribution returns 1
    :type scaletonormal: bool
    
    :returns: the interquartile range as a float
    
    """
    from scipy.stats import scoreatpercentile
    from scipy.special import erfinv
    
    x = np.array(values,copy=False).ravel()
    res = scoreatpercentile(x,75) - scoreatpercentile(x,25)
    if scaletonormal:
        nrm = 8**0.5*erfinv(.5)
        return res/nrm
    else:
        return res
    
def median_absolute_deviation(values,scaletonormal=False,cenestimator=np.median):
    """
    Computes the median_absolute_deviation for the provided sequence of values, 
    a more robust estimator than the variance.
    
    :param values: the values for which to compute the MAD
    :type values: array-like, will be treated as 1D
    :param scaletonormal: Rescale the MAD so that a normal distribution is 1
    :type scaletonormal: bool
    :param cenestimator: 
        A function to estimate the center of the values from a 1D array of
        values. To actually be the "median" absolute deviation, this must be
        left as the default (median).    
    :type cenestimator: callable
    
    :returns: the MAD as a float
    
    """
    from scipy.special import erfinv
    
    x = np.array(values,copy=False).ravel()
    res = np.median(np.abs(x-np.median(x)))
    
    if scaletonormal:
        nrm = (2**0.5*erfinv(.5))
        return res/nrm
    else:
        return res
    
def biweight_midvariance(values,influencescale=9,cenestimator=np.median):
    """
    Computes the biweight midvariance of a sequence of data points, a robust 
    statistic of scale.  
    
    For normal and uniform distributions, it is typically close to, but a bit
    above the variance.
    
    :param values: the values for which to compute the biweight
    :type values: array-like, will be treated as 1D
    :param influencescale: 
        The number of MAD units away at which a data point has no weight
    :type influencescale: int
    :param cenestimator: 
        A function to estimate the center of the values from a 1D array of
        values. To be a true standard biweight midvariance, this must be left as
        the default (median).
    :type cenestimator: callable
    
    :returns: biweight,median tuple (both floats)
    
    """
       
    x = np.array(values,copy=False).ravel()
    
    Q = cenestimator(x)
    MAD = median_absolute_deviation(x,cenestimator=cenestimator)
       
    ui = (x-Q)/(influencescale*MAD)
    uisq=ui**2
    
    I = uisq <= 1
    
    numer = np.sum(I*(x-Q)**2*(1-uisq)**4)
    denom = np.sum(I*(1-uisq)*(1-5*uisq))**2
    
    
    return x.size*numer/denom,Q

del ABCMeta,abstractmethod,abstractproperty #clean up namespace
