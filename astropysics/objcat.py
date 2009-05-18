#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains objects and functions for generating catalogs of objects
where derived quantities are dynamically updated as they are changed.

The basic idea is a tree with the root always a Catalog object 

TODO: modules to also dynamically update via a web server.
*This package is currently under heavy development and subject to major change
without notice
"""

from __future__ import division
from math import pi
import numpy as np

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
    from collections import MutableSequence
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    class MutableSequence(object):
        __slots__=('__weakref__',) #support for weakrefs as in 2.6 MutableSequence objects

class _CatalogElement(object):
    __metaclass__ = ABCMeta
    __slots__=('_fieldnames','_parent','_children')
    
    @abstractmethod
    def __init__(self,parent):
        self._fieldnames = []
        self._parent = parent
        self._children = []
    
    def addField(self,field):
        if not isinstance(field,Field):
            raise ValueError('input value is not a Field')
        setattr(field.name,field)
        self._fieldsnames.append(field.name)
        
    def delField(self,fieldname):
        del self._fields[fieldname]
        delattr(self,fieldname)
        
    @property
    def parent(self):
        return self._parent
        
    #TODO: overwrite __setattr__ and __delattr__ to respond better to Field objects

class Catalog(_CatalogElement):
    """
    This class represents a catalog of objects or catalogs.
    
    A Catalog is essentially a node in the object tree that does not contain 
    fields, only children
    """
    def __init__(self):
        super(Catalog,self).__init__()
        self._fieldnames = None
    
    def addField(self,field):
        raise NotImplementedError('Catalogs cannot have Fields')
    
    def delField(self,fieldname):
        raise NotImplementedError('Catalogs cannot have Fields')
    
    
    

class _SourceMeta(type):
    #TODO: improve Source Singletons
    def __call__(cls,*args,**kwargs):
        obj = type.__call__(cls,*args,**kwargs)
        if not obj._str in Source._singdict:
            Source._singdict[obj._str] = obj
        return Source._singdict[obj._str]

class Source(object):
    __metaclass__ = _SourceMeta
    _singdict = {}
    
    def __init__(self,src):
        self._str = str(src)
        
    def __str__(self):
        return 'Source ' + self._str
    
    
    
    
class Field(MutableSequence):
    """
    This class represents an attribute/characteristic/property of the
    CatalogObject it is associated with.  It stores the current value
    as well as all the other possible values.
    
    note that while the value property will return a FieldValue 
    object (or None), calling the Field directly returns the 
    value attribute of the current FieldValue 
    """
    __slots__=('_name','_type','_vals','_currenti')
    
    def __init__(self,name,type=None):
        """
        The field must have a name, and can optionally be given a type
                
        #TODO:auto-determine name from class
        """
        self._name = name
        self._vals = []
        self._currenti = 0
        self.type = type
        
    def __call__(self):
        v = self.value
        if v is None:
            return None
        else:
            return v.value
    
    def __len__(self):
        return len(self._vals)
    
    def __str__(self):
        return 'Field %s:[%s]'%(self._name,', '.join([str(v) for v in self._vals]))
    
    def _checkValue(self,val,checkdup=True):
        if not (isinstance(val,FieldValue) or (hasattr(val,'source') and hasattr(val,'value'))):
            raise TypeError('Input not FieldValue-compatible')
        if self.type is not None:
            if isinstance(self.type,np.dtype):
                if not isinstance(val.value,np.ndarray):
                    raise TypeError('Value %s not a numpy array'%val)
                if val.value.dtype != self.type:
                    raise TypeError('Array %s does not match dtype %s'%(val,self.type))
            elif not isinstance(val.value,self.type):
                raise TypeError('Value %s is not of type %s'%(val,self.type))
            
        if checkdup:
            for v in self._vals:
                if v is val:
                    raise ValueError('value already present in Field')
        
    def __getitem__(self,key):
        if isinstance(key,Source):
            for v in self._vals:
                #TODO: == -> is performance tests
                if key==v.source:
                    return v
            raise KeyError('Could not find requested Source')
        elif isinstance(key,FieldValue):
            for v in self._vals:
                if key==v:
                    return v
        elif isinstance(key,basestring):
            if 'depends' in key.lower():
                depre = key.lower().replace('depends','').strip()
                if depre=='':
                    depnum=0
                else:
                    depnum=int(depre)
                try:
                    return [v for v in self._vals if v.depends is not None][depnum]
                except IndexError:
                    raise IndexError('dependent value key %i does not exist'%depnum)
            else:
                #TODO:test performance loss
                return self.__getitem__(Source(key))
        else:
            try:
                return self._vals[key]
            except TypeError:
                raise TypeError('Field keys must be strings or list indecies')
    def __setitem__(self,key,val):
        self._checkValue(val)
        i = self._vals.index(self[key])
        self._vals[i] = val
        
    def __delitem__(self,key):
        del self.vals[self._vals.index(self[key])]
    def insert(self,key,val):
        self._checkValue(val)
        if key == len(self._vals):
            i = len(self._vals)
        else:
            i = self._vals.index(self[key])
        self._vals.insert(i,val)
        if self._currenti >= i and len(self._vals) != 1:
            self._currenti += 1
        
    def _getValue(self):
        try:
            return self._vals[self._currenti]
        except IndexError:
            raise IndexError('Field empty')
            #return None
    def _setValue(self,val):
        try:
            self._currenti = self.index(self[val])
        except TypeError,KeyError:
            self._checkValue(val)
            
            self._vals.append(val)
            self._currenti = len(self._vals)-1
    def _delValue(self):
        try:
            del self._vals[self._currenti]
            self._currenti = 0
        except IndexError:
            raise IndexError('deleting from empty Field')
    value = property(_getValue,_setValue,_delValue,
    """
    The current value can be set by setting to a Source object,
    a string matching a Source, or a new FieldValue (adding the value and
    setting it as current.  The current value can also be deleted
    or retrieved using this property
    """)
    def _getType(self):
        return self._type
    def _setType(self,val):
        if val is None:
            self._type = None
        else:
            oldt = self._type
            self._type = val
            try:
                for v in self._vals:
                    self._checkValue(v,checkdup=False)
            except Exception,e:
                self._type = oldt
                raise e
    type = property(_getType,_setType,doc="""
    Selects the type to enforce for this field.  
    if None, no type-checking will be performed
    if a numpy dtype, the value must be an array matching the dtype
    """)
    
    @property
    def name(self):
        return self._name
    
class FieldValue(object):
    __metaclass__ = ABCMeta
    __slots__ = ('_source')
    
    @abstractmethod
    def __init__(self):
        self._source = None
    
    value = abstractproperty()
    
    def _getSource(self):
        return self._source
    def _setSource(self,val):
        if not (val is None or isinstance(val,Source)):
            try:
                val = Source(val)
            except: 
                raise TypeError('Input source is not convertable to a Source object')
        self._source = val 
    source=property(_getSource,_setSource)
    
    def __call__(self):
        return self.value
    
    def __str__(self):
        return str(self.value)
    
    
class ObservedValue(FieldValue):
    """
    This value is a observed or otherwise measured value for the field
    with the associated Source.
    """
    __slots__=('_value')
    def __init__(self,value,source):
        super(ObservedValue,self).__init__()
        if not isinstance(source,Source):
            source = Source(source)
        self.source = source
        self._value = value
        
    def __str__(self):
        return '%s:%s'%(self.value,self.source)
    
    @property    
    def value(self):
        return self._value
        
class DerivedValue(FieldValue):
    """
    This value is derived from a set of other fields (possibly on other 
    objects).  Currently it does not support cycles (where e.g. 
    DerivedValue A depends on DerivedValue B which depends on A)
    """
    __slots__=('_f','_val','dependson')
    def __init__(self,func,dependson=None,source=None):
        """
        The supplied function will be used with the sequence of Field's
        that the values for the derived value are to work with. 
        
        alternatively, if depndson is None, the  depndencies will be 
        inferred from the default values of the 
        """
        super(DerivedValue,self).__init__()
        
        from weakref import ref
        
        if dependson is None:
            from inspect import getargspec
            args,varargs,varkw,defaults = getargspec(func)
            if len(args) != len(defs) or varargs or varkw:
                raise ValueError('input function does not have all defaults \
                                  matched to dependencies or has varargs')
            dependson = defaults 
        #TODO:infer dependencies from argument names and parents/self fields
        
        
        self.dependson = []
        for dep in dependson:
            if not isinstance(dep,Field):
                raise ValueError('provided dependencies are not Fields')
            self.dependson.append(ref(dep))
        
        self._f = func
        self._val = None
        self.source = source
        #TODO:more intelligently work out argument bindings
    
    def __str__(self):
        return '%s:Derived'%self.value
        
    @property
    def value(self):
        if self._val is None:
            self._val = self._f(*(fieldwr()() for fieldwr in self.dependson))
            #TODO:check for bad wrs?
        return self._val
    
    @staticmethod
    def derivedFunc(func):
        """
        Use this as a decorator to convert a function with defaults into
        a DerivedValue instance with the defaults as the dependencies
        """
        from inspect import getargspec
        args,varargs,kwds,defs=getargspec(func)
        if len(args) != len(defs) or varargs or kwds:
            raise ValueError('input function does not have all defaults \
                              matched to dependencies or has varargs')
        return DerivedValue(func,dependson=defs)
        
    
    
    
class _CatalogObjectMeta(ABCMeta):
    def __call__(cls,*args,**kwargs):
        obj = super(_CatalogObjectMeta,cls).__call__(*args,**kwargs)
        print dir(obj)
        return obj

class CatalogObject(_CatalogElement):
    __metaclass__ = _CatalogObjectMeta
    
    def __init__(self):
        super(CatalogObject,self).__init__()
        raise NotImplementedError

    
del ABCMeta,abstractmethod,abstractproperty,MutableSequence,pi,division #clean up namespace
  
