#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains objects and functions for generating catalogs of objects
where derived quantities are dynamically updated as they are changed.

The basic idea is a tree/DAG with the root typically a Catalog object

TODO: modules to also dynamically update via a web server.
"""

from __future__ import division
from math import pi
import numpy as np

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
    from collections import Sequence,MutableSequence,MutableMapping
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    class MutableSequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
    class MutableMapping(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
    class Sequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
        
class CycleError(Exception):
    """
    This exception indicates a cycle was detected in some graph-like structure
    """
    def __init__(self,message):
        super(CycleError,self).__init__(message)

class CatalogNode(object):
    """
    This Object is the superclass for all elements/nodes of a catalog.  
    This is an abstract class that must have its initializer overriden.
    
    Subclasses must call super(Subclass,self).__init__(parent) in their __init__
    """
    
    __metaclass__ = ABCMeta
    __slots__=('_parent','_children','__weakref__')
    
    @abstractmethod
    def __init__(self,parent):
        self._children = []
        self._parent = None
        
        if parent is not None:
            self.parent = parent
        
    def _cycleCheck(self,source):
        """
        call this from a child object with the child as the Source to check 
        for cycles in the graph
        """
        if source is self:
            raise CycleError('cycle detected in graph assignment attempt')
        if self._parent is None:
            return None
        else:
            return self.parent._cycleCheck(source)
            
    def _getParent(self):
        return self._parent
    def _setParent(self,val):
        
        if val is not None:
            val._cycleCheck(self) #TODO:test performance effect/make disablable
            val._children.append(self)
            
        if self._parent is not None:
            self._parent._children.remove(self)
        self._parent = val
        
    parent=property(_getParent,_setParent)
    
    
    @property
    def children(self):
        return tuple(self._children)
    
    def reorderChildren(self,neworder):
        """
        Change the order pf the children
        
        neworder can be either a sequence of  indecies (e.g. to reorder
        [a,b,c] to [c,a,b], neworder would be [2,0,1]), the string
        'reverse', or a function like the cmp keyword as would appear
        in the sorted builtin (can be None to do default sorting). 
        """
        if neworder == 'reverse':
            self._children.reverse()
        elif callable(neworder):
            self._children.sort(cmp=neworder)
        else: #TODO:faster way to do this if necessary?
            if len(neworder) != len(self._children):
                raise ValueError('input sequence does not have correct number of elements')
            
            added = np.zeros(len(self._children),dtype=bool)
            newl = []
            for i in neworder:
                if added[i]:
                    raise ValueError('input sequence has repeats')
                newl.append(self._children[i])
                added[i] = True
                
    @property
    def nnodes(self):
        """
        this gives the number of total nodes at this point in the tree
        (including self - e.g. a leaf in the tree returns 1)
        """
        return sum([c.nnodes for c in self._children],1)
        
    
    def visit(self,func,traversal='postorder',filter=False):
        """
        This function walks through the object and all its children, 
        executing func(CatalogNode)
        
        traversal is the traversal order of the tree - can be:
        *'preorder',
        *'postorder'
        *an integer indicating at which index the root should 
        be evaluated (pre/post are 0/-1)
        * a float between -1 and 1 indicating where the root should 
        be evaluated as a fraction
        *'level'/'breathfirst' 
        
        filter can be:
        *False: process and return all values
        *a callable: is called as g(node) and if it returns False, the
        node will not be processed nor put in the (also ignores anything 
        that returns None)
        *any other: if the node returns this value on processing, it will
                    not be returned
        """
        if callable(filter):
            func = lambda *args,**kwargs:func(*args,**kwargs) if filter(args[0]) else None
            filter = None
        
        if type(traversal) is int:
            retvals = []
            doroot = True
            for i,c in enumerate(self._children):
                if i == traversal:
                    retvals.append(func(self))
                    doroot = False
                retvals.extend(c.visit(func,traversal))
            if doroot:
                retvals.append(func(self))
        elif type(traversal) is float:
            retvals = []
            doroot = True
            travi = int(traversal*self._children)
            for i,c in enumerate(self._children):
                if i == travi:
                    retvals.append(func(self))
                    doroot = False
                retvals.extend(c.visit(func,traversal))
            if doroot:
                retvals.append(func(self))
        elif traversal is None: #None means postorder
            retvals = []
            for c in self._children:
                retvals.extend(c.visit(func,traversal))
            retvals.append(func(self))    
        elif traversal == 'postorder':
            retvals = self.visit(func,None)
        elif traversal == 'preorder':
            retvals = self.visit(func,0)
        elif traversal == 'level' or traversal == 'breadthfirst':
            from collections import deque
            
            retvals=[]
            q = deque()
            q.append(self)
            while len(q)>0:
                elem = q.popleft()
                retvals.append(func(elem))
                q.extend(elem._children)
        else:
            raise ValueError('unrecognized traversal type')
        
        if filter is not False:
            retvals = [v for v in retvals if v is not filter]
        return retvals
    
class FieldNode(CatalogNode,MutableMapping,Sequence):
    """
    A node in the catalog that has Fields.  This is an abstract class that 
    must have its initializer overriden.
    
    Note that for these subclasses, attribute access (e.g. node.fieldname) 
    accesses the Field object, while mapping or sequence-style access 
    (e.g node['fieldname'] or node[1])  directly accesses the current value
    of the field (or None if there is no value)
    """
    __slots__=('_fieldnames',)
    
    @abstractmethod
    def __init__(self,parent):
        super(FieldNode,self).__init__(parent)
        self._fieldnames = []
        
    def addField(self,field):
        if not isinstance(field,Field):
            raise ValueError('input value is not a Field')
        if field.name in self._fieldnames:
            raise ValueError('Field name "%s" already present'%field.name)
        setattr(self,field.name,field)
        self._fieldnames.append(field.name)
        
    def delField(self,fieldname):
        try:
            self._fieldnames.remove(fieldname)
            if hasattr(self.__class__,fieldname):
                setattr(self,fieldname,None)
            else:
                delattr(self,fieldname)
        except ValueError:
            raise KeyError('Field "%s" not found'%fieldname)
        
    #def __iter__(self):
    #    return iter(self._fieldnames)
    
    def __len__(self):
        return len(self._fieldnames)
    
    def __contains__(self,key):
        return key in self._fieldnames
        
    def __getitem__(self,key):
        if key not in self._fieldnames:
            try:
                key = self._fieldnames[key]
            except (IndexError,TypeError):
                raise KeyError('Field "%s" not found'%key)
        try:
            return getattr(self,key)()
        except IndexError: #field empty
            return None
    
    def __setitem__(self,key,val):
        if key not in self._fieldnames:
            try:
                key = self._fieldnames[key]
            except (IndexError,TypeError):
                raise KeyError('Field "%s" not found'%key)
        field = getattr(self,key)
        field.currentobj = val
    
    def __delitem__(self,key):
        self.delField(key)
        
    @property
    def fieldnames(self):
        return tuple(self._fieldnames)
    
    def extractField(self,*args,**kwargs):
        """
        walk through the tree starting from this object
        
        see FieldNode.extractFieldFromNode for arguments
        """
        return FieldNode.extractFieldAtNode (self,*args,**kwargs)
    
    @staticmethod
    def extractFieldAtNode(node,fieldname,traversal='postorder',missing=False,dtype=None):
        """
        this will walk through the tree starting from the Node in the first
        argument and generate an array of the values for the 
        specified fieldname
        
        missing determines the behavior in the event that a field is not 
        present (or a non FieldNode is encounterd) it can be:
        *'exception': raise a KeyError if the field is missing or a 
        TypeError if  
        *'skip': do not include this object in the final array
        *'0'/False: 
        
        traversal is of an argument like that for CatalogNode.visit
        """
        #TODO: optimize with array size knowledge ?
        if missing == 'exception':
            filter = False
            def vfunc(node):
                return node[fieldname]
        elif missing == 'skip':
            filter = None
            def vfunc(node):
                try:
                    return node[fieldname]
                except (KeyError,TypeError):
                    return None
        elif not missing:
            filter = False
            def vfunc(node):
                try:
                    return node[fieldname]
                except (KeyError,TypeError):
                    return None
        else:
            raise ValueError('Unrecognized value for what to do with missing fields')
            
        lst = node.visit(vfunc,traversal=traversal,filter=filter)
        
        if dtype is None:
            try:
                #TODO: test this or be smarter
                return np.array(lst,dtype=node.type)
            except:
                pass
        
        return np.array(lst,dtype=dtype)
        
        
        
        
    
    #TODO: overwrite __setattr__ and __delattr__ to respond better to Field objects
   
#<----------------------------node attribute types----------------------------->    

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
    
#TODO: more Source information
    
class Field(MutableSequence,MutableMapping):
    """
    This class represents an attribute/characteristic/property of the
    CatalogObject it is associated with.  It stores the current value
    as well as all the other possible values.

    The values, sources, and default properties will return the actual values 
    contained in the FieldValues, while currentobj and iterating 
    over the Field will return FieldValue objects.  Calling the 
    Field (no arguments) will return the current value
    
    usedef specified if the default should be set -- if True, defaultval will 
    be used, if None, a None defaultval will be ignored but any other
    will be recognizd, and if False, no default will be set
    """
    __slots__=('_name','_type','_vals','_notify')
    
    def __init__(self,name,type=None,defaultval=None,usedef=None):
        """
        The field must have a name, and can optionally be given a type
                
        #TODO:auto-determine name from class
        """        
        self._name = name
        self._vals = []
        self._type = None
        self._notifywrs = None
        
        self.type = type
        if usedef or (usedef is None and defaultval is not None):
            self.default = defaultval
        
    def __call__(self):
        return self.currentobj.value
    def __len__(self):
        return len(self._vals)
    def __contains__(self,val):
        #TODO: decide if this performance hit is worth optimizing somehow
        if val is None:
            try:
                self[None]
                return True
            except KeyError:
                return False
        else:
            return super(Field,self).__contains__(val)
    
    def __str__(self):
        return 'Field %s:[%s]'%(self._name,', '.join([str(v) for v in self._vals]))
    
    def _checkConvInVal(self,val,dosrccheck=True):
        """
        auto-converts tuples to ObservedValues
        #TODO: auto-convert callables with necessary information to derivedvalues
        
        dosrccheck = True -> check if source is present
        dosrccheck = string/Source -> ensure that value matches specified 
        string/Source
        dosrcchecj = False -> do nothing but convert
        """
        
        if isinstance(val,tuple) and self.type is not tuple:
            val = ObservedValue(*val)   
        
        if not (isinstance(val,FieldValue) or (hasattr(val,'source') and hasattr(val,'value'))):
            raise TypeError('Input %s not FieldValue-compatible'%str(val))
        
        if dosrccheck:
            if isinstance(dosrccheck,Source):
                s = dosrccheck
                if val.source != s:
                    raise ValueError('Input %s does not match expected %s' %(val.source,s))
            elif isinstance(dosrccheck,basestring):
                s = Source(dosrccheck)
                if val.source != s:
                    raise ValueError('Input %s does not match expected %s' %(val.source,s))
            else:
                s = None
                for v in self._vals:
                    if v.source is val.source:
                        raise ValueError('value with %s already present in Field'%v.source)
            
        val.checkType(self.type)
        return val
    
    def notifyValueChange(self,oldval,newval):
        """
        notifies all registered functions that the value in this 
        field has changed
        
        (see registerNotifier)
        """
        #TODO: optimize better
        if self._notifywrs is not None:
            deadrefs=[]
            for i,wr in enumerate(self._notifywrs):
                callobj = wr()
                if obj is None:
                    deadrefs.append(i)
                else:
                    callobj(oldval,newval)
            
            if len(deadrefs) == len(self._notifywrs):
                self._notifywrs = None
            else:
                for i in reversed(deadrefs):
                    del self._notifywrs[i]
    
    def registerNotifier(self,notifier,checkargs=True):
        """
        this registers a function to be called when the value changes or is 
        otherwise rendered invalid.  The notifier will be called as
        notifier(oldvalobj,newvalobj)
        """
        from weakref import ref
        
        if not callable(notifier):
            raise TypeError('notifier not a callable')
        if checkargs:
            import inspect
            if len(inspect.getargspec(notifier)) == 2:
                raise TypeError('notifier does not have 2 arguments')
        if self._notifywrs is None:
            self._notifywrs = []
        self._notifywrs.append(ref(notifier))
    
    def __getitem__(self,key):
        if type(key) is int:
            return self._vals[key]
        else:
            if key is None:
                key = Source(None)
            elif isinstance(key,basestring):
                if 'derived' in key:
                    key = key.replace('derived','')
                    if key == '':
                        der = 0
                    else:
                        der = int(key)
                    ders = [v for v in self._vals if isinstance(v,DerivedValue)]
                    if len(ders) <= der:
                        raise IndexError('field has only %i DerivedValues' % len(ders))
                    return ders[der]
                key = Source(key)
            if isinstance(key,Source):
                for v in self._vals:
                    if v.source is key:
                        return v
                raise KeyError('Field does not have %s'%key)
            else:
                raise TypeError('key not a Source key or index')
            
    def __setitem__(self,key,val):
        if type(key) is int or key in self:
            i = key if type(key) is int else self._vals.index(self[key])
            val = self._checkConvInVal(val,self._vals[i].source)
            self._vals[i] = val
        else:
            if isinstance(key,Source):
                s = key
            elif isinstance(key,basestring):
                s = Source(key)
            elif key is None:
                s = None
            else:
                raise TypeError('specified key not a recognized Source')
            val = self._checkConvInVal(val if s is None else ObservedValue(val,s))
            self._vals.append(val)
        
    def __delitem__(self,key):
        if type(key) is not int: 
            key = self._vals.index(self[key])
        if i == 0 and self._notifywrs is not None:
            self._notifyValueChange(self._vals[0],self._vals[1] if len(self._vals)>1 else None)
        del self._vals[key]
            
    def insert(self,key,val):
        val = self._checkConvInVal(val)
        
        if type(key) is int:
            i = key
        else:
            i = self._vals.index(self[key])
        if i == 0 and self._notifywrs is not None:
            self._notifyValueChange(val,self._vals[0] if len(self._vals)>0 else None)
        self._vals.insert(i,val)
        
#change to Sequnce,Mapping type

#check tuple-base setting
    @property
    def name(self):
        return self._name
    
    def _getType(self):
        return self._type
    def _setType(self,newtype):
        if newtype is None:
            self._type = None
        else:
            for v in self._vals:
                v.checkType(newtype)
            self._type = newtype
    type = property(_getType,_setType,doc="""
    Selects the type to enforce for this field.  
    if None, no type-checking will be performed
    if a numpy dtype, the value must be an array matching the dtype
    """)
#    #TODO:default should be Catalog-level?    
    def _getDefault(self):
        return self[None].value
    def _setDefault(self,val):
        self[None] = ObservedValue(val,None)
    def _delDefault(self):
        del self[None]
    default = property(_getDefault,_setDefault,_delDefault,"""
    The default value is the FieldValue that has a
    the None Source
    """)
    
    def _getCurr(self):
        try:
            return self._vals[0]
        except IndexError:
            raise IndexError('Field %s empty'%self._name)
    def _setCurr(self,val):
        oldcurr = self._vals[0] if len(self._vals)>0 else None
        try:
            i = self._vals.index(self[val])
            valobj = self._vals.pop(i)
        except (KeyError,IndexError,TypeError):
            valobj = self._checkConvInVal(val)
        self._vals.insert(0,valobj)
        self._notifyValueChange(oldcurr,valobj)
    currentobj = property(_getCurr,_setCurr)
    
    @property
    def values(self):
        return [v() for v in self._vals]
    
    @property
    def sources(self):
        return [str(v.source) for v in self._vals]
    
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
    
    def checkType(self,type):
        """
        ensure that the value of this FieldValue is of the requested Type
        (or None to accept anything).  Any mismatches will throw
        a TypeError
        """
        self._doTypeCheck(type,self.value)
        
    def _doTypeCheck(self,type,val):
        """
        handles interpretation of types - subclasses
        should call this with the val that should be checked
        if it is not the regular value
        """
        if type is not None:
            if isinstance(type,np.dtype):
                if not isinstance(val,np.ndarray):
                    raise TypeError('Value %s not a numpy array'%val)
                if self.value.dtype != type:
                    raise TypeError('Array %s does not match dtype %s'%(val,type))
            elif not isinstance(val,type):
                raise TypeError('Value %s is not of type %s'%(val,type))
    
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
    __slots__=('_f','_value','_valid')
    def __init__(self,f):
        import inspect
        
        if callable(f):
            self._f = f
        else:
            raise ValueError('attempted to make a DerivedValue with a non-callable')
        
        self._valid = False
        self._value = None
        
        raise NotImplementedError
    
    @property
    def value(self):
        raise NotImplementedError

class DependsSource(Source):
    __slots__=('dependson')
    
    def __init__(self):
        self._str = 'depend'
        self.dependson = []
        raise NotImplementedError
        
#class DerivedValue(FieldValue):
#    """
#    This value is derived from a set of other fields (possibly on other 
#    objects).  Currently it does not support cycles (where e.g. 
#    DerivedValue A depends on DerivedValue B which depends on A)
#    """
#    __slots__=('_f','_val','dependson')
#    def __init__(self,func,dependson=None,source=None):
#        """
#        The supplied function will be used with the sequence of Field's
#        that the values for the derived value are to work with. 
        
#        alternatively, if depndson is None, the  depndencies will be 
#        inferred from the default values of the 
#        """
#        super(DerivedValue,self).__init__()
        
#        from weakref import ref
        
#        if dependson is None:
#            from inspect import getargspec
#            args,varargs,varkw,defaults = getargspec(func)
#            if len(args) != len(defs) or varargs or varkw:
#                raise ValueError('input function does not have all defaults \
#                                  matched to dependencies or has varargs')
#            dependson = defaults 
#        #TODO:infer dependencies from argument names and parents/self fields
        
        
#        self.dependson = []
#        for dep in dependson:
#            if not isinstance(dep,Field):
#                raise ValueError('provided dependencies are not Fields')
#            self.dependson.append(ref(dep))
        
#        self._f = func
#        self._val = None
#        self.source = source
#        #TODO:more intelligently work out argument bindings
    
#    def __str__(self):
#        return '%s:Derived'%self.value
        
#    @property
#    def value(self):
#        if self._val is None:
#            self._val = self._f(*(fieldwr()() for fieldwr in self.dependson))
#            #TODO:check for bad wrs?
#        return self._val
    
#<------------------------------Node types------------------------------------->

class Catalog(CatalogNode):
    """
    This class represents a catalog of objects or catalogs.
    
    A Catalog is essentially a node in the object tree that 
    must act as a root.
    """    
    def __init__(self,name='default Catalog'):
        super(Catalog,self).__init__(parent=None)
        self.name = name
    
    @property
    def parent(self):
        return None    
    
    #these methods allow support for doing uniform mapping-like lookups over a catalog
    def __contains__(self,key):
        return hasattr(self,key)
    def __getitem__(self,key):
        return getattr(self,key)
    
    

class CatalogObject(FieldNode):
    """
    This class represents a FieldNode in the catalog that follows a particular
    data structure (i.e. a consistent set of Fields).  It is meant to be
    subclassed to define generic types of objects in the catalog.
    
    The fields and names are inferred from the class definition and 
    hence the class attribute name must match the field name.  Any 
    FieldValues present in the class objects will be ignored
    """
    #__metaclass__ = _CatalogObjectMeta
    
    def __init__(self,parent):
        import inspect
        
        super(CatalogObject,self).__init__(parent)
        self._altered = False
        
        #apply Fields from class into new object as new Fields
        for k,fi in inspect.getmembers(self.__class__,lambda x:isinstance(x,Field)): #TODO:faster way than lambda?
            if fi.name != k: #TODO: figure out if this can be done at "compile-time"
                raise KeyError('Name of Field (%s) does not match name in class attribute (%s)'%(fi.name,k))
            if None in fi:
                fobj = Field(fi.name,fi.type,fi.default, True)
            else:
                fobj = Field(fi.name,fi.type)
            setattr(self,k,fobj)
            self._fieldnames.append(k)
            
         
    
    @property
    def altered(self):
        """
        If True, the object no longer matches the specification given by the 
        class.  Note that this will remain True even if the offending fields
        are returned to their correct state.
        """
        return self._altered
    
    def revert(self):
        """
        Revert this object back to the standard Fields for  the class.
        Any deleted fields will be populated with the class Default Value
        any attributes that match the names of deleted Fields will be 
        overwritten
        
        TODO:test
        """
        import inspect
        import new
        
        #replace any deleted Fields with defaults and keep track of which should be kept
        fields=[]
        for k,fi in inspect.getmembers(self.__class__,lambda x:isinstance(x,Field)):
            fields.append(k)
            if not hasattr(self,k) or not isinstance(getattr(self,k),Field):
                if None in fi:
                    fobj = Field(fi.name,fi.type,fi.default, True)
                else:
                    fobj = Field(fi.name,fi.type)
                setattr(self,k,fobj)
                
        for k,v in inspect.getmembers(self,lambda x:isinstance(x,Field)):
            if k not in fields:
                delattr(self,k)
                
        self._fieldnames = fields
        self._altered = False
        self.addField = new.instancemethod(CatalogObject.addField,self,CatalogObject)
        self.delField = new.instancemethod(CatalogObject.delField,self,CatalogObject)
    
    def addField(self,field):
        self._altered = True
        self.addField = super(CatalogObject,self).addField
        self.addField(field)
        
    def delField(self,fieldname):
        self._altered = True
        self.delField = super(CatalogObject,self).delField
        self.delField(fieldname)
        
        
#<--------------------builtin catalog types------------------------------------>

class AstronomicalObject(CatalogObject):
    from .coords import AngularPosition
    
    def __init__(self,parent=None,name='default Name'):
        super(AstronomicalObject,self).__init__(parent)
        self.name.default=name
        
    _fieldorder = ('name','loc')
    name = Field('name',basestring)
    loc = Field('loc',AngularPosition)

class Test1(AstronomicalObject):
    
    def f(name,loc):
        return '%s+%s'%(name,loc)
    
del ABCMeta,abstractmethod,abstractproperty,MutableSequence,pi,division #clean up namespace
  
