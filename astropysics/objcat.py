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
    from collections import Sequence,MutableSequence
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    class MutableSequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
    class Sequence(object):
        __slots__=('__weakref__',) #support for weakrefs as necessary
        
class CycleError(Exception):
    """
    This exception indicates a cycle was detected in some graph-like structure
    """
    def __init__(self,message):
        super(CycleError,self).__init__(message)



#<-------------------------Node/Graph objects and functions-------------------->
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
            
    def __getstate__(self):
        return {'_parent':self._parent,'_children':self._children}
    def __setstate__(self,d):
        self._parent = d['_parent']
        self._children = d['_children']
        
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
    
    def save(self,file,savechildren=True):
        """
        save the file name or file-like object
        
        savechildren means the children of this node will be saved
        
        Note that the parent and everything above this point will NOT be
        saved (when reloaded the parent will be None)
        """
        import cPickle
        
        oldpar = self._parent
        self._parent = None
        
        oldchildren = self._children
        if not savechildren:
            self._children = []
            
        try:
            if isinstance(file,basestring):
                #filename
                with open(file,'w') as f:
                    return cPickle.dump(self,f)
            else:
                return cPickle.load(file)
        finally:
            self._parent = oldpar
            self._children = oldchildren
        
    
    @staticmethod
    def load(file):
        """
        load the file name or file-like object
        """
        import cPickle
        if isinstance(file,basestring):
            #filename
            with open(file,'r') as f:
                return cPickle.load(f)
        else:
            return cPickle.load(file)
    
class FieldNode(CatalogNode,Sequence):
    """
    A node in the catalog that has Fields.  This is an abstract class that 
    must have its initializer overriden.
    
    Note that for these subclasses, attribute access (e.g. node.fieldname) 
    accesses the Field object, while mapping or sequence-style access 
    (e.g node['fieldname'] or node[1])  directly accesses the current value
    of the field (or None if there is no value).  This means that 
    iterating over the object will also give values.  To iterate over
    the Field objects, use the fields() method.
    """
    __slots__=('_fieldnames',)
    
    @abstractmethod
    def __init__(self,parent):
        super(FieldNode,self).__init__(parent)
        self._fieldnames = []
        
    def __getstate__(self):
        d = super(FieldNode,self).__getstate__()
        d['_fieldnames'] = self._fieldnames
        for n in self._fieldnames:
            val = getattr(self,n)
            if isinstance(val,DerivedValue):
                from warnings import warn
                warn("cannot pickle derived values that aren't part of a structure - skipping %s"%val._str)
            else:
                d[n] = val
        return d
    def __setstate__(self,d):
        super(FieldNode,self).__setstate__(d)
        self._fieldnames = d['_fieldnames']
        for n in self._fieldnames:
            fi = d[n]
            setattr(self,n,fi)
            fi.node = self
            
    def addField(self,field):
        if not isinstance(field,Field):
            raise ValueError('input value is not a Field')
        if field.name in self._fieldnames:
            raise ValueError('Field name "%s" already present'%field.name)
        setattr(self,field.name,field)
        if field.node is not None:
            raise ValueError('a Field can only reside in one Node')
        field.node = self
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
        
    def fields(self):
        """
        this yields an iterator over all of the Field objects (rather than 
        their values, as regular sequence access does)
        """
        for n in self._fieldnames:
            yield getattr(self,n)
            
    def __str__(self):
        return 'FieldNode with fields %s'%self._fieldnames
        
    def __cmp__(self,other):
        try:
            return cmp(list(self),list(other))
        except TypeError:
            return 1
        
    def __len__(self):
        return len(self._fieldnames)
    
    def __contains__(self,key):
        return key in self._fieldnames
        
    def __getitem__(self,key):
        if key not in self._fieldnames:
            try:
                key = self._fieldnames[key]
            except (IndexError,TypeError):
                raise IndexError('Field "%s" not found'%key)
        try:
            return getattr(self,key)()
        except IndexError: #field empty
            return None
    
    def __setitem__(self,key,val):
        if key not in self._fieldnames:
            try:
                key = self._fieldnames[key]
            except (IndexError,TypeError):
                raise IndexError('Field "%s" not found'%key)
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
        return FieldNode.extractFieldAtNode(self,*args,**kwargs)
    
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
                except (KeyError,IndexError,TypeError):
                    return None
        elif not missing:
            filter = False
            def vfunc(node):
                try:
                    return node[fieldname]
                except (KeyError,IndexError,TypeError):
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

def generate_pydot_graph(node,graphfields=True):
    """
    this function will generate a pydot.Dot object representing the supplied 
    node and its children.  Note that the pydot package must be installed
    installed for this to work
    
    graphfields includes the fields as a record style graphviz graph
    """
    import pydot
    
    nodeidtopdnode={}
    def visitfunc(node):
        pdnode = pydot.Node(id(node),label=str(node))
        if isinstance(node,FieldNode):
            if graphfields:
                pdnode.set_shape('record')
                fieldstr = '|'.join([f.strCurr().replace(':',':') for f in node.fields()])
                pdnode.set_label('"{%s| | %s}"'%(node,fieldstr))
            else:
                pdnode.set_shape('box')
        nodeidtopdnode[id(node)] = pdnode
        try:
            edge = pydot.Edge(nodeidtopdnode[id(node.parent)],pdnode)
        except KeyError:
            edge = None
        return pdnode,edge
    nelist = node.visit(visitfunc,traversal='preorder')
    
    g = pydot.Dot()
    g.add_node(nelist[0][0])
    for node,edge in nelist[1:]:
        g.add_node(node)
        g.add_edge(edge)
        
    return g

#<----------------------------node attribute types----------------------------->    
 
class Field(MutableSequence):
    """
    This class represents an attribute/characteristic/property of the
    FieldNode it is associated with.  It stores the current value
    as well as all the other possible values.

    The values, sources, and default properties will return the actual values 
    contained in the FieldValues, while currentobj and iterating 
    over the Field will return FieldValue objects.  Calling the 
    Field (no arguments) will return the current value
    
    usedef specified if the default should be set -- if True, defaultval will 
    be used, if None, a None defaultval will be ignored but any other
    will be recognizd, and if False, no default will be set
    """
    __slots__=('_name','_type','_vals','_nodewr','_notifywrs')
    
    def __init__(self,name,type=None,defaultval=None,usedef=None):
        """
        The field must have a name, and can optionally be given a type
        """        
        self._name = name
        self._vals = []
        self._type = None
        self._notifywrs = None
        self._nodewr = None
        
        self.type = type
        if usedef or (usedef is None and defaultval is not None):
            self.default = defaultval
            
    def __getstate__(self):
        return {'_name':self._name,'_type':self._type,'_vals':self._vals}
    def __setstate__(self,d):
        self._name = d['_name']
        self._type = d['_type']
        self._vals = d['_vals']
        self._notifywrs = None
        self._nodewr = None
        #notifiers should late-attach when values are first accessed, and  
        #the Node does it's own attaching
        
    def __call__(self):
        return self.currentobj.value
    def __len__(self):
        return len(self._vals)
    def __contains__(self,val):
        #TODO: optimize!
        try:
            self[val]
            return True
        except (KeyError,IndexError):
            return False
    
    def __str__(self):
        return 'Field %s:[%s]'%(self._name,', '.join([str(v) for v in self._vals]))
    
    def strCurr(self):
        """
        returns a string with the current value instead of the list of 
        values (the behavior of str(Field_obj)
        """
        try:
            return 'Field %s: %s'%(self.name,self())
        except IndexError:
            return 'Field %s empty'%self.name
    
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
                    if v.source == val.source:
                        raise ValueError('value with %s already present in Field'%v.source)
        
        val.checkType(self.type)
        
        if isinstance(val,DerivedValue):
            if val.field is not None:
                raise ValueError('DerivedValues can only reside in a single field for dependencies')
            val.field = self
        return val
    
    def notifyValueChange(self,oldval=None,newval=None):
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
                if callobj is None:
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
        notifier(oldvalobj,newvalobj) BEFORE the value change is finalized.
        """
        from weakref import ref
        
        if not callable(notifier):
            raise TypeError('notifier not a callable')
        if checkargs:
            import inspect
            if len(inspect.getargspec(notifier)[0]) == 2:
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
                    #TODO: replace with only 1 if too slow?
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
                    if v.source == key:
                        return v
                raise KeyError('Field does not have %s'%key)
            else:
                raise TypeError('key not a Source key or index')
            
    def __setitem__(self,key,val):
        if type(key) is int or key in self:
            i = key if type(key) is int else self._vals.index(self[key])
            val = self._checkConvInVal(val,self._vals[i].source)
            if i == 0:
                self.notifyValueChange(self._vals[0],val)
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
        if type(key) is int: 
            i = key
        else:
            i = self._vals.index(self[key])
            
        if i == 0 and self._notifywrs is not None:
            self.notifyValueChange(self._vals[0],self._vals[1] if len(self._vals)>1 else None)
        del self._vals[i]
            
    def insert(self,key,val):
        val = self._checkConvInVal(val)
        
        if type(key) is int:
            i = key
        else:
            i = self._vals.index(self[key])
        if i == 0 and self._notifywrs is not None:
            self.notifyValueChange(val,self._vals[0] if len(self._vals)>0 else None)
        self._vals.insert(i,val)
        
    @property
    def name(self):
        return self._name
    
    def _getNode(self):
        return None if self._nodewr is None else self._nodewr()
    def _setNode(self,val):
        if val is None:
            self._nodewr = None
        else:
            from weakref import ref
            self._nodewr = ref(val)
        for d in self.derived:
            d.sourcenode = val
    node = property(_getNode,_setNode,doc='the node to which this Field belongs')
    
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
    can also be a sequence of types (accepts all) or a function
    that will be called directly on the function that returns True if
    the type is valid
    """)
    #TODO:default should be Catalog-level?    
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
        self.notifyValueChange(oldcurr,valobj)
        self._vals.insert(0,valobj)
    currentobj = property(_getCurr,_setCurr)
    
    @property
    def values(self):
        return [v() for v in self._vals]
    
    @property
    def sources(self):
        return [v.source for v in self._vals]
    
    @property
    def sourcenames(self):
        return [str(v.source) for v in self._vals]
    
    @property
    def observed(self):
        """
        returns a list of all the ObservedValue objects except the default
        """
        return [o for o in self if (isinstance(o,ObservedValue) and o.source._str != 'None')]
        
    @property
    def derived(self):
        """
        returns a list of all the DerivedValue objects
        """
        return [o for o in self if isinstance(o,DerivedValue)]
    
class SEDField(Field):
    """
    This field represents the Spectral Energy Distribution of this object - 
    e.g. a collection of Spectra or Photometric measurements
    """
    from .spec import Spectrum
    from .phot import PhotObservation
    
    __slots__ = ['_maskedsedvals','_unit']
    
    def __init__(self,name='SED',unit='angstroms'):
        super(SEDField,self).__init__(name)
        
        self._type = tuple((Spectrum,PhotObservation))
        self._unit = unit
        self._maskedsedvals = set()
        
    def __getstate__(self):
        d = super(SEDField,self).__getstate__()
        d['_maskedsedvals'] = self._maskedsedvals
        d['_unit'] = self._unit
        return d
    
    def __setstate__(self,d):
        super(SEDField,self).__setstate__(d)
        self._maskedsedvals = d['_maskedsedvals']
        self._unit = d['_unit']
        
    def __call__(self):
        return self.getFullSED()
    
    def __setitem__(self,key,val):
        """
        allow self['source'] = (bands,values) syntax
        """
        if isinstance(val,tuple) and len(val) == 2:
            return super(SEDField,self).__setitem__(key,PhotObservation(*val))
        else:
            return super(SEDField,self).__setitem__(key,val)
        
    @property
    def type(self):
        return self._type
    
    @property
    def default(self):
        return self()
    
    def _getUnit(self):
        return self._unit
    def _setUnit(self,val):
        oldu = self._unit
        try:
            for obj in self:
                obj.unit = val
            self._unit = val
        except:
            for obj in self:
                obj.unit = oldu
            raise
    unit = property(_getUnit,_setUnit,doc="""
    The units to use in the objects of this SED - see 
    astropysics.spec.HasSpecUnits for valid units
    """)
    
    def getMasked(self):
        """
        return a copy of the values masked from the full SED
        """
        return tuple(sorted(self._maskedsedvals))
    def mask(self,val):
        """
        mask the value (either index or source name) to not appear in the
        full SED
        """
        if isinstance(val,int):
            self._maskedsedvals.add(val)
        else:
            self._maskedsedvals.add(self.index(val))
    def unmask(self,val):
        """
        unmask the value (either index or source name) to appear in the
        full SED
        """
        if isinstance(val,int):
            self._maskedsedvals.remove(val)
        else:
            self._maskedsedvals.remove(self.index(val))
    def unmaskAll(self):
        """
        unmask all values (all appear in full SED)
        """
        self._maskedsedvals.clear()
    
    def getFullSED(self):
        """
        the generates a astropysics.spec.Spectrum object that represents all 
        the information contained in this SEDField
        """
        for i,v in self.values:
            if i not in self._maskedsedvals:
                raise NotImplementedError
        
    def plotSED(self):
        """
        generates a plot of the SED of this object
        """
        raise NotImplementedError
    
class _SourceMeta(type):
    #TODO: improve Source Singleton concepts, replace with __eq__ support?
    def __call__(cls,*args,**kwargs):
        obj = type.__call__(cls,*args,**kwargs)
        if not obj._str in Source._singdict:
            Source._singdict[obj._str] = obj
        return Source._singdict[obj._str]

class Source(object):
    __metaclass__ = _SourceMeta
    _singdict = {}
    __slots__=['_str']
    
    def __init__(self,src):
        self._str = str(src)
        
    def __reduce__(self):
        return (Source,(self._str,))
        
    def __str__(self):
        return 'Source ' + self._str
    
#TODO: more Source information in subclasses

    
class FieldValue(object):
    __metaclass__ = ABCMeta
    __slots__ = ('_source')
    
    @abstractmethod
    def __init__(self):
        self._source = None
    
    def __getstate__(self):
        return {'_source':self._source}
    def __setstate__(self,d):
        self._source = d['_source']
    
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
    
    def checkType(self,typetocheck):
        """
        ensure that the value of this FieldValue is of the requested Type
        (or None to accept anything).  Any mismatches will throw
        a TypeError
        """
        if typetocheck is not None:
            from operator import isSequenceType
            if isinstance(typetocheck,type):
                self._doTypeCheck((typetocheck,),self.value)
            elif callable(typetocheck):
                if not typetocheck(self.value):
                    raise TypeError('custom function type-checking failed')
            elif isSequenceType(typetocheck):
                self._doTypeCheck(typetocheck,self.value)
            else:
                raise ValueError('invalid type to check')
        
    def _doTypeCheck(self,types,val):
        """
        handles interpretation of types - subclasses
        should call this with the val that should be checked
        if it is not the regular value
        """
        for type in types:
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
        return 'Value %s'%self.value
    
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
        
    def __getstate__(self):
        d = super(ObservedValue,self).__getstate__()
        d['_value'] = self._value
        return d
    def __setstate__(self,d):
        super(ObservedValue,self).__setstate__(d)
        self._value = d['_value']
        
    def __str__(self):
        return 'Value %s:%s'%(self.value,self.source)
    
    @property    
    def value(self):
        return self._value
        
class DerivedValue(FieldValue):
    """
    A FieldValue that derives its value from a function of other FieldValues
    
    the dependent values are initially set through the default values of 
    the function, and can be either references to fields or strings used to
    locate the dependencies in the catalog tree, using the sourcenode 
    argument as the current source (see DependentSource for 
    details of how dependent values are derefernced)
    
    fieldnotifier is a notifier function, like that of Field.notifyValueChange
    to be called when the value is invalidated
    
    failedvalueaction determines the action to take if an 
    problem is encountered in deriving the value and can be:
    *'raise':raise an exception (or pass one along)
    *'warn': issue a warning when it occurs, but continue execution with the 
    value returned as None
    *'skip':the value will be returned as None but the value
    will be left invalid
    *'ignore':the value will be returned as None and will be marked valid
    """
    __slots__=('_f','_value','_valid','_fieldwr')
    #TODO: auto-reassign nodepath from source if Field moves
    failedvalueaction = 'raise' 
    
    def __init__(self,f,sourcenode=None):
        import inspect
        
        if callable(f):
            self._f = f
            args, varargs, varkw, defaults = inspect.getargspec(f)
            
            if varargs or varkw:
                raise TypeError('DerivedValue function cannot have variable numbers of args or kwargs')
            if len(args) != len(defaults) or defaults is None:
                raise TypeError('DerivedValue function does not have defaults to provide initial linkage')
        else:
            raise TypeError('attempted to initialize a DerivedValue with a non-callable')
        
        self._valid = False
        self._value = None
        
        self._fieldwr = None
        self._source = DependentSource(defaults,sourcenode,self._invalidateNotifier)
        
    def __getstate__(self):
        d = super(DerivedValue,self).__getstate__()
        d['_value'] = self._value
        from pickle import PicklingError
        raise PicklingError("DerivedValue can't be pickled because it depends on a function")
        d['_f'] = self._f #TODO: find some way to work around the function pickling issue?
        return d
    def __setstate__(self,d):
        super(DerivedValue,self).__setstate__(d)
        self._value = d['_value']
        self._valid = False
        self._fieldwr = None
        self._f = d['_f'] #TODO: find some way to work around this?
        
    def __str__(self):
        try:
            return 'Derived value: %s'%self.value
        except:
            return 'Derived value: Underivable'
    
    @property
    def source(self):
        return self._source
    
    def _getNode(self):
        return self._source.pathnode
    def _setNode(self,val):
        self._source.pathnode = val
    sourcenode=property(_getNode,_setNode,doc='The current location in the Catalog tree')
    
    def _getField(self):
        return None if self._fieldwr is None else self._fieldwr()
    def _setField(self,val):
        if val is None:
            self._notifierwr = None
        else:
            from weakref import ref
            self._fieldwr = ref(val)
    field=property(_getField,_setField,doc='A function like Field.notifyValueChange')

    
    def checkType(self,typetocheck):
        oldaction = DerivedValue.failedvalueaction
        try:
            DerivedValue.failedvalueaction = 'skip'
            super(DerivedValue,self).checkType(typetocheck)
        finally:
            DerivedValue.failedvalueaction = oldaction
    checkType.__doc__ = FieldValue.checkType.__doc__
    
    def _invalidateNotifier(self,oldval,newval):
        return self.invalidate()
    
    
    __invcycleinitiator = None
    def invalidate(self):
        """
        This marks this derivedValue as incorrect
        """
        try:
            if DerivedValue.__invcycleinitiator is None:
                DerivedValue.__invcycleinitiator = self
            elif DerivedValue.__invcycleinitiator is self:
                raise CycleError('attempting to set a DerivedValue that results in a cycle')
            self._valid = False
            if self.field is not None:
                self.field.notifyValueChange(self,self)
        finally:
            DerivedValue.__invcycleinitiator = None
    
    @property
    def value(self):
        if self._valid:
            return self._value
        else:
            try:
                self._value = self._f(*self._source.getDeps())
                self._valid = True
            except (ValueError,IndexError),e:
                if self.failedvalueaction == 'raise':
                    if len(e.args) == 2 and isinstance(e.args[1],list):
                        fields = [self._f.func_code.co_varnames[i] for i in e.args[1]]
                        raise ValueError('Could not get dependent values for field%s %s'%('s' if len(fields)>1 else '',fields))
                    else:
                        raise
                elif self.failedvalueaction == 'warn':
                    from warnings import warn
                    if len(e.args) == 2 and isinstance(e.args[1],list):
                        fields = [self._f.func_code.co_varnames[i] for i in e.args[1]]
                        warn('Problem getting dependent values for field%s %s'%('s' if len(fields)>1 else '',fields))
                    else:
                        warn('Problem encountered while deriving value '+str(e))
                    self._value = None
                    self._valid = False
                elif self.failedvalueaction == 'skip':
                    self._value = None
                    self._valid = False
                elif self.failedvalueaction == 'ignore':
                    self._value = None
                    self._valid = True
                else:
                    raise ValueError('invalid failedvalueaction')
            
            return self._value

class DependentSource(Source):
    """
    This class holds weak references to the Field's that are
    necessary to generate values such as in DerivedValue. 
    
    This source must know the Node that it is expected to inhabit to properly
    interpret string codes for fields.  Otherwise, the input fields must be
    Field objects
    """
    
    __slots__ = ('depfieldrefs','depstrs','_pathnoderef','notifierfunc')
    _instcount = 0
    
    
    #depfieldrefs: weakrefs to FieldValue objects 
    #depstrs: strings that should be used to locate FieldValue objects 
    #_pathnoderef: weakref to the CatalogNode used to dereference 
    
    def __new__(cls,*args,**kwargs):
        obj = super(DependentSource,cls).__new__(cls)
        DependentSource._instcount += 1
        return obj
    
    def __noneer(self):
        return None
    
    def __init__(self,fields,pathnode,notifierfunc=None):
        from weakref import ref
        
        self._str = 'dependent%i'%DependentSource._instcount
        self.depfieldrefs = depfieldrefs = []
        self.depstrs = depstrs = []
        self.pathnode = pathnode
        self.notifierfunc = notifierfunc
        
        for f in fields:
            if isinstance(f,basestring):
                depstrs.append(f)
                depfieldrefs.append(self.__noneer)
            elif isinstance(f,Field):
                depstrs.append(None)
                depfieldrefs.append(ref(f))
                if notifierfunc is not None:
                    f.registerNotifier(notifierfunc)
            elif f is None:
                depstrs.append(None)
                depfieldsrefs.append(self.__noneer)
            else:
                raise ValueError('Unrecognized field code %s'%str(f))
    
    def __reduce__(self):
        #Only thing guaranteed are the strings
        return (DependentSource,(self.depstrs,None))
        
    def __len__(self):
        return len(self.depfieldrefs)
    
    def _getPathnode(self):
        return self._pathnoderef()
    def _setPathnode(self,val):
        from weakref import ref
        
        if val is None:
            self._pathnoderef = self.__noneer
        else:
            if not isinstance(val,CatalogNode):
                raise TypeError('attemted to set pathnode that is not a CatalogNode')
            self._pathnoderef = ref(val)
        #invalidadte weakrefs that are dereferenced
        for i,s in enumerate(self.depstrs):
            if s is None:
                self.depfieldrefs[i] = self.__noneer
    pathnode = property(_getPathnode,_setPathnode,doc='The CatalogNode for dereferencing source names')
    
    @staticmethod
    def _locatestr(s,node):
        """
        this method translates from a string and a location to the actual
        targretted Field
        
        ^^^ means up that elements in the catalog, while ... means down 
        ^(name) means go up until name is found 
        """
#        if s in node.fieldnames:
#            return getattr(node,s)
#        else:
#            raise ValueError('Linked node does not have requested field "%s"'%self.depstrs[i])
        
        #TODO:optimize
        upd = {}
        for i,c in enumerate(s):
            if c is '^':
                d[i] = True
            if c is '.':
                d[i] = False
        if len(upd) == 0:
            return getattr(node,s)
        if len(s)-1 in upd:
            raise ValueError('Improperly formatted field string - no field name')
        
        pairs = []
        sortk = sorted(upd.keys())
        lasti = sortk[0]
        for i in sortk[1:]:
            pairs.append((lasti,i))
            lasti = i
        if len(pairs) == 0:
            pairs.append((sortk[0]-1,sortk[0]))
            
        for i1,i2 in pairs:
            if i2-i1 == 1:
                if upd[i1]:
                    node = node.parent
                else:
                    node = node.children[0]
            else:
                subs = s[i1:i2]
                if upd[i1]:
                    try:
                        node = node.parent
                        while node.__class__!=substr and node['name']!=substr:
                            node = node.parent
                    except AttributeError:
                        raise ValueError('No parrent matching "%s" found'%substr)
                else:
                    try:
                        nchild = int(subs)
                        node = node.children[nchild]
                    except ValueError:
                        startnode = node
                        for n in node.children:
                            if node.__class__==substr or node['name']==substr:
                                node = n
                                break
                        if node is startnode:
                            raise ValueError('No child matching "%s" found'%substr)
        
    
    def populateFieldRefs(self):
        """
        this relinks all dead weakrefs using the dependancy strings and returns 
        all of the references or raises a ValueError if any of the strings
        cannot be dereferenced
        """
        from weakref import ref
        
        if self.pathnode is None:
            refs = [wr() for wr in self.depfieldrefs]
            if None in refs:
                raise ValueError('Missing/dead field(s) cannot be dereferenced without a catalog location',[i for i,r in enumerate(refs) if r is None])
        else:
            if not hasattr(self.pathnode,'fieldnames'):
                raise ValueError('Linked pathnode has no fields or does not exist')
            
            refs = []
            
            for i,wrf in enumerate(self.depfieldrefs):
                if wrf() is None:
                    f = self._locatestr(self.depstrs[i],self.pathnode)
                    refs.append(f)
                    self.depfieldrefs[i] = ref(f)
                    if self.notifierfunc is not None:
                            f.registerNotifier(self.notifierfunc)
                else:
                    refs.append(wrf())
        
        return refs
        
    def getDeps(self):
        """
        get the values of the dependent fields
        """
        fieldvals = [wr() for wr in self.depfieldrefs]    
        if None in fieldvals:
            fieldvals = self.populateFieldRefs()
        return [fi() for fi in fieldvals]
    
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
        
    def __str__(self):
        return 'Catalog %s'%self.name 
    
    @property
    def parent(self):
        return None    
    
    #these methods allow support for doing uniform mapping-like lookups over a catalog
    def __contains__(self,key):
        return hasattr(self,key)
    def __getitem__(self,key):
        return getattr(self,key)
    
class _StructuredFieldNodeMeta(ABCMeta):
    #Metaclass is used to check at class creation-time that fields all match names
    def __new__(mcs,name,bases,dct):
        cls = super(_StructuredFieldNodeMeta,mcs).__new__(mcs,name,bases,dct)
        for k,v in dct.iteritems():
            if isinstance(v,Field) and k != v.name:
                raise ValueError('StructuredFieldNode class %s has conficting field names - Node attribute:%s, Field.name:%s'%(name,k,v.name))
        return cls

class StructuredFieldNode(FieldNode):
    """
    This class represents a FieldNode in the catalog that follows a particular
    data structure (i.e. a consistent set of Fields).  It is meant to be
    subclassed to define generic types of objects in the catalog.
    
    The fields and names are inferred from the class definition and 
    hence the class attribute name must match the field name.  Any 
    FieldValues present in the class objects will be ignored
    """
    __metaclass__ = _StructuredFieldNodeMeta
    
    @staticmethod
    def __fieldInstanceCheck(x):
        return isinstance(x,Field) or (isinstance(x,tuple) and len(x) == 2 and isinstance(x[0],DerivedValue))
    
    def __init__(self,parent):
        import inspect
        
        super(StructuredFieldNode,self).__init__(parent)
        self._altered = False
       
        dvs=[]  #derived values to apply to fields as (derivedvalue,field)
        #apply Fields from class into new object as new Fields
        for k,v in inspect.getmembers(self.__class__,self.__fieldInstanceCheck):
            if isinstance(v,tuple):
                dv,fi = v
            else:
                fi = v
                dv = None
            
            if None in fi:
                fobj = Field(fi.name,fi.type,fi.default, True)
            else:
                fobj = Field(fi.name,fi.type)
            setattr(self,k,fobj)
            fobj.node = self
            
            if dv is not None:
                dvs.append((dv,fobj))
            
            self._fieldnames.append(k)
            
        for dv,fobj in dvs:
            fobj.insert(0,DerivedValue(dv._f,self))
            
    def __getstate__(self):
        import inspect
        
        currderind = {}
        for k,v in inspect.getmembers(self.__class__,self.__fieldInstanceCheck):
            if isinstance(v,tuple):
                n = v[1].name
                if n in self._fieldnames:
                    fi = getattr(self,n)
                    while 'derived' in fi:
                        dv = fi['derived']
                        if dv._f is v[0]._f:
                            currderind[n] = fi.index(dv)
                        del fi[dv._source._str]
                    
        d = super(StructuredFieldNode,self).__getstate__()
        d['_altered'] = self._altered
        d['currderind'] = currderind
        return d
    def __setstate__(self,d):
        import inspect
        
        self._altered = d['_altered']
        super(StructuredFieldNode,self).__setstate__(d)
        for k,v in inspect.getmembers(self.__class__,self.__fieldInstanceCheck):
            if isinstance(v,tuple):
                n = v[1].name
                if n in self._fieldnames:
                    fi = getattr(self,n)
                    try:
                        ind = d['currderind'][n]
                    except KeyError:
                        from warnings import warn
                        warn('missing current index for structured derived value "%s" - assuming as default'%n)
                        ind = 0
                    if ind > len(fi):
                        ind = len(fi)
                    fi.insert(ind,DerivedValue(v[0]._f,self))
    
    @property
    def alteredstruct(self):
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
        import inspect,types
        
        dvs=[]  #derived values to apply to fields as (derivedvalue,field)
        #replace any deleted Fields with defaults and keep track of which should be kept
        fields=[]
        for k,v in inspect.getmembers(self.__class__,self.__fieldInstanceCheck):
            if isinstance(v,tuple):
                dv,fi = v
            else:
                fi = v
                dv = None
                
            fields.append(k)
            if not hasattr(self,k) or not isinstance(getattr(self,k),Field):
                if None in fi:
                    fobj = Field(fi.name,fi.type,fi.default, True)
                else:
                    fobj = Field(fi.name,fi.type)
                setattr(self,k,fobj)
                fobj.node = self
                
                if dv is not None:
                    dvs.append((dv,fobj))
        for k,v in inspect.getmembers(self,lambda x:isinstance(x,Field)):
            if k not in fields:
                delattr(self,k)
        
        self._fieldnames = fields
        
        for dv,fobj in dvs:
            fobj.insert(0,DerivedValue(dv._f,self))
        
        self._altered = False
        self.addField = types.MethodType(StructuredFieldNode.addField,self,StructuredFieldNode)
        self.delField = types.MethodType(StructuredFieldNode.delField,self,StructuredFieldNode)
    
    def addField(self,field):
        self._altered = True
        self.addField = super(StructuredFieldNode,self).addField
        self.addField(field)
        
    def delField(self,fieldname):
        self._altered = True
        self.delField = super(StructuredFieldNode,self).delField
        self.delField(fieldname)
    
    @staticmethod
    def derivedFieldFunc(f=None,type=None,defaultval=None,usedef=None):
        """
        this method is to be used as a function decorator to generate a 
        field with a name matching that of the function.  Note that the
        function should NOT have self as the leading argument
        """
        if f is not None: #do actual operation
            fi = Field(name=f.__name__,type=type,defaultval=defaultval,usedef=usedef)
            dv = DerivedValue(f,None)
            return dv,fi
        else: #allow for decorator arguments
            return lambda f:StructuredFieldNode.derivedFieldFunc(f,type,defaultval,usedef)
        
        
#<--------------------builtin catalog types------------------------------------>

class AstronomicalObject(StructuredFieldNode):
    from .coords import AngularPosition
    
    def __init__(self,parent=None,name='default Name'):
        super(AstronomicalObject,self).__init__(parent)
        self.name.default=name
        
    def __str__(self):
        return 'Object %s'%self.name()
        
    _fieldorder = ('name','loc')
    name = Field('name',basestring)
    loc = Field('loc',AngularPosition)

class Test1(AstronomicalObject):
    num = Field('num',float,4.2)
    
    @StructuredFieldNode.derivedFieldFunc(defaultval='f')
    def f(num='num'):
        return num+1
    
class Test2(AstronomicalObject):
    SED = SEDField()
    
def test_cat():
    c=Catalog()
    ao=AstronomicalObject(c,'group')
    t1=Test1(ao)
    t12=Test1(ao)
    t2=Test2(ao)
    t2.SED['src'] = ('BVRI',[12,11.5,11.3,11.2])
    
    return c,locals()

    
del ABCMeta,abstractmethod,abstractproperty,Sequence,MutableSequence,pi,division #clean up namespace
  
