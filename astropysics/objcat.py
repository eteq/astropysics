#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""

======
objcat
======

The :mod:`objcat` module contains classes and functions for generating catalogs
of objects. This catalog is unique in that it stores source information and can
dynamical deriv quantities that are dynamically updated as their sources are
changed. It also provides mechanisms for safely saving these object catalogss in
a database, and someday soon will include a web interface to interact with the
catalog via the internet.

Fundamentally, the catalog should be thought of as a tree in the computer
science sense, or a directed acyclic graph of :class:`FieldNode` objects.
Normally, a :class:`Catalog` object acts as the root.

Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.objcat
   :parts: 1

Module API
----------

"""

#TODO: seperate CatalogNode into a class that has children and one that doesnt
#TODO: methods to automatically plot parts of the data sets and save visualizations
#TODO: more stable persistence options (with atomic edits and histories?)
#TODO: a gui catalog viewer w/ plotting options as impleemnted above
#TODO: modules to also dynamically update via a web server

from __future__ import division,with_statement
from .constants import pi
import numpy as np
from collections import deque

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
        
class CycleWarning(Warning):
    """
    This warning indicates a cycle was detected in some graph-like structure
    """
    def __init__(self,message):
        super(CycleWarning,self).__init__(message)
        
class SourceDataError(Exception):
    """
    This exception indicates a problem occured while trying to retrieve 
    external Source-related data
    """
    def __init__(self,message):
        super(SourceDataError,self).__init__(message)



#<-------------------------Node/Graph objects and functions-------------------->
class CatalogNode(object):
    """
    This object is the superclass for all elements/nodes of a catalog with both
    parents and children.  
    
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
        Change the order of the children
        
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
    
    def addChild(self,node):
        """
        Adds a node as a child of this node.  Note that this will replace the
        added node's current parent (if any)
        """
        if not isinstance(node,FieldNode):
            raise ValueError('children must be FieldNodes')
        node.parent = self
        
    def removeChild(self,node):
        """
        Removes the requested child from this node, leaving it an orphan (i.e. 
        it's parent is None)  
        
        If child is not present, a ValueError will be raised
        """
        if node not in self:
            raise ValueError('requested child not present in this node')
        node.parent = None
        
    
    @property
    def nnodes(self):
        """
        this gives the number of total nodes at this point in the tree
        (including self - e.g. a leaf in the tree returns 1)
        """
        return sum([c.nnodes for c in self._children],1)
    
    def idstr(self):
        """
        a string uniquely identifying the object in its catalog heirarchy
        """
        if self.parent is None:
            return 'Root@'+hex(id(self))
        else:
            obj = self
            i = j = 0
            while obj.parent is not None:
                i+=1
                obj = obj.parent
            assert obj.parent is None,'Catalog has no root!'
            for c in self.parent.children:
                if self is c:
                    break
                j+=1
            return '{0}/{1}/{2}'.format(i,j,obj.idstr()) 
        
    
    def visit(self,func,traversal='postorder',filter=False,includeself=True):
        """
        This function walks through the object and all its children, executing
        func(:class:`CatalogNode`) and returning a list of the return values
        
        traversal is the traversal order of the tree - can be:
        
        * 'preorder'
        * 'postorder'
        * an integer indicating at which index the root should be evaluated
          (pre/post are 0/-1)
        * a float between -1 and 1 indicating where the root should be evaluated
          as a fraction
        * 'level'/'breathfirst' 
        * None:only visit this Node
        
        `filter` can be:
        
        * False: process and return all values
        * a callable: is called as g(node) and if it returns False, the node
          will not be processed nor put in the list (also ignores anything that
          returns None)
        * any other: if the node returns this value on processing, it will not
          be included in the returned list
                    
        if `includeself` is not True, the function will not visit the node 
        itself (only the sub-trees)
        """       
        if callable(filter):
            oldfunc = func
            def func(*args,**kwargs):
                if filter(args[0]):
                    return oldfunc(*args,**kwargs)
                else:
                    return None
            filterval = None
        else:
            filterval = filter
        
        if type(traversal) is int:
            retvals = []
            doroot = True
            for i,c in enumerate(self._children):
                if i == traversal and includeself:
                    retvals.append(func(self))
                    doroot = False
                retvals.extend(c.visit(func,traversal))
            if doroot and includeself:
                retvals.append(func(self))
        elif type(traversal) is float:
            retvals = []
            doroot = True
            travi = int(traversal*self._children)
            for i,c in enumerate(self._children):
                if i == travi and includeself:
                    retvals.append(func(self))
                    doroot = False
                retvals.extend(c.visit(func,traversal))
            if doroot and includeself:
                retvals.append(func(self))
        elif traversal == 'postorder': #None means postorder
            retvals = []
            for c in self._children:
                retvals.extend(c.visit(func,traversal))
            if  includeself:
                retvals.append(func(self))    
        elif traversal == 'preorder':
            retvals = self.visit(func,0,filter,includeself)
        elif traversal == 'level' or traversal == 'breadthfirst':            
            retvals=[]
            q = deque()
            if includeself:
                q.append(self)
            while len(q)>0:
                elem = q.popleft()
                retvals.append(func(elem))
                q.extend(elem._children)
        elif traversal is None:
            retvals = [func(self)]
        else:
            raise ValueError('unrecognized traversal type')
        
        if filterval is not False:
            retvals = [v for v in retvals if v is not filterval]
        return retvals
    
    def save(self,file,savechildren=True,**kwargs):
        """
        save this node as a file with the given name or file-like object
        
        if savechildren is True, the entire subtree will be saved, if False,
        just this node
        
        Note that the parent and everything above this point will NOT be
        saved (when reloaded the parent will be None)
        
        extra kwargs are passed into utils.fpickle
        """
        from .utils import fpickle

        oldpar = self._parent
        self._parent = None
        
        oldchildren = self._children
        if not savechildren:
            self._children = []
          
        try:
            fpickle(self,file,**kwargs)
        finally:
            self._parent = oldpar
            self._children = oldchildren
        
    #this is a staticmethod to keep both save and load methods in the same 
    #place - they are also available at the package level
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

#place save/load at module-level, as well
def save(node,file,savechildren=True):
    """
    save the specified node as a file with the given name or file-like object
    
    if savechildren is True, the entire subtree will be saved, if False,
    just this node
    
    Note that the parent and everything above this point will NOT be
    saved (when reloaded the parent will be None)
    """      
    return node.save(file,savechildren)
load = CatalogNode.load

class ActionNode(object):
    """
    This object is the superclass for nodes of a catalog that perform an action
    but do not store data and hence are not a part of the :class:`CatalogNode`
    heirarchy. Thus, they have a parent, but no children.
    
    This class is typically used for nodes that are associated with a
    :class:`Catalog` objects, as most other :class:`CatalogNode` objects cannot 
    store an :class:`ActionNode`. 
    
    Subclasses must call ``super(Subclass,self).__init__(parent,[name])`` in
    their :meth:`__init__` and override the :meth:`__call__` method.
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self,parent,name='default action node'):
        self._parent = None
        
        if parent is not None:
            self.parent = parent
            
        self.name = name
    
    def _getParent(self):
        return self._parent 
    def _setParent(self,val):            
        if self._parent is not None:
            self._parent._actchildren.remove(self)
        if not hasattr(val,'_actchildren'):
            raise AttributeError('Attempted to assign to parent that does not support action nodes')
        val._actchildren.append(self)
        self._parent = val
    parent=property(_getParent,_setParent)
    
    @abstractmethod
    def __call__(self,*args,**kwargs):
        raise NotImplementedError
        
    
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
    
    Further, setting an attribute (e.g. node['fieldname'] = value) has different 
    behavior depending on the type of the input value:
    
    * an integer i: sets the current value for that field to the ith entry in 
      the field
    * a string: sets the current value to the entry that has the provided source
      name
    * a tuple (string,value): sets the field entry for the source named by the
      string to the value, and sets that source as the current value
    * a tuple (string,(value,lerr[,uerr]): sets the field entry for the source 
      named by the string to the value and sets the errors for that entry, and 
      sets that source as the current value
    """
    __slots__=('_fieldnames',)
    
    @abstractmethod
    def __init__(self,parent,**kwargs):
        """
        Create a new fieldNode with the provided parent.
        
        kwargs are assigned as node values as self[kwargkey] = kwargvalue
        """
        #TODO: consider using CatalogNode.__init__?
        super(FieldNode,self).__init__(parent)
        self._fieldnames = []
        
        for k,v in kwargs.iteritems():
            self[k] = v
        
    def __getstate__(self):
        d = super(FieldNode,self).__getstate__()
        d['_fieldnames'] = self._fieldnames
        for n in self._fieldnames:
            val = getattr(self,n)
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
        if isinstance(field,basestring):
            field=Field(field)
        elif not isinstance(field,Field):
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
    
    def getFieldString(self,sep='\n'):
        """
        Generate a string with field:value pairs for all fields using
        sep as the seperator between the strings
        """
        strs = []
        for f in self.fields():
            try:
                strs.append(f.name+':'+str(f()))
            except IndexError:
                strs.append(f.name+' empty')
        return sep.join(strs)

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
        iserr = issrc = False
        if key not in self._fieldnames:
            try:
                if isinstance(key,basestring):
                    if key.endswith('_src'):
                        key = key[:-4]
                        issrc = True
                    elif key.endswith('_err'):
                        key = key[:-4]
                        iserr = True
                    else:
                        raise IndexError
                else:
                    key = self._fieldnames[key]
            except (IndexError,TypeError):
                raise IndexError('Field "%s" not found'%key)
        try:
            if iserr:
                return getattr(self,key).currenterror
            elif issrc:
                return getattr(self,key).currentsource
            else:
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
    
    @property
    def values(self):
        vals = []
        for f in self.fields():
            try:
                vals.append(f.currentobj.value)
            except IndexError:
                vals.append('')
        return tuple(vals)
    
    @property
    def sourcenames(self):
        vals = []
        for f in self.fields():
            try:
                vals.append(str(f.currentobj.source))
            except IndexError:
                vals.append('')
        return tuple(vals)
    
    @property
    def fielddict(self):
        vals = {}
        for f in self.fields():
            try:
                vals[f.name] = (str(f.currentobj.source),f.currentobj.value)
            except IndexError:
                vals[f.name] = None
        return vals
    
    def setToSource(self,src,missing='skip'):
        """
        Sets a string or :class:`Source` object passed in as the `src` argument
        as the current source for the requested FieldNode and it's subtree.
        
        `missing` can be:
        
        * 'raise'/'exception'
            raise a ValueError if there is no value with that source
        * 'warn'
            give a warning if the value is missing
        * 'skip'
            do nothing 
        
        """
        if not isinstance(src,Source) and src != 'derived':
            src = Source(src)      #TODO:see if this is speedier          
            
        if callable(missing):
            action = missing
        elif missing == 'skip':
            def action(src,f,node):
                return None
        elif missing == 'raise' or missing == 'exception':
            def action(src,f,node):
                raise ValueError('could not find src %s in field %s and node %s'%(src,f,node))
        elif missing == 'warn':
            from warnings import warn
            def action(src,f,node):
                warn('could not find src %s in field %s and node %s'%(src,f,node))
        else:
            raise ValueError('invalid missing action')
        
        for f in self.fields():
            if src in f:
                f.currentobj = src
            else:
                action(src,f,self)
                    
    @staticmethod
    def setToSourceAtNode(node,src,missing='skip',traversal='postorder',siblings=False):
        """
        Sets a string or :class:`Source` object passed in as the `src` argument
        as the current source for the requested FieldNode and it's subtree.
        
        `missing` can be:
        
        * 'raise'/'exception'
            raise a ValueError if there is no value with that source
        * 'warn'
            give a warning if the value is missing
        * 'skip'
            do nothing 
        
        `traversal` is the type of traversal to perform on the subtree. See
        :meth:`FieldNode.visit` for values.
        
        If `siblings` is True, the function applies to all the children of this
        :class:`FieldNode's<FieldNode>` parent.
        """
        if isinstance(node,FieldNode):
            if siblings and node.parent is not None:
                node = node.parent
                includeself = False
            else:
                includeself = True
                
            def f(node):
                if hasattr(node,'setToSource'):
                    node.setToSource(src,missing)
                    return node
            return node.visit(f,traversal=traversal,includeself=includeself,filter=None)
        else: #assume iterable
            vals = []
            for n in node:
                vals.extend(FieldNode.setToSourceAtNode(n,src,missing,traversal,siblings))
            return vals
        
    def setToSourceAtSelf(self,*args,**kwargs):
        """
        Sets the given src string or Source object as the current source for 
        this FieldNode and it's subtree
        
        see :meth:`FieldNode.setToSourceAtNode` for details
        """
        return FieldNode.setToSourceAtNode(self,*args,**kwargs)
        
        
    def getFieldValueNodes(self,fieldname,value):
        """
        Searches through the tree and returns all nodes for which a particular 
        field name has the requested value
        """
        return FieldName.getFieldValueNodesAtNode(self,fieldname,value)
    
    @staticmethod
    def extractFieldAtNode(node,fieldnames,traversal='postorder',filter=False,
                           missing=0,converter=None,sources=False,errors=False,
                           asrec=False,includeself=True):
        """
        This will walk through the tree starting from the :class:`FieldNode`
        given by the `node` argument and generate an :class:`array
        <numpy.ndarray>` of values for the specified fieldname.
        
        `fieldnames` can be a single field name, a sequence of fieldnames, or a
        comma-seperated string of field names.
        
        `traversal` and `filter` accept arguments like :meth:`CatalogNode.visit`
        
        `missing` determines the behavior in the event that a field is not 
        present (or a non :class:`FieldNode` is encountered) it can be:
        
        * 'exception'
            raise an exception if the field is missing
        * 'skip'
            do not include this object in the final array
        * 'mask'
            put 0 in the array location and return a mask of missing values as
            the second return value
        * 'masked'
            return numpy.ma.MaskedArray objects instead of arrays, with the
            missing values masked
        * any other
            fill any missing entries with this value
        
        `converter` is a function that is applied to the data before being 
        added to the array (or None to perform no conversion).
        
        `sources` determines if a sequence of sources should be returned - if
        False, only the array is returned, if it evaluates to True, the sources
        will be returned as described below.  By default, string representations
        of the sources are returned - if `sources`=='object', the actual
        :class:`Source` objects will be returned instead.
        
        `dtypes` determines the numpy data type of the output arrays - if None, 
        they will be inferred from the data. Otherwise, it must be either a 
        sequence of size matching the number of fields requested, an equivalent 
        comma-seperated string, or a dictionary mapping fieldnames to dtypes
        
        If `errors` is True, the errors will be returned (see return values
        below).
        
        Is `asrec` is True, a :class:`record array <numpy.recarray>` is
        generated instead of a regular :class:`array <numpy.ndarray>`. It can
        optionally specify the numpy data type of the values, if it is a
        sequence of size matching the number of fields, an equivalent
        comma-seperated string, or a dictionary mapping fieldnames to dtypes.
        Otherwise, the dtypes are inferred from the arrays themselves.
        
        
        If `includeself` is True, the node itself will be included.
        
        *returns*
        
        * a :class:`record array<numpy.recarray>` if `asrec` is True
        * an f x N :class:`array<numpy.ndarray>` of values if `sources` and 
          `errors` are False
        * (values,errors,sources) if `sources` and `errors` are True
        * (values,sources) if `sources` are True and `errors` are False
        * (values,errors) if `errors` are True and `sources` are False
        
        If `missing` is 'mask', the above return values are wrapped in a
        (returnval,mask) tuple.
        """
        from functools import partial
        
        if missing in ('exception','raise','skip','mask','masked'):
            missingval = 0
        else:
            missingval = missing
            
        
        if isinstance(fieldnames,basestring):
            if ',' in fieldnames:
                fieldnames = fieldnames.split(',')
            else:
                fieldnames = [fieldnames]
                
        
        if converter is None:
            if missing == 'exception' or missing == 'raise':
                def visitfunc(node,fieldname):
                    try:
                        return node[fieldname]
                    except AttributeError,e:
                        args = list(e.args)
                        args[0] = "Node %s has no field '%s'"%(node.idstr(),fieldname)
                        e.args = tuple(args)
                        e.message = args[0]
                        raise
            else:
                def visitfunc(node,fieldname):
                    try:
                        return node[fieldname]
                    except (KeyError,IndexError,TypeError,AttributeError):
                        return missingval
        else:
            if missing == 'exception' or missing == 'raise':
                def visitfunc(node,fieldname):
                    try:
                        return converter(node[fieldname])
                    except AttributeError,e:
                        args = list(e.args)
                        args[0] = "Node %s has no field '%s'"%(node.idstr(),fieldname)
                        e.args = tuple(args)
                        e.message = args[0]
                        raise
            else:
                def visitfunc(node,fieldname):
                    try:
                        return converter(node[fieldname])
                    except (KeyError,IndexError,TypeError,AttributeError):
                        return missingval
                    
        def maskfunc(node,fieldname):
            try:
                node[fieldname]
                return True
            except (KeyError,IndexError,TypeError,AttributeError):
                return False 
            
        if filter is not False and not callable(filter):
            def maskfilter(node,fieldname):
                return visitfunc(node,fieldname)!=filter
        else:
            maskfilter = filter
            
        lsts = []
        masks = []
        maskfilterp = maskfilter
        for fn in fieldnames:
            lsts.append(node.visit(partial(visitfunc,fieldname=fn),traversal=traversal,filter=filter,includeself=includeself))
            if filter is not False and not callable(filter):
                maskfilterp = partial(maskfilter,fieldname=fn)                
            masks.append(node.visit(partial(maskfunc,fieldname=fn),traversal=traversal,filter=maskfilterp,includeself=includeself))
        
        if missing=='skip':
            lsts = [np.array(l)[np.array(m)] for l,m in zip(lsts,masks)]
        else:
            lsts = [np.array(l) for l in lsts]
        

        if sources:
            def srcfunc(node,fieldname):
                try:
                    return getattr(node,fieldname).currentobj.source
                except (KeyError,IndexError,TypeError,AttributeError):
                    return None 
            if filter is not False and not callable(filter):   
                srcs = [node.visit(partial(srcfunc,fieldname=fn),traversal=traversal,filter=partial(maskfilter,fieldname=fn),includeself=includeself) for fn in fieldnames]
            else:
                srcs = [node.visit(partial(srcfunc,fieldname=fn),traversal=traversal,filter=filter,includeself=includeself) for fn in fieldnames]
            
            if sources != 'object':
                srcs = [[str(s) for s in f]  for f in srcs] 
            
        if errors:
            def errfunc(node,fieldname):
                try:
                    return getattr(node,fieldname).currentobj.errors
                except (KeyError,IndexError,TypeError,AttributeError),e:
                    return (0,0)
            if filter is not False and not callable(filter):      
                errs = [node.visit(partial(errfunc,fieldname=fn),traversal=traversal,filter=partial(maskfilter,fieldname=fn),includeself=includeself) for fn in fieldnames]
            else:
                errs = [node.visit(partial(errfunc,fieldname=fn),traversal=traversal,filter=maskfilter,includeself=includeself) for fn in fieldnames]

        if asrec:
            from operator import isMappingType,isSequenceType
            
            if asrec is True:
                asrec = {}
            elif isinstance(asrec,basestring):
                asrec = asrec.split(',')
            
            if isMappingType(asrec):
                for fn in fieldnames:
                    if fn not in asrec:
                        asrec[fn] = None
            elif isSequenceType(asrec):
                if len(asrec) != len(fieldnames):
                    raise ValueError('asrec sequence does not match fieldnames')
                asrec = dict(zip(fieldnames,asrec))
                    
            else:
                raise ValueError('invalid asrec entry')
            
            
            ralists = []
            newfnms = []
            for i,(fnm,lst,msk) in enumerate(zip(fieldnames,lsts,masks)):
                newfnms.append(fnm)
                ralists.append(np.array(lst,dtype=asrec[fnm]))
                if errors:
                    errarr = np.array(errs[i],dtype=asrec[fnm])
                    newfnms.append(fnm+'_lerr')
                    ralists.append(errarr[:,0])
                    newfnms.append(fnm+'_uerr')
                    ralists.append(errarr[:,1])
                if sources:
                    newfnms.append(fnm+'_src')
                    ralists.append(srcs[i])
                if missing == 'mask':
                    newfnms.append(fnm+'_mask')
                    ralists.append(msk)
                    
            res = np.rec.fromarrays(ralists,names=newfnms)
            if missing == 'masked':
                return np.ma.MaskedArray(res,~np.array(msk))
            else:
                return res
            
        else:
            arr = np.array(lsts)
            if errors:
                errs = np.array(errs)
            
            if len(arr)==1:
                arr = arr[0]
                masks = masks[0]
                if sources:
                    srcs = srcs[0]
                if errors:
                    errs = errs[0]
                    
            if missing == 'mask':
                res = arr,np.array(masks)
            elif missing == 'masked':
                res = np.ma.MaskedArray(arr,~np.array(masks))
            else:
                res = arr
            
            if sources and errors:
                return (res,errs,srcs)    
            elif sources:
                return (res,srcs)
            elif errors:
                return (res,errs)
            else:
                return res
            
    def extractField(self,*args,**kwargs):
        """
        Walk through the tree starting from this object.
        
        If the kwarg `siblings` is set to True, this will also extract the
        siblings (e.g. the parent subtree, except for the parent itself), and
        this overwrites `includeself`.
        
        The other arguments are identical to
        :meth:`FieldNode.extractFieldAtNode`.
        """
        sib = kwargs.pop('siblings',False)
        if sib and self.parent is not None:
            if len(args)<6:
                kwargs['includeself'] = False
            else:
                args[5] = False
            return FieldNode.extractFieldAtNode(self.parent,*args,**kwargs)
        else:
            return FieldNode.extractFieldAtNode(self,*args,**kwargs)
    #extractField.__doc__ += extractFieldAtNode.__doc__
    
    @staticmethod
    def getFieldValueNodesAtNode(node,fieldname,value,visitkwargs={}):
        """
        Searches through the tree and returns all nodes for which a particular 
        field name has the requested value
        """
        def visitfunc(node):
            if hasattr(node,fieldname):
                v = getattr(node,fieldname)
                if callable(v) and v() == value:
                    return node
        visitkwargs['filter'] = None    
        return node.visit(visitfunc,**visitkwargs)
    
    

def generate_pydot_graph(node,graphfields=True):
    """
    this function will generate a pydot.Dot object representing the supplied 
    node and its children.  Note that the pydot package must be installed
    installed for this to work
    
    if `graphfields` is True, this includes the fields as a record style 
    graphviz graph
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
    Field (no arguments) will return the current value itself
    
    usedef specified if the default should be set -- if True, defaultval will 
    be used, if None, a None defaultval will be ignored but any other
    will be recognizd, and if False, no default will be set
    """
    __slots__=('_name','_type','_vals','_nodewr','_notifywrs','_units')
    
    def __init__(self,name,type=None,defaultval=None,usedef=None,units=None):
        """
        The field must have a name, and can optionally be given a type
        """        
        self._name = name
        self._vals = []
        self._type = type
        self._notifywrs = None
        self._nodewr = None
        self._units = units
        
        if usedef or (usedef is None and defaultval is not None):
            self.default = defaultval
    
    _okDVs = set() #DerivedValues that are safe to ignore, used by StructuredFieldNode
    def __getstate__(self):
        prunedvals = []
        for v in self._vals:
            if isinstance(v,DerivedValue):
                if v in Field._okDVs:
                    Field._okDVs.remove(v)
                else:
                    from warnings import warn
                    warn("can't pickle DerivedValue in %s"%self)
            else:
                prunedvals.append(v)
        return {'_name':self._name,'_type':self._type,'_vals':prunedvals}
        
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
    
    @property
    def currentsource(self):
        return self.currentobj.source
        
    @property
    def currenterror(self):
        o = self.currentobj
        if hasattr(o,'errors'):
            return o.errors
        else:
            return (0,0)
    
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
        "check and convert input value
        
        auto-converts tuples to ObservedValues or ObservedErroredValues
        taking the first element to be the source, and if the second
        is a tuple, it will be passed into ObservedErroredValue init
        #TODO: auto-convert callables with necessary information to derivedvalues
        
        * dosrccheck = True -> check if source is present and raise a  
          ValueError if so
        * dosrccheck = string/Source -> ensure that value matches specified 
          string/Source and if not raise a ValueError
        * dosrccheck = False -> do no checking - just convert
        """
        from operator import isSequenceType
        
        if not isinstance(val,FieldValue):
            if hasattr(val,'source') and hasattr(val,'value'):
                pass
            elif isSequenceType(val) and len(val)==2:
                if isinstance(val[1],tuple) and self.type is not tuple:
                    val = ObservedErroredValue(val[0],*val[1])
                else:
                    val = ObservedValue(val[0],val[1])
            else:
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
    
    def notifyValueChange(self,oldvalobject=None,newvalobject=None):
        """
        notifies all registered functions that the value in this 
        field has changed
        
        (see registerNotifier)
        """
        if self._notifywrs is not None:
            deadrefs=[]
            for i,wr in enumerate(self._notifywrs):
                callobj = wr()
                if callobj is None:
                    deadrefs.append(i)
                else:
                    callobj(oldvalobject,newvalobject)
                    
            
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
            if not isinstance(val,FieldValue) and key in self:
                val = (key,val)
            val = self._checkConvInVal(val,dosrccheck=self._vals[i].source)
            if i == 0:
                self.notifyValueChange(self._vals[0],val)
            self._vals[i] = val
        else:
            if isinstance(key,Source):
                s = key
            elif isinstance(key,basestring) or key is None:
                s = Source(key)
            else:
                raise TypeError('specified key not a recognized Source')
            if not isinstance(val,FieldValue):
                val = (s,val)
            val = self._checkConvInVal(val,dosrccheck=True)
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
        val = self._checkConvInVal(val,dosrccheck=True)
        
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
    
    #TODO: some sort of unit support beyond just names?
    def _getUnits(self):
        return self._units
    def _setUnits(self,val):
        self._units = val
    units = property(_getUnits,_setUnits)
    
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
        self[None] = val#ObservedValue(None,val)
    def _delDefault(self):
        del self[None]
    default = property(_getDefault,_setDefault,_delDefault,"""
    The default value is the FieldValue that has a the None Source -
    this property is the default value itself, not the fieldValue object
    """)
    
    def _getCurr(self):
        try:
            return self._vals[0]
        except IndexError:
            raise IndexError('Field %s empty'%self._name)
    def _setCurr(self,val):
        oldcurr = self._vals[0] if len(self._vals)>0 else None
        try:
            #TODO: actually check for cases instead of handling exceptions -> speed
            i = self._vals.index(self[val])
            valobj = self._vals.pop(i)
        except (KeyError,IndexError,TypeError):
            #convert input
            valobj = self._checkConvInVal(val,dosrccheck=False)
            #if it's already present, remove already existing value
            if valobj.source in self:
                del self[valobj.source]
        
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
    
    @property
    def errors(self):
        return [v.errors if hasattr(v,'errors') else (0,0) for v in self._vals]
    
class SEDField(Field):
    """
    This field represents the Spectral Energy Distribution of this object - 
    e.g. a collection of Spectra or Photometric measurements
    """
    
    
    __slots__ = ['_maskedsedvals','_unit']
    
    def __init__(self,name='SED',unit='angstroms',type=None,defaultval=None, usedef=None):
        from .spec import Spectrum,HasSpecUnits
        from .phot import PhotObservation
        
        super(SEDField,self).__init__(name)
        
        
        
        self._type = tuple((Spectrum,PhotObservation))
        self._maskedsedvals = set()
        unittuple = HasSpecUnits.strToUnit(unit)
        self._unit = unittuple[0]+'-'+unittuple[1]
        
        if type is not None and set(type) != set(self._type):
            raise ValueError("SEDFields only accept Spectrum and PhotObservation objects - can't set type")
        
        if defaultval:
            self[None] = defaultval
        
        
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
        from .phot import PhotObservation
        
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
    
    @property
    def specs(self):
        from .spec import Spectrum
        return [o for i,o in enumerate(self.values) if isinstance(o,Spectrum) if i not in self._maskedsedvals]
    @property
    def specsources(self):
        from .spec import Spectrum
        return [self.sources[i] for i,o in enumerate(self.values) if isinstance(o,Spectrum) if i not in self._maskedsedvals]
    
    @property
    def phots(self):
        from .phot import PhotObservation
        return [o for i,o in enumerate(self.values) if isinstance(o,PhotObservation) if i not in self._maskedsedvals] 
    @property
    def photsources(self):
        from .phot import PhotObservation
        return [self.sources[i] for i,o in enumerate(self.values) if isinstance(o,PhotObservation) if i not in self._maskedsedvals] 
    
    def _getUnit(self):
        return self._unit
    def _setUnit(self,val):
        from .spec import HasSpecUnits
        #this checks to make sure the unit is valid
        val = HasSpecUnits.strToUnit(val)
        val = val[0]+'-'+val[1]
        
        oldu = self._unit
        try:
            for obj in self:
                if hasattr(obj,'unit'):
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
        
    def getBand(self,bands,asflux=False,asdict=False):
        """
        determines the magnitude or flux in the requested band.
        
        The first photometric measurement in this SEDField that has
        the band will be used - if not present, it will be computed 
        from  the first Spectrum with appropriate overlap.  If none
        of these are found, a ValueError will be raised
        
        if asflux is True, the result will be returned as a flux - 
        otherwise, a magnitude is returned.
        
        asdict returns a dictionary of results - otherwise, a sequence will
        be returned (or a scalar if only one band is requested)
        """
        from .phot import str_to_bands
        bands = str_to_bands(bands)
        vals = []
        for b,bn in zip(bands,[b.name for b in bands]):
            v = None
            for p in self.phots:
                if bn in p.bandnames:
                    i = p.bandnames.index(bn)
                    v = p.flux[i] if asflux else p.mag[i]
                    break
            if v is None:
                for s in self.specs:
                    if b.isOverlapped(s):
                        v = b.computeFlux(s) if asflux else b.computeMag(s)
                        break
            if v is None:
                raise ValueError('could not locate value for band '%bn)
            vals.append(v)
        
        if asdict:
            return dict([(b.name,v) for b,v in zip(bands,vals)])
        elif len(vals)==1:
            return vals[0]
        else:
            return vals
    
    def getFullSED(self):
        """
        the generates a astropysics.spec.Spectrum object that represents all 
        the information contained in this SEDField
        """
        from .spec import Spectrum
        
        x = np.array([])
        f = np.array([])
        e = np.array([])
        for s in self.specs:
            x = np.r_[x,s.x]
            f = np.r_[f,s.flux]
            e = np.r_[e,s.err]
        
        for p in self.phots:
            pf,pe = p.flux,p.fluxerr
            px,pw = p.getBandInfo()
            px = np.tile(px/pw,np.prod(px.shape)/len(p))
            px = px.reshape((len(p),np.prod(px.shape)/len(p)))
            pf = pf.reshape((len(p),np.prod(pf.shape)/len(p)))
            pe = pe.reshape((len(p),np.prod(pe.shape)/len(p)))
            x = np.r_[x,px.ravel()]
            f = np.r_[f,pf.ravel()]
            e = np.r_[e,pe.ravel()]
            
        return Spectrum(x,f,e,unit=self.unit)
        
    def plotSED(self,specerrs=True,photerrs=True,plotbands=True,colors=('b','g','r','r','k'),log='',clf=True):
        """
        Generates a plot of the SED of this object.
        
        colors is a tuple of colors as (spec,phot,specerr,photerr,other)
        """       
        from matplotlib import pyplot as plt
        from .spec import HasSpecUnits
        
        specs = self.specs
        phots = self.phots
        
        mxy1 = np.max([np.max(s.flux) for s in specs]) if len(specs) > 0 else None
        mxy2 = np.max([np.max(p.getFluxDensity(self.unit)[0]) for p in phots]) if len(phots) > 0 else None
        mxy = max(mxy1,mxy2)
        mny1 = np.min([np.min(s.flux) for s in self.specs]) if len(specs) > 0 else None
        mny2 = np.min([np.min(p.getFluxDensity(self.unit)[0]) for p in self.phots]) if len(phots) > 0 else None
        if mny1 is None:
            mny = mny2
        elif mny2 is None:
            mny = mny1
        else:
            mny = min(mny1,mny2)
        mxx1 = np.max([np.max(s.x) for s in self.specs]) if len(specs) > 0 else None
        mxx2 = np.max([np.max(p.getBandInfo(self.unit)[0]) for p in self.phots]) if len(phots) > 0 else None
        mxx = max(mxx1,mxx2)
        mnx1 = np.min([np.min(s.x) for s in self.specs]) if len(specs) > 0 else None
        mnx2 = np.min([np.min(p.getBandInfo(self.unit)[0]) for p in self.phots]) if len(phots) > 0 else None
        if mnx1 is None:
            mnx = mnx2
        elif mnx2 is None:
            mnx = mnx1
        else:
            mny = min(mnx1,mnx2)
        
        preint = plt.isinteractive()
        try:
            plt.interactive(False)

            if clf:
                plt.clf()
                
            if 'x' in log and 'y' in log:
                plt.loglog()
            elif 'x' in log:
                plt.semilogx()
            elif 'y' in log:
                plt.semilogy()
            
            c = (colors[0],colors[2],colors[4],colors[4])
            for s in specs:
                s.plot(fmt='-',ploterrs=specerrs,colors=c,clf=False)
                
            lss = ('--',':','-.','-')
            for i,p in enumerate(phots):
                if plotbands:
                    if plotbands is True:
                        plotbands = {'bandscaling':(mxy - mny)*0.5,'bandoffset':mny}
                    plotbands['ls'] = lss[i%len(lss)]
                p.plot(includebands=plotbands,fluxtype='fluxden',unit=self.unit,clf=False,fmt='o',c=colors[1],ecolor=colors[3])
                
            rngx,rngy=mxx-mnx,mxy-mny
            plt.xlim(mnx-rngx*0.1,mxx+rngx*0.1)
            plt.ylim(mny-rngy*0.1,mxy+rngy*0.1)
            
            xl = '-'.join(HasSpecUnits.strToUnit(self.unit)[:2])
            xl = xl.replace('wavelength','\\lambda')
            xl = xl.replace('frequency','\\nu')
            xl = xl.replace('energy','E')
            xl = xl.replace('angstrom','\\AA')
            xl = xl.replace('micron','\\mu m')
            xl = tuple(xl.split('-'))
            plt.xlabel('$%s/{\\rm %s}$'%xl)
            
            plt.ylabel('$ {\\rm Flux}/({\\rm erg}\\, {\\rm s}^{-1}\\, {\\rm cm}^{-2} {\\rm %s}^{-1})$'%xl[1])
            
            if preint:
                plt.show()
                plt.draw()
        finally:
            plt.interactive(preint)
    
class _SourceMeta(type):
    def __call__(cls,*args,**kwargs):
        obj = type.__call__(cls,*args,**kwargs)
        if obj._str in Source._singdict:
            singobj = Source._singdict[obj._str]
            ol = obj.location
            sl = singobj.location
            if ol is not None and ol != sl:
                from warnings import warn
                warn('overwriting location %s with %s in %s'%(sl,ol,singobj))
                singobj._adscode = obj._adscode
        else:
            Source._singdict[obj._str] = obj
            
        return Source._singdict[obj._str]

class Source(object):
    """
    A source for an observation/measurement/value.  Note that there is always 
    only one instance if a source at a given time - any two Sources with the
    same source string are the same object
    
    The source can optionally include a URL location to look up metadata 
    like authors, publication date, etc (location property)
    
    the constructor string can be of the form 'str/loc' in which case
    loc will be interpreted as the location, if it is not specified
    in the argument.  If it is 'str//loc', the loc will not be validated
    (e.g. it is assumed to be a correct ADS abstract code)
    """
    __metaclass__ = _SourceMeta
    __slots__=['_str','_adscode','__weakref__']
    
    from weakref import WeakValueDictionary
    _singdict = WeakValueDictionary()
    del WeakValueDictionary
    
    def __init__(self,src,location=None):
        src = src._str if hasattr(src,'_str') else str(src)
        
        if location is None and '/' in src:
            srcsp = src.split('/')
            
            if srcsp[-2] == '': #don't do checking for // case
                self._str = '/'.join(srcsp[:-2]).strip()
                self._adscode = srcsp[-1].strip()

            else:
                self._str = '/'.join(srcsp[:-1]).strip()
                self.location = srcsp[-1].strip()
        else:
            self._str = src
            self.location = location
        
    def __reduce__(self):
        return (Source,(self._str+('' if self._adscode is None else ('//'+self._adscode)),))
        
    def __str__(self):
        return self._str + ((' @' + self.location) if self._adscode is not None else '')
    
    adsurl = 'adsabs.harvard.edu'
    
    def _getLoc(self):
        return self._adscode
    def _setLoc(self,val):
        if val is not None:
            if val == '':
                val = None
            else:
                val = self._findADScode(val)
        self._adscode = val
    location = property(_getLoc,_setLoc)
    
    @staticmethod
    def _findADScode(loc):
        from urllib2 import urlopen,HTTPError
        from contextlib import closing


        lloc = loc.lower()
        if 'arxiv' in lloc:
            url = 'http://%s/abs/arXiv:%s'%(Source.adsurl,lloc.replace('arxiv:','').replace('arxiv','').strip())
        elif 'astro-ph' in lloc:
            url = 'http://%s/abs/arXiv:%s'%(Source.adsurl,lloc.replace('astro-ph:','').replace('astro-ph','').strip())
        elif 'doi' in lloc: #TODO:check if doi is case-sensitive
            url = 'http://%s/doi/%s'%(Source.adsurl,lloc.replace('doi:','').replace('doi',''))
        elif 'http' in lloc:
            url = loc
        else: #assume ADS abstract code
            url = 'http://%s/abs/%s'%(Source.adsurl,loc)
            
        url += '?data_type=PLAINTEXT'
        try:
            with closing(urlopen(url)) as page: #raises exceptions if url DNE
                for l in page:
                    if 'Bibliographic Code:' in l: 
                        return l.replace('Bibliographic Code:','').strip()
        except HTTPError:
            raise SourceDataError('Requested location %s does not exist at url %s'%(loc,url))
        raise SourceDataError('Bibliographic entry for the location %s had no ADS code, or parsing problem'%loc)
    
    
    _adsxmlcache = {}
    
    @staticmethod
    def clearADSCache(adscode=None, disable=False):
        """
        this clears the cache of the specified adscode, or everything, if 
        the adscode is None
        """
        if adscode is None:
            Source._adsxmlcache.clear()
        else:
            del Source._adsxmlcache[adscode]
    
    @staticmethod
    def useADSCache(enable=True):
        """
        This is used to disable or enable the cache for ADS lookups - if the 
        enable argument is True, the cache is enable (or unaltered if it
        is already active)) and if it is False, it will be disabled
        
        note that if the cache is disabled, all entries are lost
        """
        if enable:
            if Source._adsxmlcache is None:
                Source._adsxmlcache = {}
        else:
            Source._adsxmlcache = None
            
    def _getADSXMLRec(self):
        adscode = self._adscode
        if adscode is None:
            raise SourceDataError('No location provided for additional source data')
        
        if Source._adsxmlcache is not None and adscode in Source._adsxmlcache:
            xmlrec = Source._adsxmlcache[adscode]
        else:
            from urllib2 import urlopen
            from contextlib import closing
            from xml.dom.minidom import parseString
            
            with closing(urlopen('http://%s/abs/%s>data_type=XML'%(Source.adsurl,adscode))) as f:
                xmld = parseString(f.read())
            
            recs = xmld.getElementsByTagName('record')
            if len(recs) > 1:
                raise SourceDataError('Multiple matching ADS records for code %s'%adscode)
            
            xmlrec = recs[0]
            if Source._adsxmlcache is not None: 
                Source._adsxmlcache[adscode] = xmlrec
        
        return xmlrec
        
    def getBibEntry(self):
        """
        returns a string with the BibTeX formatted entry for this source, retrieved from ADS
        (requires network connection)
        """
        from urllib2 import urlopen
        
        if self._adscode is None:
            raise SourceDataError('No location provided for additional source data')
        
        with closing(urllib2.urlopen('http://%s/abs/%s>data_type=BIBTEX'%(Source.adsurl,self._adscode))) as xf:
            return xf.read()
        
    @property
    def authors(self):
        """
        The author list for this Source as a list of strings
        """
        rec = self._getADSXMLRec()
        return [e.firstChild.nodeValue for e in rec.getElementsByTagName('author')]
    
    @property
    def title(self):
        """
        The publication title for this Source
        """
        rec = self._getADSXMLRec()
        es = rec.getElementsByTagName('title')
        if len(es) != 1:
            raise SourceDataError('Title not found for %s'%self)
        return es[0].firstChild.nodeValue
    
    @property
    def abstract(self):
        """
        The abstract for this Source
        """
        rec = self._getADSXMLRec()
        es = rec.getElementsByTagName('abstract')
        if len(es) != 1:
            raise SourceDataError('Abstract not found for %s'%self)
        return es[0].firstChild.nodeValue
    
    @property
    def date(self):
        """
        The publication date of this Source
        """
        rec = self._getADSXMLRec()
        es = rec.getElementsByTagName('pubdate')
        if len(es) != 1:
            raise SourceDataError('Publication date not found for %s'%self)
        return es[0].firstChild.nodeValue
    
    @property
    def adsabs(self):
        """
        The URL for the ADS abstract of this Source
        """
        rec = self._getADSXMLRec()
        for e in rec.getElementsByTagName('link'):
            if e.attributes['type'].value == 'ABSTRACT':
                urlnodes = e.getElementsByTagName('url')
                if len(urlnodes)==1:
                    return urlnodes[0].firstChild.nodeValue
                else:
                    return [n.firstChild.nodeValue for n in urlnodes]
        raise SourceDataError('Abstract URL not found for %s'%self)
    
    @property
    def keywords(self):
        """
        The keywords for this source, and the type of the keywords, if 
        present
        """
        rec = self._getADSXMLRec()
        kwn = rec.getElementsByTagName('keywords')[0]
        kws = [n.firstChild.nodeValue for n in kwn.getElementsByTagName('keyword')]
        try:
            type = kwn.attributes['type'].value
        except KeyError:
            type = None
        return kws,type
    
class FieldValue(object):
    """
    Superclass for values that belong to a field.
    """
    __metaclass__ = ABCMeta
    __slots__ = ('_source')
    
    @abstractmethod
    def __init__(self):
        self._source = None
    
    def __getstate__(self):
        return {'_source':self._source}
    def __setstate__(self,d):
        self._source = d['_source']
    
    value = abstractproperty(doc='The value stored by this :class:`FieldValue`')
    """
    The value stored by this :class:`FieldValue`
    """
    
    def _getSource(self):
        return self._source
    def _setSource(self,val):
        if not (val is None or isinstance(val,Source)):
            try:
                val = Source(val)
            except: 
                raise TypeError('Input source is not convertable to a Source object')
        self._source = val 
    source=property(_getSource,_setSource,doc='The :class:`Source` for this :class:`FieldValue`')
    
    def checkType(self,typetocheck):
        """
        ensure that the value of this FieldValue is of the requested Type
        (or None to accept anything).  Any mismatches will throw
        a TypeError
        """
        from .utils import check_type
        check_type(typetocheck,self.value)
            
    def __call__(self):
        return self.value
    
    def __str__(self):
        return 'Value %s'%self.value
    
class ObservedValue(FieldValue):
    """
    This value is an observed or otherwise measured value for the field
    with the associated Source.
    """
    __slots__=('_value')
    def __init__(self,source,value):
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
        return '%s:Value %s'%(self.source,self.value)
    
    @property    
    def value(self):
        return self._value

class ObservedErroredValue(ObservedValue):
    """
    This value is an observed value with errorbars
    
    The value is either (value,err) or (value,uperr,lowerr)
    """
    __slots__=('_upperr','_lowerr')
    def __init__(self,source,value,upperr=None,lowerr=None):
        super(ObservedErroredValue,self).__init__(source,value)
        if upperr is None and lowerr is not None:
            #symmetric errors should be stored internally on upperr
            upperr = lowerr
            lowerr = None
        self._upperr = upperr
        self._lowerr = lowerr
    def __getstate__(self):
        d = super(ObservedErroredValue,self).__getstate__()
        d['_upperr'] = self._upperr
        d['_lowerr'] = self._lowerr
        return d
    def __setstate__(self,d):
        super(ObservedErroredValue,self).__setstate__(d)
        self._upperr = d['_upperr']
        self._lowerr = d['_lowerr']
        
    def __str__(self):
        if self._lowerr is not None:
            return '%s:Value %s +%s -%s'%(self.source,self.value,self._upperr,self._lowerr)
        else:
            return '%s:Value %s +/- %s'%(self.source,self.value,self._upperr)
        
    def checkType(self,typetocheck):
        """
        ensure that the value and errors are of the requested Type
        (or None to accept anything).  Any mismatches will raise
        a TypeError
        """
        if typetocheck is not None: #fast skip
            from .utils import check_type
            check_type(typetocheck,self.value)

            if self._upperr is not None:
                check_type(typetocheck,self._upperr)
            if self._lowerr is not None:
                check_type(typetocheck,self._lowerr)
    
    @property
    def errors(self):
        u,l = self._upperr,self._lowerr
        if u is None and l is None:
            return 0,0
        elif l is None:
            return u,u
        #elif u is None: #impossible case due to __init__
        #    return l,l
        else:
            return u,l
    
    def hasSymmetricErrors(self):
        return self._upperr==self._lowerr
        
class DerivedValue(FieldValue):
    """
    A FieldValue that derives its value from a function of other
    :class:`FieldValues<FieldValue>`.
    
    The values to use as arguments to the function are initially set through the
    default values of the function, and can be either references to fields or
    strings used to locate the dependencies in the catalog tree, using the
    `sourcenode` argument as the current source (see :class:`DependentSource`
    for details of how dependent values are dereferenced)
    
    Alternatively, the initializer argument `flinkdict` may be a dictionary of
    link values that overrides the defaults from the function.
    
    if `ferr` is True, the return value of `f` should be a 3-sequence of the
    form (value,uppererror,lowererror), and the arguments of `f` will be passed
    as (depvalue,uerr,lerr) from the source :class:`FieldValues<FieldValue>`.
    
    The class attribute :attr:`failedvalueaction` determines the action to take
    if a an exception is encountered while deriving the value and can be:
    
    * 'raise'
        raise an exception (or pass one along)
    * 'warn'
        issue a warning when it occurs, but continue execution with the value
        returned as None
    * 'skip'
        the value will be returned as None but the value will be left invalid
    * 'ignore'
        the value will be returned as None and will be marked valid
    
    """
    __slots__=('_f','_ferr','_value','_errs','_valid','_fieldwr')
    #TODO: auto-reassign nodepath from source if Field moves
    failedvalueaction = 'raise' 
    
    def __init__(self,f,sourcenode=None,flinkdict=None,ferr=False):
        """
        `f` is the function to use for deriving the value - must have default
        field links for every argument if `flinkdict` is None. 
        
        `sourcenode` specifies the current node to allow for links to be strings
        instead of explicit references to :class:`Fields<Field>`.
        
        `flinkdict` maps argument names to links, overriding the default values
        of the arguments of `f`
        """
        import inspect
        
        if callable(f):
            self._f = f
            self._ferr = ferr
            
            args, varargs, varkw, defaults = inspect.getargspec(f)
            
            if varargs or varkw:
                raise TypeError('DerivedValue function cannot have variable numbers of args or kwargs')
            if flinkdict:
                #populate any function defaults if not given already
                if defaults is not None:
                    for a,d in zip(reversed(args),defaults):
                        if a not in flinkdict:
                            flinkdict[a] = d
                defaults = [flinkdict[a] for a in args if a in flinkdict]
                
            if defaults is None or len(args) != len(defaults) :
                raise TypeError('DerivedValue does not have enought initial linkage items')
        else:
            raise TypeError('attempted to initialize a DerivedValue with a non-callable')
        
        self._valid = False
        self._value = None
        self._errs = None
        
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
            return 'Derived Value: %s'%self.value
        except:
            return 'Derived Value: Underivable'
    
    @property
    def source(self):
        return self._source
    
    @property
    def flinkdict(self):
        from inspect import getargspec
        
        args = getargspec(self._f)[0]
        return dict([t for t in zip(args,self.source.depstrs)])
    
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
    
    def _invalidateNotifier(self,oldvalobj,newvalobj):
        return self.invalidate()
    
    
    __invcyclestack = deque()
    def invalidate(self):
        """
        This marks this derivedValue as incorrect
        """
        try:
            if self in DerivedValue.__invcyclestack:
                from warnings import warn
                if self.sourcenode is None:
                    cycleloc = self.idstr()
                else:
                    if 'name' in self.sourcenode:
                        cycleloc = 'Node '+ self.sourcenode['name']
                    else:
                        cycleloc = 'Node '+ self.sourcenode.idstr()
                warn('Setting a DerivedValue that results in a cycle at '+cycleloc,CycleWarning)
                DerivedValue.__invcyclestack.append(self)
                return
            else:
                DerivedValue.__invcyclestack.append(self)
                
            self._valid = False
            if self.field is not None:
                self.field.notifyValueChange(self,self)
        finally:
            DerivedValue.__invcyclestack.pop()
    
    @property
    def value(self):
        if self._valid:
            return self._value
        else:
            try: 
                deps = self._source.getDeps(geterrs = True)
                if self._ferr:
                    valres = self._f(*deps)
                    self._value = valres[0]
                    self._errs = (valres[1],valres[2])
                else:
                    depvals = [d[0] for d in deps]
                    self._value = val = self._f(*depvals)
                    if np.all([(d[1]==0 and d[2]==0) for d in deps]):
                        
                        #if errors are all zero in dependencies, set to zero straight
                        self._errs = (0,0)
                    else:
                        #compute (df/dx)dx for all variables and add in quadrature
                        uerrs = []
                        lerrs = []
                        for i,d in enumerate(deps):
                            depvals[i] += d[1]
                            uerrs.append(self._f(*depvals) - val)
                            depvals[i] -= d[1]+d[2]
                            lerrs.append(val - self._f(*depvals))
                            depvals[i] += d[2]
                        uerr = np.sum(np.power(uerrs,2))**0.5
                        lerr = np.sum(np.power(lerrs,2))**0.5
                        self._errs = (uerr,lerr)
                self._valid = True
            except (ValueError,IndexError,CycleError),e:
                if isinstance(e,CycleError) and ' at ' not in e.args[0]:
                    #TODO: remove this node locating if it is too burdensome?
                    if self.sourcenode is None:
                        cycleloc = ' at '+self.idstr()
                    else:
                        if 'name' in self.sourcenode:
                            cycleloc = ' at Node '+ self.sourcenode['name'] 
                        else:
                            cycleloc = ' at Node '+ self.sourcenode.idstr()
                        if self.field is not None:
                            cycleloc += ' Field ' + self.field.name
                    e.args = (e.args[0]+cycleloc,)
                    e.message = e.args[0]
                
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
    @property
    def errors(self):
        v = self.value #need to compute value to make sure errors are up-to-date
        return self._errs

class DependentSource(Source):
    """
    This class holds weak references to the :class:`Field`s that are necessary
    to generate values such as in :class:`DerivedValue`.
    
    This source must know the :class:`CatalogNode` that it is expected to
    inhabit to properly interpret string codes for fields. Otherwise, the input
    fields must be :class:`Field` objects.
    
    References to parent objects or (first) child objects may be performed
    achieved by putting '^' or '.', respectively in front of the field name
    string (this can be done for multiple layers by chaining characters).
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
    
    def __init__(self,depfields,pathnode,notifierfunc=None):
        from weakref import ref
        
        self._str = 'dependent%i'%DependentSource._instcount
        self._adscode = None
        self.depfieldrefs = depfieldrefs = []
        self.depstrs = depstrs = []
        self.pathnode = pathnode
        self.notifierfunc = notifierfunc
        
        for f in depfields:
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
    
    @property
    def location(self):
        return None
    
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
        targeted Field
        
        see :class:`DependentSource` docs for a description of the mini-language
        """
        #TODO:REIMPLEMENT!
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
                        raise ValueError('No parent matching "%s" found'%substr)
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
        all of the references or raises a :exc:`ValueError` if any of the
        strings cannot be dereferenced
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
        
    __depcyclestack = deque()
    def getDeps(self,geterrs=False):
        """
        gets the values of the dependent field.
        
        the return value is a list of the values in the fields if geterrs is
        False, or (value,upperr,lowerr) if it is True
        """
        if self in DependentSource.__depcyclestack:
            DependentSource.__depcyclestack.append(self)
            raise CycleError('Attempting to compute a DerivedValue that results in a cycle')
        else:
            DependentSource.__depcyclestack.append(self)
        try:
            fieldvals = [wr() for wr in self.depfieldrefs]    
            if None in fieldvals:
                fieldvals = self.populateFieldRefs()
            if geterrs:
                fvs = []
                for fi in fieldvals:
                    if hasattr(fi,'errors'):
                        es = fi.currenterror
                        if es is None:
                            fvs.append((fi(),0,0))
                        else:
                            fvs.append((fi(),es[0],es[1]))
                    else:
                        fvs.append((fi(),0,0))
                return fvs
            else:
                return [fi() for fi in fieldvals]
        finally:
            DependentSource.__depcyclestack.pop()
    
#<------------------------------Node types------------------------------------->
class Catalog(CatalogNode):
    """
    This class represents a catalog of objects or catalogs.
    
    A Catalog is essentially a node in the object tree that 
    is typically (although not necessarily) the root of the 
    catalog tree.
    
    Iterator access is over the children, as is container testing and 
    index access.  
    
    Attributes can also be accessed as catobj[attrname] (this allows 
    access to Catalog attributes in the same way as FieldNodes)
    """    
    def __init__(self,name='default Catalog',parent=None):
        super(Catalog,self).__init__(parent)
        self.name = name
        self._actchildren = []
        
    def __str__(self):
        return 'Catalog %s'%self.name  
    
    def __contains__(self,key):
        if isinstance(key,basestring):
            return hasattr(self,key)
        else:
            return key in self._children
    def __getitem__(self,key):
        if isinstance(key,int):
            return self._children[key]
        return getattr(self,key)
    def __len__(self):
        return len(self._children)
    def __iter__(self):
        return iter(self._children)
    
    def isRoot(self):
        return self._parent is None
    
    @property
    def actionchildren(self):
        return tuple(self._actchildren)
    
    def mergeNode(self,node,skipdup=False,testfunc=None):
        """
        merges the given FieldNode into the catalog by reassigning all of 
        it's children to this Catalog

        testfunc is a function that will be called on the node to
        be transferred, and if it is True, the node will be moved, or
        if not it will be left alone.  If it is None, this 
        test will be skipped
        
        if skipdup is True, any nodes that already have an object of the same
        name will not be transferred
        
        returns the number of objects added
        """
        if testfunc is None:
            testfunc = lambda n:True
        
        nms = [n['name'] for n in self.children]
        
        i = -1
        countoffset = 1
        for i,c in enumerate(node.children):
            if (skipdup and c['name'] in nms) or not testfunc(c):
                countoffset -= 1    
            else:
                c.parent = self
            
        return i+countoffset
    
    #<----------------------FieldNode children-specific ----------------------->
    def extractField(self,*args,**kwargs):
        """
        Walk through the tree starting from this object and generate an array 
        of the values for the requested fieldname if any :class:`FieldNode`s are in this
        :class:`Catalog`.
        
        The other arguments are identical to
        :meth:`FieldNode.extractFieldAtNode`.
        """
        kwargs['includeself']=False
        return FieldNode.extractFieldAtNode(self,*args,**kwargs)
    
    def getFieldNames(self):
        """
        Searches the Catalog and finds all field names present in the children
        """
        s = set()
        def visfunc(node):
            if hasattr(node,'fieldnames'):
                s.update(node.fieldnames)
        self.visit(visfunc)
        return s
    
    def getFieldValueNodes(self,fieldname,value):
        """
        Searches the Catalog and finds all objects with the requested name
        """
        return FieldNode.getFieldValueNodesAtNode(self,fieldname,value,{'includeself':False})
    
    def locateName(self,name):
        """
        Searches the Catalog and finds all objects with the requested name
        """
        return FieldNode.getFieldValueNodesAtNode(self,'name',name,{'includeself':False})
    
    
class _StructuredFieldNodeMeta(ABCMeta):
    #Metaclass is used to check at class creation-time that fields all match names
    def __new__(mcs,name,bases,dct):
        cls = super(_StructuredFieldNodeMeta,mcs).__new__(mcs,name,bases,dct)
        for k,v in dct.iteritems():
            if isinstance(v,Field) and k != v.name:
                raise ValueError('StructuredFieldNode class %s has conficting \
                      field names - Node attribute:%s, Field.name:%s'%(name,k,v.name))
            elif isinstance(v,tuple) and len(v)==2 and \
                 isinstance(v[1],Field) and k != v[1].name: 
                 #derived values in the class should also be checked
                raise ValueError('StructuredFieldNode class %s has conficting \
                      derived field names - Node attribute:%s, Field.name:%s'%(name,k,v[1].name))
        return cls


class StructuredFieldNode(FieldNode):
    """
    This class represents a :class:`FieldNode` in the catalog that follows a
    particular data structure (i.e. a consistent set of :class:`Fields<Field>`).
    It is meant to be subclassed to define generic types of objects in the
    catalog.
    
    The fields and names are inferred from the class definition and hence the
    class attribute name must match the field name. Any
    :class:`FieldValues<FieldValue>` present in the class objects will be
    ignored.
    
    The :attr:`checkonload` class attribute determines the action to take if an
    unaltered :class:`StructuredFieldNode` is loaded and found to be
    inconsistent with the current class definition. It may be:
    
    * False/None/0
        Do no checking.
    * 'warn'
        show a warning when the loaded object is inconsistent with the
        definition
    * 'raise'
        raise a ValueError when the loaded object is inconsistent with the
        definition
    * 'revert'
        silently revert to the current definition if inconsistencies are found
    * 'revertwarn'
        issue a warning, then revert
    
    """
    __metaclass__ = _StructuredFieldNodeMeta
    
    @staticmethod
    def __fieldInstanceCheck(x):
        return isinstance(x,Field) or (isinstance(x,tuple) and len(x) == 2 and isinstance(x[0],DerivedValue))
    
    def __init__(self,parent,**kwargs):
        import inspect
        
        super(StructuredFieldNode,self).__init__(parent) #kwargs processed below
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
                fobj = fi.__class__(fi.name,type=fi.type,defaultval=fi[None], usedef=True)
            else:
                fobj = fi.__class__(fi.name,type=fi.type)
            setattr(self,k,fobj)
            fobj.node = self
            
            if dv is not None:
                dvs.append((dv,fobj))
            
            self._fieldnames.append(k)
            
        for dv,fobj in dvs:
            fobj.insert(0,DerivedValue(dv._f,sourcenode=self,flinkdict=dv.flinkdict,ferr=dv._ferr))
            
        for k,v in kwargs.iteritems():
            self[k] = v
            
    def __getstate__(self):
        import inspect
        
        #locate derived values and store where they should be re-inserted
        currderind = {}
        for k,v in inspect.getmembers(self.__class__,self.__fieldInstanceCheck):
            if isinstance(v,tuple):
                n = v[1].name
                if n in self._fieldnames:
                    fi = getattr(self,n)
                    for dv in fi:
                        if hasattr(dv,'_f') and dv._f is v[0]._f:
                            currderind[n] = fi.index(dv)
                            Field._okDVs.add(dv)
                    
        d = super(StructuredFieldNode,self).__getstate__()
        d['_altered'] = self._altered
        d['currderind'] = currderind
        return d
    
    
    checkonload = 'warn'
    def __setstate__(self,d):
        #function spends a lot of time in getmembers - optimize if possible?
        import inspect
        
        self._altered = d['_altered']
        super(StructuredFieldNode,self).__setstate__(d)
        
        if StructuredFieldNode.checkonload:
            inconsistent = False
            fields = []
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
                    fi.insert(ind,DerivedValue(v[0]._f,self,v[0].flinkdict))
            
                    
            #consistency check to ensure that the object has the necessary fields
            if StructuredFieldNode.checkonload:
                fields.append(k)
                if (not inconsistent) and getattr(self,k) is getattr(self.__class__,k) or not isinstance(getattr(self,k),Field):
                    inconsistent = k
                
        #do consistency check to ensure that there are no extra fields that shouldn't be present
        if StructuredFieldNode.checkonload:
            for k,v in inspect.getmembers(self,lambda x:isinstance(x,Field)):
                if k not in fields:
                    inconsistent = str(k)+' not in '+str(fields)
                    break
                
            if inconsistent is not False:
                if 'raise' in StructuredFieldNode.checkonload:
                    raise ValueError("object %s is not consistent with it's class definition due to %s"%(self,inconsistent))
                
                if 'warn' in StructuredFieldNode.checkonload or StructuredFieldNode.checkonload is True:
                    from warnings import warn
                    warn("object %s is not consistent with it's class definition due to %s"%(self,inconsistent))
                    
                if 'revert' in StructuredFieldNode.checkonload:
                    self.revert()
                
    
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
            #check to make sure the field is present, and if not, add it back in
            if hasattr(self,k) is hasattr(self.__class__,k) or not hasattr(self,k) or not isinstance(getattr(self,k),Field):
                if None in fi:
                    fobj = fi.__class__(fi.name,fi.type,fi.default, True)
                else:
                    fobj = fi.__class__(fi.name,fi.type)
                setattr(self,k,fobj)
                fobj.node = self
                
                if dv is not None:
                    dvs.append((dv,fobj))
                    
        #check to make sure there are no extra fields
        for k,v in inspect.getmembers(self,lambda x:isinstance(x,Field)):
            if k not in fields:
                delattr(self,k)
        
        self._fieldnames = fields
        
        for dv,fobj in dvs:
            fobj.insert(0,DerivedValue(dv._f,self,dv.flinkdict))
        
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
    def derivedFieldFunc(f=None,name=None,type=None,defaultval=None,
                         usedef=None,units=None,ferr=False,**kwargs):
        """
        this method is to be used as a function decorator to generate a 
        field with a name matching that of the function.  Note that the
        function should NOT have self as the leading argument
        
        Arguments are the same as for the Field constructor, except that 
        leftover kwargs are passed in as links that override the 
        defaults of the function
        """
        if f is None: #allow for decorator arguments
            return lambda f:StructuredFieldNode.derivedFieldFunc(f,name,type,defaultval,usedef,units,ferr,**kwargs)
        else: #do actual operation
            if name is None:
                name = f.__name__
            fi = Field(name=name,type=type,defaultval=defaultval,usedef=usedef,units=units)
            
            dv = DerivedValue(f,sourcenode=None,flinkdict=kwargs,ferr=ferr)
            return dv,fi
        
        
def arrayToNodes(values,source,fields,nodes,errors=None,matcher=None,
                 converters=None,namefield=None,setcurr=True):
    """
    Applies values from an array to CatalogNodes.
    
    values must be a 2-dimensional array or a 1-dimensional structured array.
    If a 2D array, the first dimension should be the fields, while the
    second dimension should be the nodes to be created.
    
    source is a Source object or a string that will be converted to a Source
    object
    
    nodes is either a container of nodes to visit or a CatalogNode 
    in which all FieldNodes within will be used as the nodes to visit
    
    fields is either a sequence of field names (must be at least as long as
    the array's second index dimension) or a mapping from indecies to field
    names or structured array field names to FieldNode field names
    
    if errors is not None, it should either be an object matching values 
    in shape, or a tuple (upper,lower) and will be applied as the errors
    on the values
    
    matcher is a callable called as matcher(arrayrow,node) - if it returns 
    False, the array will be matched to the next node - if True, the array
    values will be applied to the node
    
    converters  are either a sequence of callables or a mapping from indecies 
    to callables or structured array field names to callables that will be
    called on each array element before being added to the FieldNode as 
    converter(value).  If errors are provided, the converter will be called
    as converter((value,uerr,lerr)). If converters is None, no converting 
    is performed.
    
    
    namefield is a field name that will be set to a unique code of the form
    'src-i' where i is element index of the array row (or None to apply 
    no name).  It can also be a 2-tuple (name,converter) where converter is
    a callable of the form converter(i) that should return the value to 
    apply to the field.
    """
    from operator import isSequenceType,isMappingType
    from inspect import getargspec
    
    if not isinstance(source,Source):
        source = Source(source)
     
    if isinstance(nodes,CatalogNode):
        if isinstance(matcher,basestring):
            traversal = matcher
            matcher = None
        else:
            traversal = 'preorder'
        nodes = nodes.visit(lambda n:n,traversal,lambda n:isinstance(n,FieldNode))
    else:
        nodes = list(nodes)
    
    if isinstance(fields,basestring):
            fields={0:fields}    
    elif not isMappingType(fields):
        fields = dict([t for t in enumerate(fields)])
    else:
        fields = dict(fields) #copy
       
    if converters is None:
        pass
    elif not isMappingType(converters):
        converters = dict([t for t in enumerate(converters)])
    else:
        converters = dict(converters) #copy
    
    array = np.array(values,copy=False).T #now first dimension is along nodes and second is fields
    if errors is None:
        lerrors = uerrors = None
    elif isinstance(errors,tuple):
        uerrors,lerrors = errors
    else:
        uerrors = errors
        lerrors = None
    
    if converters is None:
        converters = {}
    
    if array.dtype.names is not None and len(array.shape) == 1:
        #structured/record array
        if len(array.shape) != 1:
            raise ValueError('structured array must be 1d')
        
        nms = array.dtype.names
        for k in fields.keys():
            if k in nms:
                fields[nms.index(k)] = fields[k]
                del fields[k]
        for k in converters.keys():
            if k in nms:
                converters[nms.index(k)] = converters[k]
    elif len(array.shape)==1:
        #1d array
        if len(fields)!=1:
            raise ValueError('input array is 1D but multiple fields were given')
        array = array.reshape((array.size,1)) #make 2D for the rest
        if uerrors is not None:
            uerrors = np.array(uerrors,copy=False)
            uerrors = uerrors.reshape((1,uerrors.size)) #errors are dimension (field,node)
        if lerrors is not None:
            lerrors = np.array(lerrors,copy=False)
            lerrors = lerrors.reshape((1,lerrors.size)) #errors are dimension (field,node)
    elif array.dtype.names is None and len(array.shape) == 2:
        #2d array
        if array.shape[1] != len(fields):
            raise ValueError('first dimension of array does not match number of fields')
    else:
        raise ValueError('invalid input array')
    
    fieldseq = []
    for i in range(len(array[0])):
        if i not in fields:
            fieldseq.append(None)
        else:
            fi = fields[i]
            fieldseq.append(fi)
            del fields[i]
    fieldseq = tuple(fieldseq)
    
    if fields:
        raise ValueError('fields %s were not found in the array'%fields)
            
    convseq = []
    for i in range(len(array[0])):
        if i not in converters:
            convseq.append(lambda val:val)
        else:
            convseq.append(converters[i])
    
    convseq = tuple(convseq)
    
    twoargseq =[]
    for i,conv in enumerate(convseq):
        if conv is None:
            twoargseq.append(None)
        else:
            if not callable(conv):
                raise ValueError('non-callable converter for field %s'%fieldseq[i])
            else:
                aspec = getargspec(conv)
                if not aspec[1]:
                    if len(aspec[0])==1:
                        twoargseq.append(False)
                    elif len(aspec[0])==2:
                        twoargseq.append(True)
                    else:
                        raise ValueError('converter for field %s has wrong number of arguments'%fieldseq[i])
    
    if matcher is None and len(array) != len(nodes):
        raise ValueError('with no matcher, the number of nodes must match the size of the array')
            
    
    for i,a in enumerate(array):
        if matcher is None:
            n = nodes[0]
            del nodes[0]
        else:
            #find the node that a matcher matches, and remove it if it is matched to get log(N) run time
            for j,n in enumerate(nodes):
                if matcher(a,n):
                    break
            else:
                j = None
            if j is not None:
                del nodes[j]

        for fi in fieldseq:
            if fi is not None and fi not in n.fieldnames:
                n.addField(fi)
        
        for fiind,v in enumerate(a):
            if fieldseq[fiind] is not None:
                fi = getattr(n,fieldseq[fiind])
                if uerrors is not None and uerrors[fiind] is not None:
                    ue = uerrors[fiind][i]
                    v = (v,ue,ue if lerrors is None else lerrors[fiind][i])
                if twoargseq[fiind]:
                    fi[source] = convseq[fiind](v,a)
                else:
                    fi[source] = convseq[fiind](v)
                if setcurr:
                    fi.currentobj = source
        
def arrayToCatalog(values,source,fields,parent,errors=None,nodetype=StructuredFieldNode,
                   converters=None,filter=None,namefield=None,nameconv=None,
                   setcurr=True):
    """
    Generates a catalog of nodes from the array of data.  
    
    See ``arrayToNodes`` for values,source,fields, and covnerter arguments.
    
    parent is the object to use as the parent for all of the nodes (which will
    be returned, a string (in which case a Catalog object will be created
    and returned), or None (return value will be a sequence of nodes)
    
    nodetype is the class to use to create the nodes (usually a subclass of 
    ``StructuredFieldNode``)
    
    filter is a function that will be called on the array row and if it 
    returns True, a node will be created, and if False, that row will be 
    skipped
    
    namefield is the field to use to apply the name of the node from the 
    index of the array row.  nameconv can be:
    
    * None: name field will be set to '<sourcestr>-i'
    * a callable: name field will be set to nameconv(i)
    * a sequence of mapping: name field will be set as nameconv[i]
    
    if namefiled is None, no name will be applied

    """
    if isinstance(parent,basestring):
        parent = Catalog(parent)
    
    source = Source(source) 
    if namefield:
        if nameconv is None:
            srcstr = source._str.split('/')[0]
            nameconv = lambda i:srcstr+'-'+str(i)
        elif callable(nameconv):
            pass
        else:
            nameseq = nameconv
            nameconv = lambda i:nameseq[i]
        
    if filter:  
        def matcher(arrayrow,node):
            return filter(arrayrow)
    else:
        filter = lambda a:True
        matcher = None

    array = np.array(values,copy=False)
    
    nodes = []
    for i,a in enumerate(array.T):
        if filter(a):
            n = nodetype(parent=parent)
            nodes.append(n)
            if namefield:
                if namefield not in n:
                    n.addField(namefield)
                nfi = getattr(n,namefield)
                if len(nfi) == 1 and None in nfi:
                    del nfi[None]
                nfi[source] = nameconv(i)
    
    arrayToNodes(array,source,fields,nodes,errors=errors,converters=converters,
                 matcher=matcher,setcurr=setcurr)
    
    if parent is None:
        return nodes
    else:
        return parent
    
    
#<--------------------builtin catalog types------------------------------------>

class AstronomicalObject(StructuredFieldNode):
    from .coords import EquatorialCoordinates
    
    def __init__(self,parent=None,name='default Name'):
        super(AstronomicalObject,self).__init__(parent)
        self.name.default=name
        
    def __str__(self):
        return 'Object %s'%self.name()
        
    _fieldorder = ('name','loc')
    name = Field('name',basestring)
    loc = Field('loc',EquatorialCoordinates)
    sed = SEDField('sed')

class Test1(StructuredFieldNode):
    num = Field('num',(float,int),(4.2,1,2))
    num2 = Field('num2',(float,int),(5.6,1.5,0.3))
    
    @StructuredFieldNode.derivedFieldFunc
    def f(num='num',num2='num2'):
        return 2*num+1+num2
    
    @StructuredFieldNode.derivedFieldFunc(ferr=True)
    def f2(num='num',num2='num2'):
        num,unum,lnum = num
        val = 2/num
        uerr = unum*2*num**-2
        lerr = lnum*2*num**-2
        return val,uerr,lerr
    
    @StructuredFieldNode.derivedFieldFunc
    def f3(num='num',num2='num2'):
        return num2*num
    
class Test2(StructuredFieldNode):
    val = Field('val',float,4.2)
    
    @StructuredFieldNode.derivedFieldFunc(num='num')
    def d1(val='val',d2='d2'):
        return val+np.exp(d2)
    
    @StructuredFieldNode.derivedFieldFunc(num='num')
    def d2(d1='d1'):
        if d1 is not None:
            return np.log(d1)
    
def test_cat():
    c = Catalog()
    t1 = Test1(c)
    t2 = Test1(c)
    t2['num'] = ('testsrc1',7)
    
    subc = Catalog('sub-cat',c)
    ts1 = Test1(subc)
    ts1['num'] = ('testsrc2',8.5)
    ts2 = Test1(subc,num=('testsrc2',12.7),f=('testows',123.45))
    ts3 = Test1(subc)
    ts3['num'] = ('testsrc2',10.3)
    
    return c

def test_sed():
    from numpy.random import randn,rand
    from numpy import linspace
    from .spec import Spectrum
    from .phot import PhotObservation
    from .models import BlackbodyModel
    
    f = SEDField()
    scale = 1e-9
    
    f['s1'] = Spectrum(linspace(3000,8000,1024),scale*(randn(1024)/4+2.2),scale*rand(1024)/12)
    m = BlackbodyModel(T=3300)
    m.peak = 2.2
    x2 = linspace(7500,10000,512)
    err = randn(512)/12
    f['s2'] = Spectrum(x2,scale*(m(x2)+err),scale*err)
    
    f['o1'] = PhotObservation('BVRI',[13,12.1,11.8,11.5],.153)
    f['o2'] = PhotObservation('ugriz',randn(5,12)+12,rand(5,12)/3)
    
    return f

    
del ABCMeta,abstractmethod,abstractproperty,MutableSequence,pi,division #clean up namespace
  
