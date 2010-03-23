#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""

=====
utils
=====

The ``utils`` module contains classes and funtions of general utility used in
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
    Unpickle a file specified by ``filename``, which can be either a
    file name string or a file object, and return the pickled contents.
    
    ``number`` specifies the number of objects to return -- if greater than 0,
    a list of length equal to number is returned.  If exactly 0, just a single 
    call to pickle.load is performed and the result is returned (e.g. the first 
    object in the file).  If <0, all objects in the file are retreived and a 
    list of all of them is returned.
    
    ``usecPickle`` determines if the cPickle module is to be used in place of
    pickle (generally cPickle is faster)
    """
    if Pickle:
        import cPickle as pickle
    else:
        import pickle
        
    if isinstance(fileorname,basestring):
        f = open(filename,'r')
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
    Pickle the ``object`` into the file given by ``fileorname``, which can be
    either a file name string or a file object
    
    if ``append`` is True, the object is appended to the end of the file, 
    otherwise the file will be replaced with a new one that the object will be
    pickled into.
    
    ``protocol`` determines which pickle protocol to use - see the pickle module
    for details on these options.  If None, the most recent protocol will be 
    used.
    
    ``usecPickle`` determines if the cPickle module is to be used in place of
    pickle (generally cPickle is faster)
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
    
    note that None will always be accepted if acceptnone is True
    
    types can be a single type, a sequence of types, or a callable that 
    accepts one argument and will be called with the value - if it 
    returns True, the value is the correct type. Otherwise, a TypeError
    will be raised.
    
    if the type is a numpy dtype, the input array type is also checked
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
    (Not thread safe?)
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
    mapping (Not thread safe?)
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
    A class to register data sets used throughout a module  (e.g. photometric 
    bands)
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
        register a set of objects, possibly in a group
        
        `objects` is a dictionary mapping the datum name to the object
        
        `groupname` is the name of the group to apply to the dictionary
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
        adds the registered object with the provided name to a group with the 
        given name
        """
        obj = self[key] 
        group = self._groupdict.setdefault(group,[])
        if key not in group:
            group.append(key)
            
    def removeGroup(group):
        """
        removes the selected group
        """
        del self._groupdict[group]
        
    def getGroupData(self,groupname):
        return self._groupdict[groupname][:]
    
    def getObjects(self,objectstrs,addmissing=False):
        """
        this translates a string, sequence of strings, or mixed sequence of data 
        objects and strings into a sequence of data objects 
        
        if addmissing is True, new objects will be added to the registry -- if the
        object has a `name` attribute, it will be used as the name, otherwise a name
        will be generated of the form "dataname##".  if addmissing is False, the 
        object will be returned but not added to the registry
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
        generates and returns  list of names of the objects in this registry, 
        extracted from the `name` attribute of the object if present, otherwise 
        from the registry key.  The order is the same as that generated by
        keys() or values()
        
        if retdict is True, the return value is a dictionary mapping keys from
        the registry to the object names
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
    


#<---------------------------data pipelining----------------------------------->
class PipelineError(Exception):
    """
    This exception is used to indicate problems with the pipeline.  It
    should NOT be used to indicate data matching problems within 
    PipelineElements -- ValueError or TypeError should be used for 
    that purpose
    """
    pass


class Pipeline(object):
    """
    This class represents a pipeline, composed of PipelineElements joined
    together as a data analysis or reduction pipeline
    
    joins should be a sequence of 3-tuples of the form 
    (elementindex,attribute,Pipelineobject) that are used to populate the 
    specified element's attribute by extracting the output of the 
    supplied pipeline
    
    dochecks validates the inputs as the correct type
    """
    def __init__(self,elements,joins=None,dochecks=True):
        from collections import deque
        
        if dochecks:
            for e in elements:
                if not isinstance(e,PipelineElement):
                    raise ValueError('object %s is not a PipelineElement'%e)
            if joins is not None:
                for j in joins:
                    elements[j[0]]
                    if not isinstance(j[1],basestring):
                        raise TypeError('joins %s does not have a proper attachment point'%j)
                    if not isinstance(joins[2],Pipeline):
                        raise ValeuError('joins %s does not have a Pipeline'%j)
                
        self.elements = list(elements)
        self.cycles = [0 for i in range(len(elements))]
        self.datadeques = [deque() for i in range(len(elements)+1)]
        if joins is None:
            self.joins = None
        else:
            self.joins = jdict = defaultdict(list)
            for j in joins:
                jdict[j[0]].append((j[1],j[2]))
                
    def feed(self,data):
        """
        feed an initial datum into the first stage of the pipeline.
        """
        self.datadeques[0].appendleft(data)
    
    def feedIter(self,dataiter):
        """
        feed initial data into the first stage of the pipeline.  The input
        will be treated as an iterator and each item will be added to
        the pipeline
        """
        for dat in dataiter:
            self.datadeques[0].appendleft(dat)
    
    
    def extract(self,extractall=False,autoprocess=False):
        """
        extract from the last stage of the pipeline.  
        
        If autoprocess is True, the stages will be processed 
        so that at least one item is in the final stage
        
        if extractall is true, a sequence of objects will be extracted
        from the final stage (combined with autoprocess, this will
        clear the pipeline)
        
        raises an IndexError if there is no data left to extract
        """              
        d = self.datadeques[-1]
          
        if autoprocess:
            self.processToStage(-1,True)
        
        if len(d) == 0:
            extractall = False
        
        if extractall:
            alld = list(reversed(d))
            d.clear()
            return alld
        else:
            return d.pop()
        
    def processStage(self,stagenum):
        st = self.elements[stagenum]
        
        #process the items to join into the stage if necessary
        if self.joins is not None:
            for attrname,joinpl in self.joins[stagenum]:
                try:
                    joinres = joinpl.extract(True,True)
                    setattr(st,attrname,joinres)
                except IndexError,e:
                    if e.message=='pop from an empty deque':
                        pass #empty pipeline
                    else:
                        raise
        
        data = self.datadeques[stagenum].pop()
        try:
            check_type(st._plintype,data) #does nothing if st._plintype is None
            newdata = st.plProcess(data,self,stagenum)
            if newdata is None:
                newdata = st.plInteract(data,self,stagenum)
                if newdata is None:
                    self.cycles[stagenum] += 1
                else:
                    self.datadeques[stagenum+1].appendleft(newdata)
                    self.cycles[stagenum] = 0
            else:
                self.datadeques[stagenum+1].appendleft(newdata)
                self.cycles[stagenum] = 0
        except TypeError,e:
            #if type-checking fails, let the data disppear
            raise TypeError('TypeError in stage %i(%s), removing invalid data '%(stagenum,e))
        except:
            self.datadeques[stagenum].append(data)
            raise
    
    def process(self,repeat=False):
        """
        process the pipeline to the end
        
        if repeat is True, the pipeline will by processed continually, 
        otherwise, only one step through the pipeline will be run
        """
        return self.processToStage(-1,repeat)
    
    def processSingle(self,input,processinglimit=10):
        """
        processes the input value through the pipeline, skipping over anything
        in the queue and repeating a processing stage until complete
        
        processinglimit is the maximum number of times to run the 
        processing of any given stage - if the processing is still 
        incomplete, a PipelineError will be raised.  If it is 0,
        no limit will be used (e.g. processing continues until 
        completion)
        
        return value is the result of the full pipeline (automatically 
        extracted)
        """
        self.datadeques[0].append(input)
        for stage in range(len(self.elements)):
            self._processStage(stage)
            while self.cycles[stage]:
                if processinglimit and self.cycles[stage] >= processinglimit:
                    raise PipelineError('hit processing limit at stage %i'%stage)
                self._processStage(stage)
            self.datadeques[stage+1].append(self.datadeques[stage].popleft())
                
        return self.datadeques[-1].pop()
        
    def processToStage(self,stagenum,repeat=False):
        """
        stage numbers are 0-indexed or can be negative to go from end
        
        if repeat is True, processing will continue until all earlier
        stages are empty.  If it is an integer, it will be taken as a 
        maximum number of times to attempt any given stage before
        a PipelineError is raised
        
        returns the number of times all of the stages were fed and processed
        """
        if stagenum >= 0:
            stages = range(stagenum+1)
        else:
            stages = range(len(self.datadeques)+stagenum)
        counts = []    
        for st in stages:
            dd = self.datadeques[st]
            count=0
            if repeat:
                while len(dd)>0:
                    self._processStage(st)
                    while self.cycles[st]:
                        if repeat is not True and self.cycles[st] >= repeat:
                            raise PipelineError('hit processing limit at stage %i'%st)
                        self._processStage(st)
                    count+=1   
            else:
                if len(dd)>0:
                    self._processStage(st)
                    if not self.cycles[st]:
                        count=1
            counts.append(count)
        return counts
    
    def clear(self,stages=None):
        """
        clears the inputs of the stage(s) requested (all of them if None)
        """
        if stages is None:
            stages = range(len(self.elements)+1)
        for st in stages:
            self.datadeques[st].clear()
    
class PipelineElement(object):
    """
    This class represents an element in a Pipeline.  Implementing classes
    must override the following method:
    
    * :meth:`plProcess`
        Process the data in whatever manner is appropriate.
    
    This method may be overridden: 
    * :meth:`plInteract`    
        Performs interactive stage. The :attr:`_plintype` attribute may be set
        (on classes or instances) as a processing indicator. See
        :func:`check_type` for valid types to use (it it is None, no checking is
        performed)
    """
    __metaclass__ = ABCMeta
    
    _plintype = None
    
    @abstractmethod
    def plProcess(self,data,pipeline,elemi):
        """
        this method should perform the data processing steps that 
        this element of the pipeline is supposed to do.  
        
        pipeline is the object that called the processing of this element,
        and elemi is the index of this stage of the pipeline
        
        this method should return None if processing was incomplete for any 
        reason - otherwise, the return value will be passed to the next 
        stage.
        
        Note that if this returns None and the data provided is not saved, it 
        will disappear - the next processing attempt will feed in the next
        piece of data.  If this is not the desired behavior, this method
        should do the following:
        pipeline.datadeques[elemi].append(data)
        """
        raise NotImplementedError
    
    def plInteract(self,data,pipeline,elemi):
        """
        called if :meth:`plProcess` returns None. Should return the data to be
        passed into the next stage (if None, :meth:`plProcess` will be called
        again)
        
        arguments are the same as :meth:`plProcess`
        """
        return None
        
    

class PipelineAccumulate(PipelineElement):
    """
    This pipeline element stores up everything that is fed into it until
    a condition is met, and then emits a list of objects with all objects
    previously stored.
    
    If accumsize is not None, if the number of objects stored ever 
    reaches accumsize, the list will be emitted
    
    if stoponobj is True, and the input matches stopobj
    (comparison uses "is" rather than "=="), the list will be emitted.  
    If stoponobj is True, accumsize determines when the sequence is emitted 
    
    onempty determines the behavior when inputs run out: if it is 'exception',
    an IndexError will be raised by the pipeline.  if 'block', nothing will be 
    emitted but the pipeline will continue to process.  if 'pass', a shortened
    sequence will be emitted (possibly empty if there are no inputs)
    """
    
    def __init__(self,accumsize=None,stopobj=None,stoponobj=True,onempty='exception'):
        self.accumsize = accumsize
        self.stopobj = stopobj
        self.stoponobj = stoponobj
        self.onempty = onempty
        
        self._accum = []
        
        
    def _getonempty(self):
        return self._onempty
    def _setonempty(self,val):
        if val not in ('exception','block','pass'):
            raise ValueError('invalid onempty "%s"'%val)
        self._onempty = val
    onempty = property(_getonempty,_setonempty)
    
    def plProcess(self,data,pipeline,elemi):
        self._accum.append(data)
        
        if self.stoponobj and data == self.stopobj:
            accumres = self._accum
            self._accum = []
            return accumres
        
        if self.accumsize is not None and len(self._accum) >= self.accumsize:
            accumres = self._accum
            self._accum = []
            return accumres
        
        if len(pipeline.datadeques[elemi]) == 0:
            if self._onempty == 'block':
                pipeline.datadeques[elemi].append(self._accum.pop())
            elif self._onempty == 'pass':
                accumres = self._accum
                self._accum = []
                return accumres
            elif self._onempty == 'exception':
                pass
            else:
                raise RuntimeError('impossible value for onempty')


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
    elif np.isscalar(method):
        res = method
    elif callable(method):
        res = method(arr)
    if method == 'median':
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
    Compute moments of the provided n-d array arr.
    
    The desired order of moment is given by the ms argument, which can
    be a scalar, or a sequence of length equal to the input's dimensions.
    If scalar, a sequence of moments (size equal to the number of input 
    dimensions)
    
    If axes are None, the output axes will be in 0-based pixels.  Otherwise,
    this should be a seqence of arrays where each element is of the same 
    length as the corresponding input array dimension.
    
    If bgmethod is not None, a background estimate will be subtracted before
    the moments are computed - see `estimate_background` method keywork 
    for useable values.
    
    if norm is True, the moment will be normalized - i.e.
    sum(x^m*arr)/sum(arr) instead of just sum(x^m*arr)
    
    if std is True, the output will be standardized (mth moment divided 
    by standard deviation to the mth power)
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
    Convinience function calling `moments`  to get just the first
    normalized moment (e.g. the centroid)
    """
    return moments(val,1,None if axes is None else (axes,),offset,True,False)

def sigma_clip(data,sig=3,iters=1,center='median',maout=False):
    """
    This performs the sigma clipping algorithm - i.e. the data will
    be iterated over `iters` times, each time rejecting points 
    that are more than `sig` standard deviations discrepant 
    
    center is the estimation technique to compute the  assumed 
    center of the clipping - see `estimate_background` for options
    
    if maout is True, returns a numpy.ma.Maskedarray with the filtered points
    masked.  Otherwise, returns a tuple (filtereddata,mask) 
    (mask is False for clipped points, shaped as original) 
    """
    
    data = np.array(data,copy=False)
    oldshape = data.shape
    data = data.ravel()
    
    mask = np.ones(data.size,bool)
    for i in range(iters):
        do = data-estimate_background(data[mask],center)
        mask = do*do <= np.var(data[mask])*sig
        
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
    
    
    Example:
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

def lin_to_log_rescale(val,lower=1,upper=1000,base=10):
    """
    maps linear values onto the range [lower,upper] and then applies the
    requested base of logarithm
    """
    if lower > upper:
        raise ValueError('lower must be less than upper')
    if lower <= 0:
        raise ValueError('lower must be greater than 0')
    
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
    Returns a boolean mask :class:`numpy.ndarray` for the point where the
    input (interpreted as a 1D array) crosses or has complted crossing (if
    between) from below (or above) a threshold.
    
    If `belowtoabove` is True, the returned masks is for where the input
    transitions from below to above.  Otherwise, from above to below.
    
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

#<--------------------------Robust statistics---------------------------------->    
def interquartile_range(values,scaletonormal=False):
    """
    Computes the interquartile range for the provided sequence of values, a
    more robust estimator than the variance.
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
    
def median_absolute_deviation(values,scaletonormal=False):
    """
    computes the median_absolute_deviation for the provided sequence of values, 
    a more robust estimator than the variance.
    """
    from scipy.special import erfinv
    
    x = np.array(values,copy=False).ravel()
    res = np.median(np.abs(x-np.median(x)))
    
    if scaletonormal:
        nrm = (2**0.5*erfinv(.5))
        return res/nrm
    else:
        return res
    
def biweight_midvariance(values,influencescale=9):
    """
    Computes the biweight midvariance of a sequence of data points, a robust 
    statistic of scale.
    
    `influence` sets the number of MAD units away at which a data point has no
    weight
    
    returns bmv,median
    """    
    x = np.array(values,copy=False).ravel()
    
    Q = np.median(x)
    MAD = median_absolute_deviation(x)
       
    ui = (x-Q)/(influencescale*MAD)
    uisq=ui**2
    
    I = uisq <= 1
    
    numer = np.sum(I*(x-Q)**2*(1-uisq)**4)
    denom = np.sum(I*(1-uisq)*(1-5*uisq))**2
    
    
    return x.size*numer/denom,Q

del ABCMeta,abstractmethod,abstractproperty #clean up namespace