#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains various astropysics or python-related utility
functions that have no obvious place in some other module, as well
as a variety of general-purpose algorithms that are not particular
to a module.
"""

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
            print 'at',key,val
            oldval = self[key] if key in self else None
            print 'ov',oldval
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
        
    def _processStage(self,stagenum):
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
            newdata = st._plProcess(data,self,stagenum)
            if newdata is None:
                newdata = st._plInteract(data,self,stagenum)
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
    *_plProcess(self,data,last,pipeline): process the data
    
    This method may be overridden: 
    *_plinteract(self,data,last,pipeline): perform interactive stage     
    this attribute may be set (on classes or instances) as a processing 
    indicator:
    *_plintype: the type expected for data fed into this PipelineElement - 
        see `check_type` for valid types to use (None means no checking)
    """
    __metaclass__ = ABCMeta
    
    _plintype = None
    
    @abstractmethod
    def _plProcess(self,data,pipeline,elemi):
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
    
    def _plInteract(self,data,pipeline,elemi):
        """
        called if _plProcess returns None.  Should return the data to be 
        passed into the next stage (if None, _plProcess will be called again)
        
        arguments are the same as _plProcess
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
    
    def _plProcess(self,data,pipeline,elemi):
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
def centroid(val,axes=None,offset=None):
    """
    computes the centroid of an n-dimensional array
    
    axes can either be a list of axis values for the array (must be
    same length as the number of dimensions) or None to report just the 
    pixel value
    
    offset is an amount to subtract from the input array, can be:
    
    * 'median'
    * 'mean'
    * '32' : 3*median - 2*mean
    * '2515' : 2.5*median - 1.5*mean
    * '21' : 2*median - 1*mean
    * None
    """
    if offset:
        if offset == 'median':
            val = val - np.median(val)
        elif offset == 'mean':
            val = val - np.mean(val)
        elif offset == '32':
            val = val - (3*np.median(val) - 2*np.mean(val))
        elif offset == '2515':
            val = val - (2.5*np.median(val) - 1.5*np.mean(val))
        elif offset == '21':
            val = val - (2*np.median(val) - np.mean(val))
        else:
            raise ValueError('unrecognized offset type')
    
    shp = val.shape
    if axes is None:
        axes=[np.arange(s) for s in shp]
    if len(shp) != len(axes):
        raise ValueError("input dimensions don't match axes")
    
    cens = []
    for i,(s,x) in enumerate(zip(shp,axes)):
        if s != len(x):
            raise ValueError("axis #%i doesn't match array size"%i)
        rng = range(len(shp))
        del rng[i]
        vsum = val
        for j in reversed(rng):
            vsum = val.sum(axis = j)
        cens.append(np.sum(x*vsum)/np.sum(vsum))
    
    return tuple(cens)

def sigma_clip(data,sig=3,iters=1,center='median',maout=False):
    """
    This performs the sigma clipping algorithm - i.e. the data will
    be iterated over `iters` times, each time rejecting points 
    that are more than `sig` standard deviations discrepant 
    
    center can be any of:
    * 'median'
    * 'mean'
    * '32' : 3*median - 2*mean
    * '2515' : 2.5*median - 1.5*mean
    * '21' : 2*median - 1*mean
    * a callable that takes a 1D array and outputs a scalar
    
    if maout is True, returns a numpy.ma.Maskedarray with the filtered points
    masked.  Otherwise, returns a tuple (filtereddata,mask) 
    (mask is False for clipped points, shaped as original) 
    """
    if center == 'median':
        center = np.median
    elif center =='mean':
        center = np.mean
    elif callable(center):
        pass
    elif center == '32':
        center = lambda x:3*np.median(x) - 2*np.mean(x)
    elif center =='2515':
        center = lambda x:2.5*np.median(x) - 1.5*np.mean(x)
    elif enter == '21':
        center = lambda x:2*np.median(x) - np.mean(x)
    else:
        raise ValueError('unrecognized center value')
    
    data = np.array(data,copy=False)
    oldshape = data.shape
    data = data.ravel()
    
    mask = np.ones(data.size,bool)
    for i in iter:
        do = data-center(data[mask])
        mask = do*do <= np.var(data[mask])*sig
        
    if maout:
        return np.ma.MaskedArray(data,~mask,copy=False)
    else:
        return data[mask],mask.reshape(oldshape)

del ABCMeta,abstractmethod,abstractproperty #clean up namespace