#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains various astropysics or python-related utility
functions that have no obvious place in some other module
"""

from __future__ import division,with_statement
import numpy as np

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    
    
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


#<---------------------------data pipelining----------------------------------->

class Pipeline(object):
    """
    This class represents a pipeline, composed of PipelineElements joined
    together as a data analysis or reduction pipeline
    
    dochecks validates the components as PipelineElements before accepting them
    into the pipeline
    """
    def __init__(self,components,joins=None,dochecks=True):
        from collections import deque
        
        if joins is not None:
            raise NotImplementedError('non-linear pipelines not yet ready')
        
        if dochecks:
            for c in components:
                if not isinstance(c,PipelineElement):
                    raise ValueError('object %s is not a PipelineElement'%str(c))
                
        self.components = list(components)
        self.processcycles = [0 for i in range(len(components))]
        self.datadeques = [deque() for i in range(len(components)+1)]
    
    def feed(self,data,iterate=False):
        """
        feed initial data into the first stage of the pipeline
        
        if iterate is True, the input will be taken as an iterator to feed
        each element into the first stage
        """
        if not iterate:
            data = [data]
            
        for dat in data:
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
        st = self.components[stagenum]
        lastst = self.components[:(stagenum-1)] if stagenum != 0 else None
        
        data = self.datadeques[stagenum].pop()
        try:
            check_type(st._plintype,data) #does nothing if st._plintype is None
            newdata = st._plProcess(data,lastst,self)
            if newdata is None:
                newdata = st._plInteract(data,lastst,self)
                if newdata is None:
                    self.processcycles[stagenum] += 1
                else:
                    self.datadeques[stagenum+1].appendleft(newdata)
                    self.processcycles[stagenum] = 0
            else:
                self.datadeques[stagenum+1].appendleft(newdata)
                self.processcycles[stagenum] = 0
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
        for stage in range(len(self.components)):
            self._processStage(stage)
            while self.processcycles[stage]:
                if processinglimit and self.processcycles[stage] >= processinglimit:
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
                    while self.processcycles[st]:
                        if repeat is not True and self.processcycles[st] >= repeat:
                            raise PipelineError('hit processing limit at stage %i'%stage)
                        self._processStage(st)
                    count+=1   
            else:
                if len(dd)>0:
                    self._processStage(st)
                    if not self.processcycles[st]:
                        count=1
            counts.append(count)
        return counts
    
    def clear(self,stages=None):
        """
        clears the inputs of the stage(s) requested (all of them if None)
        """
        if stages is None:
            stages = range(len(self.components)+1)
        for st in stages:
            self.datadeques[st].clear()
    
class PipelineElement(object):
    """
    This class represents an element in a Pipeline.  Implementing classes
    must override the following method:
    *_plProcess(self,data,last,pipeline): process the data
    
    This method may be overridden: 
    *_plinteract(self,data,last,pipeline): perform interactive stage - 
        typical intended calling order is _plProcess,_plInteract
    
    this attribute may be set (on classes or instances) as a processing 
    indicator:
    *_plintype: the type expected for data fed into this PipelineElement - 
        see `check_type` for valid types to use (None means no checking)
    """
    __metaclass__ = ABCMeta
    
    _plintype = None
    
    @abstractmethod
    def _plProcess(self,data,srcstages,pipeline):
        """
        this method should perform the data processing steps that 
        this element of the pipeline is supposed to do.  
        
        srcstages is a sequence of all the previous stages in the pipeline
        
        pipeline is the object that called the processing of this element
        
        this method should return None if processing was incomplete for any 
        reason - otherwise, the return value will be passed to the next 
        stage
        """
        raise NotImplementedError
    
    def _plInteract(self,data,srcstage,pipeline):
        """
        called if _plProcess retruns None and plinteract is True.  Should 
        return the data to be passed into the next stage (if None, _plProcess
        will be called again)
        """
        return None
        
    
class PipelineError(Exception):
    """
    This exception is used to indicate problems with the pipeline.  It
    should NOT be used to indicate data matching problems within 
    PipelineElements -- ValueError or TypeError should be used for 
    that purpose
    """
    pass



#<--------------------Analysis/simple numerical functions---------------------->
def centroid(val,axes=None,offset='median'):
    """
    computes the centroid of an n-dimensional array
    
    axes can either be a list of axis values for the array (must be
    same length as the number of dimensions) or None to report just the 
    pixel value
    
    offset is an amount to subtract from the input array, can be:
    *'median'
    *'mean'
    *None
    """
    if offset:
        if offset == 'median':
            val = val - np.median(val)
        elif offset == 'mean':
            val = val - np.mean(val)
        else:
            raise ValueError('unrecognized offset type')
    
    shp = val.shape
    if axes is None:
        axes=[np.arange(s) for s in shp]
    if len(shp) != len(axes):
        raise ValueError("input dimensions don't match axes")
    
    cens = []
    for i,(s,a) in enumerate(zip(shp,axes)):
        if s != len(a):
            raise ValueError("axis #%i doesn't match array size"%i)
        rng = range(len(s))
        del rng[i]
        vsum = val
        for j in reversed(rng):
            vsum = val.sum(axis = j)
        cens.append(np.sum(x*vsum)/np.sum(vsum))
    
    return tuple(cens)

del ABCMeta,abstractmethod,abstractproperty #clean up namespace