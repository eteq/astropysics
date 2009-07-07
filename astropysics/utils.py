#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains various astropysics or python-related utility
functions that have no obvious place in some other module
"""

from __future__ import division,,with_statement
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

#<---------------------------data pipelining----------------------------------->

class Pipeline(object):
    """
    This class represents a pipeline, composed of PipelineElements joined
    together as a data analysis or reduction pipeline
    
    dochecks validates the components as PipelineElements before accepting them
    into the pipeline
    
    TODO: conditional interrupts
    """
    def __init__(self,components,joins=None,dochecks=True):
        if joins is not None:
            raise NotImplementedError('non-linear pipelines not yet ready')
        
        if dochecks:
            for c in components:
                if not isinstance(c,PipelineElement):
                    raise ValueError('object %s is not a PipelineElement'%str(c))
                
        self.components = list(components)
    
    def feed(self,data,iterate=False):
        """
        feed initial data into the first stage of the pipeline
        
        if iterate is True, the input will be taken as an iterator to feed
        each element into the first stage
        """
        if not iterate:
            data = [data]
            
        for dat in data:
            self.components[0]._plFeed(dat,None)
    
    def extract(self,autoprocess=True,extractall=False):
        """
        extract from the last stage of the pipeline.  
        
        If autoprocess is True, the stages will be processed 
        so that at least one item is in the final stage
        
        if extractall is True, this will continue until None is returned
        """                
        
        if extractall:
            results = []
            if autoprocess:
                self.processToStage(-1,True)
            results.append(self.components[-1]._plExtract())
            while results[-1] is not None:
                results.append(self.components[-1]._plExtract())
            return results
        
        else:
            extracted = self.components[-1]._plExtract()
            if autoprocess and extracted is None:
                self.processToStage(-1,False)
                extracted = self.components[-1]._plExtract()
            return extracted
    
    def processToStage(stagenum,fullyprocess=False):
        """
        stage numbers are 0-indexed or can be negative to go from end
        
        if fullyprocess is True, processing continue until all earlier
        stages are empty
        
        returns the number of times the stage was fed and processed
        """
        prestages = self.components[:stagenum]
        finalstage = self.components[stagenum]
        
        count = 0
        fed = True
        while fed:
            tofeed = lastst = fed = None
            for st in prestages:
                if tofeed is not None:
                    st._plFeed(tofeed,lastst)
                st._plProcess()
                tofeed = st._plextract()
                if fullyprocess and tofeed is not None:
                    fed = True
                lastst = st
                    
            if tofeed is not None:
                count+=1
                finalstage._plFeed(tofeed,lastst)
                finalstage._plProcess()
        return count
    
    def clear(self,stages=None):
        """
        clears the stage(s) requested (all of them if None)
        """
        from operator import isSequenceType
        if stages is None:
            stages=range(len(self.components))
        if not isSequenceType(stages):
            stages = [stages]
        for st in stages:
            self.components[stages]._plClear()
    
class PipelineElement(object):
    """
    This class represents an element in a Pipeline.  Implementing classes
    must override the following methods:
    *_plFeed(self,data,src): load data into the pipeline element.  Minimal 
    processing should occur here (just whatever is necessary to set
    parameters)
    *_plProcess(self): process the data 
    *_plExtract(self): return the processed data 
    optionally, the following methods can be overriden:
    *_plClear(self): clear any stored data associated with this element
    """
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def _plFeed(self,data,src):
        """
        this method loads data from an earlier pipeline stage into this object
        
        data is the data from the earlier stage (interpretation is up to the
        object), while src is a refernce to the object that is feeding in the 
        data (or None if it is the first stage)
        """
        raise NotImplementedError
    
    @abstractmethod
    def _plProcess(self):
        """
        this method should perform the data processing steps that 
        this element of the pipeline is supposed to do.  
        
        It should return None when processing is completed, otherwise some 
        sort of informational object that is up to the Pipeline to interpret.
        The pipeline will continue processing until None is returned, although 
        there are pipeline mechanisms to interrupt at various points (see
        Pipeline.setInterrupt)
        """
        raise NotImplementedError
    
    @abstractmethod
    def _plExtract(self):
        """
        This method should return one piece of data from this stage of the pipeline 
        for use in the next stage.  Generally, information about the data should be 
        completely encapsulated in the object that is returned from this 
        method.  If there is no data left to be extracted, this should return
        None
        """
        raise NotImplementedError
        return data
    
    def _plClear(self):
        """
        This method should clear any stored data in the pipeline. by default 
        it processes and then extracts data until the pipeline returns None,
        but overriding is recommended if processing will be expensive per-unit
        """
        while _plProcess() is not None:
            pass
        while _plExtract() is not None:
            pass
    
class PipelineError(Exception):
    """
    This exception is used to indicate problems with the pipeline.  It
    should NOT be used to indicate data matching problems within 
    PipelineElements --ValueError or TypeError should be used for 
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