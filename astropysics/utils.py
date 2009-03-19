#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains various astropysics or python-related utility
functions that have no obvious place in some other module, including
pipeline interfaces
"""

from __future__ import division
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
    """
    def __init__(self,components):
        for c in components:
            if not isinstance(c,PipelineElement):
                raise ValueError('object %s is not a PipelineElement'%str(c))
        raise NotImplementedError
    
    def feed(self,data):
        raise NotImplementedError
    
    def extract(self):
        raise NotImplementedError
    
    def setInterrupt(self,val):
        raise NotImplementedError
    
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
        data
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
    
del ABCMeta,abstractmethod,abstractproperty #clean up namespace