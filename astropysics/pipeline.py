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

=============================================================
pipeline -- classes for data reduction and analysis pipelines
=============================================================

The :mod:`pipeline` module contains classes and utilities for constructing data
pipelines -- linear constructs of operations that process input data, passing it
through all pipeline stages.

Pipelines are represented by the :class:`Pipeline` class, which is composed of a
sequence of :class:`PipelineElement` objects representing the processing stages.
Astropysics includes some pipeline elements built-in, but
:class:`PipelineElement` can also easily be subclassed to provide cusotmized
pipeline functionality.  

Once the pipeline is constructed, data can be queued into the first stage with
the :meth:`Pipeline.feed` method (or :meth:`Pipeline.feedMany`), where the type
of the data expected is determined by the :class:`PipelineElement` objects.

Special objects (subclasses of :class:`PipelineMessage`) can also be fed into
pipelines to indicate to various pipeline elements to change how they behave.
This allows a pipeline to adapt based on the data. e.g. If two different
exposure times are being processed in a image reduction pipeline, after the
images from the first set of exposure times are processed, a
:class:`AccumulateMessage` can be added to provide a new set of darks before the
second set of images are fed into the pipeline.

.. todo:: replace above with example or add a stand-alone example

.. todo:: add list of pipelineelements to Sphinx docs with types






Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.pipeline
   :parts: 1

Module API
----------

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
    This class represents a pipeline, composed of a linear sequence of
    :class:`PipelineElements <PipelineElement>` joined together.
    
    Note that where stage numbers are mentioned in the documentation and inputs,
    they are 0-indexed (e.g. the 1st stage is stage 0).
    
    """
    
    
    def __init__(self,elements=None):
        """
        :param elements: Initial pipeline elements for this pipeline.
        :type elements: 
            Sequence of :class:`PipelineElements <PipelineElement>` or None
            
        :except TypeError: If the input is not a :class:`PipelineElement`
        """
        from collections import deque
        
        if elements is None:
            elements = tuple()
        
        for e in elements:
            if not isinstance(e,PipelineElement):
                raise TypeError('object %s is not a PipelineElement'%e)
                
        self._elements = list(elements)
        self._cycles = [0 for i in range(len(elements))]
        self._datadeques = [deque() for i in range(len(elements)+1)]
        
    
    @property
    def elements(self):
        """
        A tuple containing the pipeline elements.
        """
        return tuple(self._elements)
    
    def addElement(self,stagenum,element):
        """
        Insert a :class:`PipelineElement` to the pipeline at the specified 
        stage number. Later stages are pushed forward along with their current 
        data.
        
        :param int stagenum: The stage number for the new element.
        :param PipelineElement element: The :class:`PipelineElement` object to insert.
        
        :except TypeError: If the input is not a :class:`PipelineElement`
        """
        from collections import deque
        
        if not isinstance(element,PipelineElement):
            raise TypeError('object %s is not a PipelineElement'%element)
        
        self._elements.insert(stagenum,element)
        self._cycles.insert(stagenum,0)
        self._datadeques.insert(stagenum,deque())
    
    def removeElement(self,stagenum,keepdata=False):
        """
        Remove the requested stage number. All later stages will move one stage
        left to fill the gap.
        
        :param int stagenum: The stage number to remove.
        :param bool keepdata: 
            If True, any data waiting to be processed in the removed stage is
            added to the next stage. Otherwise, the data is lost.
        
        
        """
        elem = self._elements.pop(stagenum)
        data = self._datadeques.pop(stagenum)
        if keepdata:
            self._datadeques[stagenum].extend(data)
            self._datadeques[stagenum].rotate(len(data))
        del self._cycles[stagenum]
        
        return elem
                
    def feed(self,data):
        """
        Feed an initial datum into the first stage of the pipeline.
        """
        self._datadeques[0].appendleft(data)
    
    def feedMany(self,dataeq):
        """
        Feed initial data into the first stage of the pipel
        """
        for dat in dataeq:
            self._datadeques[0].appendleft(dat)
    
    
    def extract(self,extractall=False,autoprocess=False):
        """
        Extract from the last stage of the pipeline.
        
        :param extractall:
            If True, a sequence of objects will be extracted from the final
            stage (combined with `autoprocess`, this will clear the pipeline).
        :type extractall: bool
        :param autoprocess:
            If True, the stages will be repeatedly processed so that at least
            one item is in the final stage. (combined with `extractall`, this
            will clear the pipeline).
        :type autoprocess: bool
        
        :raises IndexError: If there is no data left to extract
        """              
        d = self._datadeques[-1]
          
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
        """
        Process a particular stage once, possibly processing earlier stages
        if there is nothing in the request stage's input.
        
        :param stagenum: Stage number to process
        :type stagenum: int
        :returns: 
            False if the stage didn't complete, True if it did, or the string
            'message' if it delivered a message instead of data.
        """
        from .utils import check_type
        
        st = self._elements[stagenum]
        
        data = self._datadeques[stagenum].pop()
        if isinstance(data,PipelineMessage):
            msg = data
            istarg = msg.isTarget(st)
            retcode = ''
            
            if istarg:
                msg(st)
                retcode += 'message_delivered'
                
            if not istarg or istarg=='continue':
                if stagenum == len(self._elements)-1:
                    if not msg.nprocessed:
                        raise PipelineError('message %s was never delivered to any element'%msg)
                else:
                    self._datadeques[stagenum+1].appendleft(msg)
                retcode += ('' if retcode=='' else ',')+'message_passed'
            return retcode
        else:
            try:
                check_type(st._plintype,data) #no-op if st._plintype is None
            except TypeError,e:
                #if type-checking fails, let the data disppear
                raise TypeError('TypeError in stage %i(%s), removing invalid data '%(stagenum,e))
            try:
                newdata = st.plProcess(data,self,stagenum)
                if newdata is None:
                    newdata = st.plInteract(data,self,stagenum)
                    if newdata is None:
                        self._cycles[stagenum] += 1
                        return False
                    else:
                        self._datadeques[stagenum+1].appendleft(newdata)
                        self._cycles[stagenum] = 0
                        return True
                else:
                    self._datadeques[stagenum+1].appendleft(newdata)
                    self._cycles[stagenum] = 0
                    return True
            
            except:
                self._datadeques[stagenum].append(data)
                raise
    
    def process(self,repeat=False):
        """
        Process the pipeline to the end.
        
        :param repeat:
            
            If True, the pipeline will by processed continually. If False, only
            one step through the pipeline will be run. If it is an integer the
            pipeline will by processed continually, it will be taken as a
            maximum number of times to attempt any given stage before a
            PipelineError is raised.
            
        :type repeat: bool or int
        
        :except PipelineError: If a stage is repeated `repeat` or more times.
        
        :returns:
            A list with the return value of :func:`processStage` if `repeat` is
            False, or a list of lists if `repeat` is True.
        """
        return self.processToStage(-1,repeat)
    
    def processSingle(self,input,processinglimit=10,removefailed=True):
        """
        Processes the input value through the pipeline, skipping over anything
        in the queue and repeating a processing stage until complete.
        
        :param input: The input value to be feed into the 1st stage
        :param processinglimit: 
            The maximum number of times to attempt to process any stage. If 0,
            processing continues indefinitely (or until completion).
        :type processinglimit: int
        :param removefailed: 
            If True and a stage reaches the `processinglimit`, the data object
            will be removed from the pipeline.  Otherwise, the object will be
            left in place where it got stuck.
        :type removefailed: bool
        
        :except PipelineError: 
            If `processinglimit` is reached and the object has still not been
            processed.  The data will be left in the front of the failed stage's
            queue unless `removefailed` is True.  
        
        :returns:
            The result of the full pipeline (it will *not* be left in the final
            stage).
        """
        self._datadeques[0].append(input)
        for stage in range(len(self._elements)):
            datapreprocess = self._datadeques[stage][-1]
            self.processStage(stage)
            while self._cycles[stage]:
                if processinglimit and self._cycles[stage] >= processinglimit:
                    if not removefailed:
                        self._datadeques[stage].append(datapreprocess)
                    raise PipelineError('hit processing limit at stage %i'%stage)
                self.processStage(stage)
            #TODO:add checks?
                
        return self._datadeques[-1].pop()
        
    def processToStage(self,stagenum,repeat=False):
        """
        Processes the pipeline from the beginning up to the requested stage.
        
        :param stagenum: The stage to process to.
        :type stagenum: int
        :param repeat:
            If True, processing will continue until all earlier stages are
            empty. If it is an integer, it will be taken as a maximum number of
            times to attempt any given stage before a PipelineError is raised
        :type repeat: bool or int
        
        :returns: 
            A list with the return value of :func:`processStage` if `repeat` is
            False, or a list of lists if `repeat` is True.
        
        :except PipelineError: 
            If an element has not completed after `repeat` processing cycles
            have been run.
        """
        if stagenum >= 0:
            stages = range(stagenum+1)
        else:
            stages = range(len(self._datadeques)+stagenum)
        results = []    
        for st in stages:
            dd = self._datadeques[st]
            
            if repeat:
                result = []
                while len(dd)>0:
                    self.processStage(st)
                    while self._cycles[st]:
                        if repeat is not True and self._cycles[st] >= repeat:
                            raise PipelineError('hit processing limit at stage %i'%st)
                        result.append(self.processStage(st))
            elif len(dd)>0:
                result = self.processStage(st)
            else:
                result = False
            results.append(result)
        return results
    
    def clear(self,stages=None):
        """
        Clears the inputs of the stage(s) requested.
        
        :param stages: 
            The stage(s) to be cleared.  If None, all will be cleared.
        :type stages: integer, sequence of integers or None
        """
        if isinstance(stages,int):
            stages=(stages,)
        
        if stages is None:
            stages = range(len(self._elements)+1)
        for st in stages:
            self._datadeques[st].clear()
    
class PipelineElement(object):
    """
    This class represents an element in a Pipeline. The :attr:`_plintype`
    attribute may be set (on classes or instances) as a processing indicator.
    See :func:`check_type` for valid types to use (it it is None, no checking is
    performed)
    
    
    **Subclassing**
    
    Implementing classes must override the following method:
    
    * :meth:`plProcess`
    
    This method should be overridden if this stage of the pipeline involves an
    interactive step: 
    
    * :meth:`plInteract`    
        
        
    """
    __metaclass__ = ABCMeta
    
    _plintype = None
    
    @abstractmethod
    def plProcess(self,data,pipeline,elemi):
        """
        This method performs whatever data processing steps that this element of
        the pipeline is supposed to do.
        
        :param data: 
            The data massed in by the previous element of the pipeline. Type
            and interpretation is left to this method to test.
        :param pipeline: 
            The :class:`Pipeline` object that called this element.
        :type pipeline: :class:`Pipeline`
        :param int elemi: The index of this stage of the pipeline.
        
        :returns:
            None if processing was incomplete for any reason - otherwise, the
            return value will be passed to the next stage.
        
        .. note:: 
        
            If this returns None and the data provided is not saved, it will
            disappear - the next processing attempt will feed in the next piece
            of data. If this is not the desired behavior, this method should 
            call ``resaveData(data,pipeline,elemi)``
                
            
        """
        raise NotImplementedError
    
    def resaveData(self,data,pipeline,elemi):
        """
        This returns the data so that the next time this element is processed,
        it will receive the same data. This is intended to be used inside
        :meth:`plProcess` if it returns None and the data should be reprocessed
        the next time the pipeline stage is run.
        """
        pipeline._datadeques[elemi].append(data)
    
    def plInteract(self,data,pipeline,elemi):
        """
        This method should be implemented if an interactive step is to be
        performed by this element. In that case, :meth:`plProcess` should return
        None, and this will be called after plProcess completes.
        
        :param data: 
            The data massed in by the previous element of the pipeline. Type
            and interpretation is left to this method to test.
        :param pipeline: 
            The :class:`Pipeline` object that called this element.
        :type pipeline: :class:`Pipeline`
        :param int elemi: The index of this stage of the pipeline.
        
        :returns:
            The data to be passed into the next stage or None if the interaction
            was not completed. In this case, the input will be run through
            :meth:`plProcess` again
        
        """
        return None
    
class PipelineMessage(object):
    """
    This class represents a message that is passed through the pipeline as
    though it were data, but is never actually processed by the elements -
    instead, it performs an action when it reaches a particular target(s).
    
    **Subclassing**
    
    * Subclasses *must* implement :meth:`deliverMessage` (see
      :meth:`deliverMessage` docs for details)
    * Subclasses should call :meth:`PipelineMessage.__init__` if overridding 
      :meth:`.__init__`.
    * Subclasses can override :meth:`isTarget` to change how the message
      determines what its target is. This can also be set to deliver the message
      to multiple elements.
    
    """
    __metaclass__ = ABCMeta
    
    def __init__(self,target):
        """
        :param target: 
            Either a specific object that this :class:`PipelineMessage` should
            be delivered to, or a class. If it is a class, the message will be
            delivered to all of that class in the pipeline, otherwise it will
            only reach the actual object specified as the target.
            
        """
        self.target = target
        self.nprocessed = 0
    
    def __call__(self,elem):
        """
        Perform whatever action this message is supposed to perform on the
        supplied pipeline element.  
        """
        self.nprocessed += 1
        self.deliverMessage(elem)
        
    def isTarget(self,elem):
        """
        Identifies whether or not the supplied element is an intended target of
        this message.
        
        :param elem: The element to be tested.
        :returns: 
            bool indicating whether or not `elem` is the target. It can also be
            the string 'continue' if the message should be delivered to this
            target but also passed further down the pipeline.
        """
        from inspect import isclass
        
        if isclass(self.target):
            return 'continue' if elem.__class__ == self.target else False
        else:
            return elem is self.target
            
        
        
    @abstractmethod
    def deliverMessage(self,elem):
        """
        This method must be overridden in subclasses, defining the action to be
        performed on the pipeline element this message is targeted at. 
        
        :param elem: 
            The element this message is to be delivered to. It is guaranteed to
            be an object :meth:`isTarget` method accept.
        :type elem: :class:`PipelineElement`
        
        """
        raise NotImplementedError
    
class SetAttributeMessage(PipelineMessage):
    """
    This object is a :class:`PipelineMessage` that sets a specified set of
    attributes when it reaches its target.
    """
    
    def __init__(self,target,**attrs):
        """
        :param target: The target object to receive this message.
        
        Other keyword arguments are taken as the attributes to be set on the
        target.  Thus, ::
            
            pipeline = Pipeline([apipelineelement])
            msg = SetAttributeMessage(apipelineelement,attr1=val1,attr2=val2)
            pipeline.feed(msg)
            pipeline.process()
            
        is equivalent to ::
            
            pipeline = Pipeline([apipelineelement])
            apipelineelement.attr1 = val1
            apipelineelement.attr2 = val2
        
        A dictionary can also be passed in as kwargs::
            
            pipeline = Pipeline([apipelineelement])
            attrs = {'attr1':val1,'attr2':val2}
            msg = SetAttributeMessage(apipelineelement,**attrs)
            pipeline.feed(msg)
            pipeline.process()
            
        
        """
        PipelineMessage.__init__(self,target)
        self.attrstoset = attrs
        
    def deliverMessage(self,elem):
        for name,value in self.attrstoset.iteritems():
            setattr(elem,name,value)
            
class CallMethodMessage(PipelineMessage):
    """
    This object is a :class:`PipelineMessage` that calls a method on the target
    pipeline element.
    
    The :attr:`retval` attribute stores the return value for the method call as 
    a list (in the same order that the calls occur, if the message is passed to 
    multiple elements).
    
    """
    
    def __init__(self,target,methodname,*args,**kwargs):
        """
        :param target: The target pipeline element to receive this message.
        :param string methodname: The name of the method that should be called.
        
        Further arguments and keyword arguments will be passed in as the 
        arguments and keyword for the method call. Thus, ::
            
            pipeline = Pipeline([apipelineelement])
            msg = CallMethodMessage(apipelineelement,'meth',arg1,kw2=val2)
            pipeline.feed(msg)
            pipeline.process()
            
        is equivalent to ::
            
            pipeline = Pipeline([apipelineelement])
            apipelineelement.meth(arg1,kw2=val2)
        
        """
        PipelineMessage.__init__(self,target)
        self.methodname = methodname
        self.methargs = args
        self.methkwargs = kwargs
        self.retval = []
        
    def deliverMessage(self,elem):
        self.retval.append(getattr(elem,self.methodname)(*args,**kwargs))
      
class _Accumulator(object):
    """
    A helper class left by an :class:`AccumulateMessage` at its target to do the 
    actual accumulating
    """
    
    def __init__(self,elem,naccum,setattrname,filterfunc):
        self.naccum = naccum
        self.setattrname = setattrname
        self.filterfunc = filterfunc
        
        self.accum = []
    
        self.realPlProcess = elem.plProcess
        self.realPlInteract = elem.plInteract
        
        #inject this in place of the real pipeline element - the finish method
        #will put the right methods back in.
        elem.plProcess = self.plProcess
        elem.plInteract = self.plInteract
        elem.accumulator = self
        self.elem = elem
    
    def plProcess(self,data,pipeline,elemi):
        assert len(self.accum) < self.naccum,'Accumulator plProcess is being called after it should have finished'
        self.accum.append(data)
        return None #unneccessary, but explicit is better than implicit...
        
    def plInteract(self,data,pipeline,elemi):
        #Finishing is done here so the pipeline doesn't try to do the real interaction step
        if len(self.accum) == self.naccum:
            self.finish(data,pipeline,elemi)
        return None
    
    def finish(self,data,pipeline,elemi):
        try:
            accum = self.accum
            if self.filterfunc is not None:
                accum = self.filterfunc(accum)
            if self.setattrname is None:
                self.elem.resaveData(accum,pipeline,elemi)
            else:
                setattr(self.elem,self.setattrname,accum)
        finally:
            self.elem.plProcess = self.realPlProcess 
            self.elem.plInteract = self.realPlInteract
            del self.elem.accumulator 
            
        
        
        
class AccumulateMessage(PipelineMessage):
    """
    This object is a :class:`PipelineMessage` that tells the target
    :class:`PipelineElement` to accumulate the next `naccum` inputs from the
    earlier stage as a list, and then either pass this list in as data to the
    next stage, or set it to an attribute on the element.
    
    """
    def __init__(self,target,naccum=1,setattrname=None,filterfunc=None):
        """
        :param target: The target pipeline element to receive this message.
        :param int naccum: 
            The number of inputs to accumulate before continuing normal
            operation of the pipeline element.
        :param setattrname: 
            The name of the attribute that the result of the accumulation should
            be set to. Alternatively, it can be None, in which case the
            accumulated list will be passed in as the next set of input data.
        :param filterfunc: 
            A callable that is called when as ``filterfunc(accum)`` where
            `accum` is the accumulated list of data, and should return the
            object to be set to the attribute `setattrname`. If None, this step
            is not performed (e.g. a list with the accumulated data is the
            result).
        """
        
        PipelineMessage.__init__(self,target)
        self.naccum = naccum
        self.setattrname = setattrname
        self.filterfunc = filterfunc
        
    def deliverMessage(self,elem):
        _Accumulator(elem,self.naccum,self.setattrname,self.filterfunc)
        
del ABCMeta,abstractmethod,abstractproperty #clean up namespace
