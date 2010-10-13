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
Astropysics includes a number of elements builtin, but :class:`PipelineElement`
can also easily be subclassed to provide new pipeline stages.

.. todo:: add list of pipelineelements to Sphinx docs

.. todo:: examples/tutorials




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
    This class represents a pipeline, composed of :class:`PipelineElements
    <PipelineElement>` joined together as a data analysis or reduction pipeline.
    
    Note that whenever stage numbers are used, they are 0-indexed 
    (e.g. the 1st stage is stage 0)
    
    """
    def __init__(self,elements,joins=None,dochecks=True):
        """
        :param elements: Pipeline elements for this pipeline.
        :type elements: :class:`PipelineElement` objects
        
        :param joins:
            Used to populate the specified element's attribute by extracting
            the output of the supplied pipeline.
        :type joins: 
            Sequence of 3-tuples (elementindex,attribute,Pipelineobject) 
        
        :param dochecks: If True, validates the inputs as the correct type.
        :type dochecks: bool
        """
        
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
        Feed an initial datum into the first stage of the pipeline.
        """
        self.datadeques[0].appendleft(data)
    
    def feedIter(self,dataiter):
        """
        Feed initial data into the first stage of the pipeline.  The input
        will be treated as an iterator and each item will be added to
        the pipeline.
        """
        for dat in dataiter:
            self.datadeques[0].appendleft(dat)
    
    
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
        """
        Process a particular stage once, possibly processing earlier stages
        if there is nothing in the request stage's input.
        
        :param stagenum: Stage number index to process
        :type stagenum: int
        """
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
        Process the pipeline to the end.
        
        :param repeat: 
            If True, the pipeline will by processed continually. Otherwise, only
            one step through the pipeline will be run.
        :type repeat: bool
        
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
        self.datadeques[0].append(input)
        for stage in range(len(self.elements)):
            self.processStage(stage)
            while self.cycles[stage]:
                if processinglimit and self.cycles[stage] >= processinglimit:
                    if removefailed:
                        raise NotImplementedError
                    raise PipelineError('hit processing limit at stage %i'%stage)
                self.processStage(stage)
            self.datadeques[stage+1].append(self.datadeques[stage].popleft())
                
        return self.datadeques[-1].pop()
        
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
        
        :returns: The number of times all of the stages were fed and processed.
        
        :except PipelineError: 
            If an element has not completed after `repeat` processing cycles
            have been run.
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
                    self.processStage(st)
                    while self.cycles[st]:
                        if repeat is not True and self.cycles[st] >= repeat:
                            raise PipelineError('hit processing limit at stage %i'%st)
                        self.processStage(st)
                    count+=1   
            else:
                if len(dd)>0:
                    self.processStage(st)
                    if not self.cycles[st]:
                        count=1
            counts.append(count)
        return counts
    
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
            
del ABCMeta,abstractmethod,abstractproperty #clean up namespace