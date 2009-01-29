#Erik Tollerud (etolleru@uci.edu) 2008

"""
This module contains objects and functions for generating catalogs of objects
where derived quantities are dynamically updated as they are changed.

The eventual plan is to include modules to also dynamically update via a web
server.

(this package is currently not in a useable state - heavy development is under 
way)
"""

from __future__ import division
from math import pi
import numpy as np
from .coords import AngularPosition,AngularCoordinate

class Datum(object):
    """
    Base for all forms of data
    """
    def __init__(self):
        raise NotImplementedError
    
class DerivedDatum(Datum):
    """
    This is a data point derived from other datum.
    """
    def __init__(self):
        raise NotImplementedError
    
class MeasuredDatum(Datum):
    """
    This is a measured data point with a source.
    """
    def __init__(self):
        raise NotImplementedError
    
class DataItem(object):
    """
    This represents a particular property of an object as a set of Datum objects.
    """
    def __init__(self):
        raise NotImplementedError
    

    
class AstronomicalObject(object):
    def __init__(self):
        raise NotImplementedError


#class MeasuredValue(object):    
#    __slots__ = ('__val','__err','__src')
        
        
#    def __getVal(self):
#        if callable(self.__val):
#            return self.__val()
#        else:
#            return self.__val
#    def __setVal(self,val):
#        self.__val = val
#    value=property(__getVal,__setVal)
    
#    def __getErr(self):
#        if callable(self.__err):
#            return self.__err()
#        else:
#            return self.__err
#    def __setErr(self,val):
#        self.__err = val
#    error=property(__getErr,__setErr)
    
#    def __getSrc(self):
#        return self.__src
#    def __setSrc(self,val):
#        self.__src = str(val)
#    source=property(__getSrc,__setSrc)
    
#    def __init__(self,value,source=None,error=None):
#        from weakref import ref
#        from operator import isSequenceType
        
#        self.value = value
#        if error is None:
#            self.error=0
#        else:
#            self.error=error
#        if source is None:
#            self.source = ''
#        else:
#            self.source = source
                
#class DerivedValue(object):
#    __slots__=('__val','__err','__srcs')
#    #TODO:auto-derive
    
#    def _wrdelcb(self,wr):
#        try:
#            self.__srcs.remove(wr)
#        except ValueError:
#            pass #already removed
        
#    def __getVal(self):
#        return self.__val
#    def __setVal(self,val):
#        if not callable(val):
#            self.__val = val
#        else:
#            raise NotImplementedError
#    value=property(__getVal,__setVal)

#    def __getErr(self):
#        return self.__err
#    def __setErr(self,val):
#        if not callable(val):
#            self.__err = val
#        else:
#            raise NotImplementedError
#    error=property(__getErr,__setErr)

#    def __getSrcs(self):
#        return [wr() for wr in self.__srcs]
#    def __setSrcs(self,val):
#        from weakref import ref
#        from operator import isSequenceType
        
#        if not isSequenceType(sources):
#            sources=[sources]
#        sources=[s if type(s) == ref else ref(s) for s in val]
#    sources=property(__getSrcs,__setSrcs)

#    def __init__(self,value,sources=None,error=None):
#        """
#        value may either be a number/object, in which case sources and error
#        can be undefined or an arbitrary sequence of sources and an arbitrary
#        error value
#        OR
#        value can be a callable, in which case the number of arguments must
#        equal the number of sources, and error must also be a callable
#        """
#        if source is None:
#            self.sources=[]
#        else:
#            self.sources = [s for s in sources]
            
#        self.value = value
            
#        if error is None:
#            self.error=0
#        else:
#            self.error = error
            
#class MeasuredValues(object):
#    __slots__=('__vals','__listeners','__curri')
    
    
#    def __findIndex(self,iorsrc):
#        if type(val)!=int:
#            for i,v in enumerate(self.__vals):
#                if v.source == iorsrc:
#                    iorsrc=j
#                    break
#            if i==len(self.__vals):
#                raise IndexError('Source %s not found'%iorsrc)
#        return iorsrc
#    def __getitem__(self,i):
#        return self.__vals[self.__findIndex(i)]
#    def __setitem__(self,i,val):
#        i = self.__findIndex(i)
#        self.__vals[i] = MeasuredValue(val)
#    def __delitem__(self,i):
#        i = self.__findIndex(i)
#        del self.__vals[i]
#    def __len__(self):
#        return len(self.__vals)
#    def __iter__(self):
#        for v in self.__vals:
#            yield v.value,v.error,v.source
    
#    def __getCurr(self):
#        return self.__vals[self__curri].value
#    def __setCurr(self,val):
#        self.__curri=self.__findIndex(val)
#        self.__notify()
#    currentvalue=property(__getCurr,__setCurr)
#    def __getCurrErr(self):
#        return self.__vals[self__curri].error
#    currenterror=property(__getCurrErr)
    
#    def __getCurrSrc(self):
#        return self.__vals[self__curri].source
#    currentsource=property(__getCurrSrc)
    
#    def __notify(self):
#        for l in self.__listeners:
#            l(self)
    
#    def addValue(self,*args):
#        """
#        args can be a MeasuredValue, or arguments will be passed into
#        MeasuredValue constructor
#        """
#        if len(args)==1 and type(args[0]) == MeasuredValue:
#            self.__vals.append(args[0])
#        else:
#            self.__vals.append(MeasuredValue(*args))
#    def updateValue(self,*args):
#        self.addValue(*args)
#        self.__curri=len(self.__vals)-1
#        self.__notify()
        
#    def addListener(self,listener,callim=True):
#        """
#        add a function to call when the current value changes
#        listener signature is f(MeasuredValues)
#        if callim is True, the listener will be called once, immediately
#        """
#        if not callable(listener):
#            raise ValueError('supplied listener not callable')
#        if callim:
#            listener(self)
#        self.__listeners.append(listener)
#    def removeListener(self,listener):
#        self.__listeners.remove(listener)
        
#    def listValues(self):
#        return [mv.value for mv in self.__vals]
#    def listErrors(self):
#        return [mv.error for mv in self.__vals]
#    def listSources(self):
#        return [mv.source for mv in self.__vals]
    
    
#    def __init__(self,*args,**kwargs):
#        """
#        possible forms:
#        MeasuredValues(MVs)
#        MeasuredValues(vals)
#        MeasuredValues(vals,sources)
#        MeasuredValues(vals,sources,errors)
#        """
#        from operator import isSequenceType
#        if len(args) > 3 or len(args) < 1:
#            raise ValueError('Unexpected number of arguments')
        
#        if len(args) <= 1:
#            if not isSequenceType(args[0]):
#                vals=[args[0]]
#            else:
#                vals=args[0]
#            bmvs=[type(v) == MeasuredValue for v in vals]
#            if np.any(bmvs):
#                if not np.all(bmvs):
#                    raise ValueError('mixed MeasuredValues with other values')
#                kwargs['MVs']=vals
#            elif 'vals' in kwargs:
#                raise ValueError('double-definition of vals')
#            else:
#                kwargs['vals']=vals
            
#        if len(args) <= 2:
#            if not isSequenceType(args[1]):
#                vals=[args[1]]
#            else:
#                vals=args[1]
#            if 'sources' in kwargs:
#                raise ValueError('double-definition of sources')
#            else:
#                kwargs['sources']=vals
                
#        if len(args) <= 3:
#            if not isSequenceType(args[2]):
#                vals=[args[2]]
#            else:
#                vals=args[2]
#            if 'errors' in kwargs:
#                raise ValueError('double-definition of errors')
#            else:
#                kwargs['errors']=vals
        
        
#        self.__vals=vals=[]
#        self.__listeners=[]
#        self.__curri=0
#        if 'MVs' in kwargs:
#            if 'sources' in kwargs or 'errors' in kwargs or 'values' in kwargs:
#                raise ValueError('mixed MeasuredValues with other values')
#            for mv in kwargs['MVs']:
#                vals.append(mv)
#        else:
#            if 'values' not in kwargs:
#                raise ValueError('no values defined')
#            if 'sources' not in kwargs or kwargs['sources'] is None:
#                kwargs['sources']=['' for v in kwargs['values']]
#            if 'errors' not in kwargs or kwargs['errors'] is None:
#                kwargs['sources']=[0 for v in kwargs['values']]
                
#            for v,s,e in zip(kwargs['values'],kwargs['sources'],kwargs['errors']):
#                vals.append(MeasuredValue(v,s,e))
                
#        self.__curri=len(vals) - 1

            
#class AstronomicalObject(object):
#    def __get_pos(self):
#        return self.__pos.getCurrent()
#    def __set_pos(self,val):
#        if type(val)==tuple and len(val)==2:
#            val=AngularPosition(val[0]),val[1]
#        else:
#            val=AngularPosition(val),None
#        self.__pos.setCurrent(val)
#    def __del_pos(self):
#        self.__pos.delCurrent()
#    position=property(__get_pos,__set_pos,__del_pos)
#    def __get_dist(self):
#        return self.__dist[self.__currdisti].val
#    def __set_dist(self,val):
#        self.__pos.setCurrent(val)
#    def __del_dist(self):
#        self.__dist.delCurrent()
#    distance=property(__get_dist,__set_dist,__del_dist)
    
#    def choosePosition(self,source):
#        self.__pos.choose(source)       
#    def chooseDistance(self,source):
#        self.__dist.choose(source)
#    def listPositions(self):
#        return [str(p) for p in self.__pos]
#    def listDistances(self):
#        return self.__dist.vals[:]
#    def listPositionSources(self):
#        return self.__pos.sources[:]
#    def listDistanceSources(self):
#        return self.__dist.sources[:]
#    def addPosition(self,val,source='unknown',updatecurr=True):
#        return self.__pos.add(val,source,updatecurr)
#    def addDistance(self,val,source='unknown',updatecurr=True):
#        return self.__dist.add(val,source,updatecurr)
    
#    def __init__(self,*args,**kwargs):
#        """
#        possible forms:
#        *AstronomicalObject()
#        *AstronomicalObject(position)
#        *AstronomicalObject(position,distance)
#        *AstronomicalObject(ra=val,dec=val | l=val,b=val):
        
#        optional keywords:
#        'epoch': J2000 or B1950
#        'distance' if AstronomicalObject(position,distance) is not used
#        'positionsource' or 'ps':source for angular position
#        'distancesource' or 'ds':source for distance measurement
#        """
#        if 'epoch' in kwargs:
#            epoch=kwargs.pop('epoch')
#        else:
#            epoch='J2000'
#        if 'distance' in kwargs:
#            dist=kwargs.pop(distance)
#        else:
#            dist=None
                
#        if len(args)==0:
#            if 'ra' in kwargs and 'dec' in kwargs:
#                ra=AngularCoordinate(kwargs.pop('ra'),epoch=epoch)
#                dec=AngularCoordinate(kwargs.pop('dec'),epoch=epoch)
#            elif 'l' in kwargs and 'b' in kwargs:
#                ra,dec=galactic_to_equatorial(l,b)
#            else:
#                ra,dec=AngularCoordinate(epoch=epoch),AngularCoordinate(epoch=epoch)
#            pos=AngularPosition(ra,dec)
                
                
#        elif len(args)==1:
#            pos=AngularPosition(args[0])
#        elif len(args)==2:
#            pos=AngularPosition(args[0])
#            if dist is not None:
#                raise ValueError('Got two values for distance')
#            dist=args[1]
#        else:
#            raise TypeError('Cannot pass more than 2 arguments to AstronomicalObject constructor')
        
#        if 'positionsource' in kwargs:
#            ps=kwargs.pop('positionsource')
#        elif 'ps' in kwargs:
#            ps=kwargs.pop('ps')
#        else:
#            ps=None
#        self.__pos = MeasuredValues(pos,source=ps)
        
#        if 'distancesource' in kwargs:
#            ds=kwargs.pop('distancesource')
#        elif 'ds' in kwargs:
#            ds=kwargs.pop('ds')
#        else:
#            ds=None
            
#        self.__dist = MeasuredValues(dist,source=ds)
        
#        if len(kwargs) != 0:
#            raise TypeError('Unrecognized keyword'+('s ' if len(kwargs)>1 else ' ')+','.join(kwargs.keys()))    
