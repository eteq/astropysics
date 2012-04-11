#Copyright 2008 Erik Tollerud
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
This module contains classes representing coordinates in spatial, celestial, and
terrestrial coordinate systems, as well as implementations of transformations
between many of the coordinate systems. There are also utility functions for 
easy-of-use of some of these classes.

The coordinate system framework is designed to allow users to add their own
coordinate systems easily, if desired, and implement transformations between
theirs and the builtin coordinate systems. This is implemented through the
transformation registry static methods of the :class:`CoordinateSystem` class
(e.g. :meth:`CoordinateSystem.registerTransform`).

Examples
^^^^^^^^
A common task might be to convert coordinates from one standard coordinate
system to another. The coords module makes this simple::

    >>> from astropysics.coords import ICRSCoordinates,GalacticCoordinates
    >>> gcoords = ICRSCoordinates('2h34m12.32s',12.3).convert(GalacticCoordinates)
    >>> print gcoords
    GalacticCoordinates: l=158.558650,b=-43.350066
    >>> print '%.3f'%gcoords.l.degrees
    158.559
    >>> print '%.3f'%gcoords.l.radians
    2.767
    >>> print gcoords.b.getDmsStr(canonical=True)
    -43:21:00.24
    
Note the initial input composed of an hours,minutes,seconds string input for the
RA, and a float for the dec -- :class:`EquatorialCoordinate` objects contain a
powerful parsing system that accepts most standard astronomical representations.
The resulting :class:`GalacticCoordinates` object's coordinates can then be
accessed in any of the various ways supported by the :class:`AngularCoordinate`
object.

.. warning::
    Errors are not currently propogated in all coordinate transforms - this
    will be corrected eventually, but for now you should check to make sure 
    errors propogate for any coordinate transform you want to perform.


{transformdiagram}
Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.coords.coordsys
   :parts: 1
   
Module API
^^^^^^^^^^

"""

#TODO: implement polar motion lookup techniques

from __future__ import division,with_statement

from ..constants import pi
from ..utils import add_docs

import numpy as np
_twopi = 2*pi
_pio2 = pi/2

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


class AngularCoordinate(object):
    """
    A class representing an angular value.
    
    Arithmetic operators can be applied to the coordinate, and will be applied 
    directly to the numerical value in radians.  For + and -, two angular 
    coordinates may be used, although for -, an AngularSeparation object will
    be returned.
    
    """
    import re as _re
    __slots__=('_decval','_range')
    
    #this disturbingly complex RE matches anything that looks like a standard sexigesimal or similar string
    __acregex = _re.compile(r'(?:([+-])?(\d+(?:[.]\d*)?)(hours|h|degrees|d|radians|rads|rad|r| |:(?=\d+:\d+[.]?\d*$)))?(?:(\d+(?:[.]\d*)?)(m|\'|[:]| ))?(?:(\d+(?:[.]\d*)?)(s|"|$))?$')
    #and this one matches all things that look like raw numbers
    __decregex = _re.compile(r'[+-]?\d+([.]\d*)?$') 
    
    def __init__(self,inpt=None,sghms=None,range=None,radians=False):
        """
        The input parser is very adaptable, and can be in any of the following 
        forms for `inpt`:
        
        * A float value
            if `radians` is True, this will be interpreted as decimal radians,
            otherwise, it is in degrees.
        * An :class:`AngularCoordinate` object
            A copy of the input object will be created.
        * None
            The default of 0 will be used.
        * A 3-tuple
            If `sghms` is True, the tuple will be interpreted as
            (hours,min,sec), otherwise, (degrees,min,sec).
        * A string of the form ##.##
            If `radians` is True, this will be cast to a float and used as
            decimal radians, otherwise, it is in degrees.
        * A string of the form ##.##d or ##.##degrees
            The numerical part will be cast to a float and used as degrees.
        * A string of the form ##.##h or ##.##hours
            The numerical part will be cast to a float and used as hours.
        * A string of the form ##.##radians,##.##rads, or ##.##r
            The numerical part will be cast to a float and used as radians.
        * A string of the form (+/-)##h##m##.##s
            The numerical parts will be treated as hours,minutes, and seconds.
        * A string of the form (+/-)##d##m##.##s or (+/-)##d##'##.##"
            The numerical parts will be treated as degrees,minutes, and seconds.
        * A string of the form (+/-)##:##:##.## or (+/-)## ## ##.##
            Sexigesimal form. If `sghms` is None the presence of a a + or - sign
            idicates that it should be interpreted as degrees, minutes, and
            seconds. If the sign is absent, the numerical portions will be
            treated as hours,min,sec. thewise, if `sghms` evaluates to True, the
            numerical parts will be treated as hours,minutes, and seconds, and
            if `sghms` evaluates to False, degrees,minutes, and seconds.
        
        :param inpt: The coordinate value -- valid forms are described above.            
        :param sghms: 
            If True, ambiguous sexigesimal inputs should be hours, minutes, and
            seconds instead of degrees,arcmin, and arcsec
        :type sghms: boolean
        :param range: 
            Sets the valid range of coordinates.  Either a
            2-sequence (lowerdegrees,upperdegrees) or None (for no limit)
        :param radians:
            If True, ambiguous inputs are treated as radians rather than
            degrees.
        :type radians: boolean
        
        **Examples**
        
        >>> from math import pi
        >>> ac = AngularCoordinate(2.5)
        >>> print ac
        +2d30'00.00"
        >>> print AngularCoordinate(ac)
        +2d30'00.00"
        >>> print AngularCoordinate(pi,radians=True)
        +180d00.00"
        >>> print AngularCoordinate('1.1')
        +1d6'00.00"
        >>> print AngularCoordinate('1.1',radians=True)
        +63d1'31.29"
        >>> print AngularCoordinate('12d25m12.5s')
        +12d25'12.50"
        >>> print AngularCoordinate('3:30:30',sghms=True)
        +52d37'30.00"
        >>> print AngularCoordinate('3:30:30',sghms=False)
        +3d30'30.00"
        >>> print AngularCoordinate('-3:30:30',sghms=None)
        -3d30'30.00"
        >>> print AngularCoordinate('+3:30:30',sghms=None)
        +3d30'30.00"
        >>> print AngularCoordinate('3:30:30',sghms=None)
        +52d37'30.00"
        
        """
        from operator import isSequenceType
        
        self._range = None
        
        if isinstance(inpt,AngularCoordinate):
            self._decval = inpt._decval
            self._range = inpt._range
            return
        elif inpt is None:
            self._decval = 0
        elif isinstance(inpt,basestring):
            sinpt = inpt.strip()
            
            decm = self.__decregex.match(sinpt)
            if decm:
                if radians:
                    self.radians = float(decm.group(0))
                else:
                    self.degrees = float(decm.group(0))
            else:
                acm = self.__acregex.match(sinpt)
                if acm:
                    sgn,dec1,mark1,dec2,mark2,dec3,mark3 = acm.group(1,2,3,4,5,6,7)
                    val = (0 if dec1 is None else float(dec1)) + \
                          (0 if dec2 is None else float(dec2)/60) + \
                          (0 if dec3 is None else float(dec3)/3600)
                    if sgn == '-':
                        val *= -1
                    if mark1 == ':' or mark1 == ' ':
                        if sghms is None:
                            if sgn is None:
                                self.hours = val
                            else: #'+' or '-'
                                self.degrees = val
                        elif sghms:
                            self.hours = val
                        else:
                            self.degrees = val
                    elif mark1 == 'hours' or mark1 == 'h':
                        self.hours = val
                    elif mark1 == 'degrees' or mark1 == 'd':
                        self.degrees = val
                    elif mark1 == 'radians' or mark1 == 'rad' or mark1 == 'rads' or mark1=='r':
                        self.radians = val
                    else:
                        try:
                            if radians:
                                self.radians = float(val)
                            else:
                                self.degrees = float(val)
                        except ValueError:
                            raise ValueError('invalid string input for AngularCoordinate')
                else:
                    raise ValueError('Invalid string input for AngularCoordinate: '+inpt)
            
        elif isSequenceType(inpt) and len(inpt)==3:
            if sghms:
                self.hrsminsec = inpt
            else:
                self.degminsec = inpt
        elif radians:
            self._decval = float(inpt)
        else:
            from math import radians
            self._decval = radians(inpt)
            
        self.range = range
    
    def _setDegminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set degminsec as a length-3 iterator')
        self.degrees = abs(dms[0])+abs(dms[1])/60.+abs(dms[2])/3600.
        if dms[0]<0:
            self._decval*=-1
    def _getDegminsec(self):
        fulldeg = abs(self.degrees)
        deg = int(fulldeg)
        fracpart = fulldeg-deg
        min = int(fracpart*60.)
        sec = fracpart*3600.-min*60.
        return -deg if self.degrees < 0 else deg,min,sec
    degminsec = property(_getDegminsec,_setDegminsec,doc="""
    The value of this :class:`AngularCoordinate` as an (degrees,minutes,seconds)
    tuple, with degrees and minutes as integers and seconds as a float.    
    """)
    dms = degminsec
    
    def _setHrsminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set hrsminsec as a length-3 iterator')
        self.degrees = 15*(dms[0]+dms[1]/60.+dms[2]/3600.)
    def _getHrsminsec(self):
        factorized = self.degrees/15.
        hrs = int(factorized)
        mspart = factorized - hrs
        min = int(mspart*60.)
        sec = mspart*3600.-min*60.
        return hrs,min,sec
    hrsminsec = property(_getHrsminsec,_setHrsminsec,doc="""
    The value of this :class:`AngularCoordinate` as an (hours,minutes,seconds)
    tuple, with hours and minutes as integers and seconds as a float.    
    """)
    hms = hrsminsec
    
    def _setDecdeg(self,deg):
        rads = deg*pi/180.
        if self.range is not None:
            rads = self._checkRange(rads)
        self._decval = rads
    def _getDecdeg(self):
        return self._decval*180/pi
    degrees = property(_getDecdeg,_setDecdeg,doc="""
    The value of this :class:`AngularCoordinate` in decimal degrees.   
    """)
    d = degrees    
    
    def _setRad(self,rads):
        if self.range is not None:
            rads = self._checkRange(rads)
        self._decval = rads
    def _getRad(self):
        return self._decval
    radians = property(_getRad,_setRad,doc="""
    The value of this :class:`AngularCoordinate` in decimal radians.   
    """)
    r = radians
    
    def _setDechr(self,hr):
        rads = hr*pi/12
        if self.range is not None:
            rads = self._checkRange(rads)
        self._decval = rads
    def _getDechr(self):
        return self._decval*12/pi
    hours = property(_getDechr,_setDechr,doc="""
    The value of this :class:`AngularCoordinate` in decimal hours.   
    """)
    h = hours
    
    def _checkRange(self,rads):
        """
        Checks if the input value is in range - returns the new value, or raises
        a :exc:`ValueError`.
        """
        if self._range is not None:
            low,up,cycle = self._range
            if cycle is None:
                if low <= rads <= up:
                    return rads 
                else:
                    raise ValueError('Attempted to set angular coordinate outside range')
            else:
                if cycle > 0:
                    #this means use "triangle wave" pattern with the given quarter-period 
                    from math import sin,asin
                    offset = low/(low-up)-0.5
                    return (up-low)*(asin(sin(pi*(2*rads/cycle+offset)))/pi+0.5)+low
                else:
                    return (rads-low)%(up-low)+low
        else:
            return rads
    def _setRange(self,newrng):
        oldrange = self._range        
        try:
            if newrng is None:
                self._range = None
            else:
                from math import radians
                newrng = tuple(newrng)
                if len(newrng) == 2:
                    if newrng[1]-newrng[0] == 360:
                        newrng = (newrng[0],newrng[1],0)
                    else:
                        newrng = (newrng[0],newrng[1],None)
                elif len(newrng)==3:
                    pass
                else:
                    raise TypeError('range is not a 2 or 3-sequence')
                if newrng[0] > newrng[1]:
                    raise ValueError('lower edge of range is not <= upper')
                
                newrng = ( radians(newrng[0]),radians(newrng[1]), \
                           None if newrng[2] is None else radians(newrng[2]) )
            self._range = newrng
            self._decval = self._checkRange(self._decval)
        except ValueError,e:
            self._range = oldrange
            if e.args[0] == 'lower edge of range is not <= upper':
                raise e
            else:
                raise ValueError('Attempted to set range when value is out of range')
    def _getRange(self):
        if self._range is None:
            return None
        else:
            from math import degrees
            if self._range[2] is None:
                return degrees(self._range[0]),degrees(self._range[1])
            else:
                return degrees(self._range[0]),degrees(self._range[1]),degrees(self._range[2])
    range = property(_getRange,_setRange,doc="""
    The acceptable range of angles for this :class:`AngularCoordinate`.  This can
    be set as a 2-sequence (lower,upper), or as a 3-sequence (lower,upper,cycle), 
    where cycle can be :
    
        * 0: Angle values are coerced to lie in the range (default for 
          2-sequence if upper-lower is 360 degrees)
        * None: A :exc:`ValueError` will be raised if out-of-range (default for 
          2-sequence otherwise)
        * A positive scalar: Values are coerced in a triangle wave scheme, with
          the scalar specifying the period. e.g. for the latitude, (-90,90,360)
          would give the correct behavior)
    """)
   
    def __str__(self):
        return self.getDmsStr(sep=('d',"'",'"'))
            
    def __eq__(self,other):
        if hasattr(other,'_decval'):
            return self._decval==other._decval
        else:
            return self._decval==other
        
    def __ne__(self,other):
        return not self.__eq__(other)
        
    def __add__(self,other):
        if hasattr(other,'_decval'):
            res = self.__class__()
            res._decval = self._decval + other._decval
        else:
            res = self.__class__()
            res._decval = self._decval + other
        return res
        
    def __sub__(self,other):
        
        if isinstance(other,AngularCoordinate):
            from math import degrees
            res = AngularSeparation(degrees(other._decval),degrees(self._decval))
        else:
            res = AngularCoordinate()
            res._decval = self._decval - other
        return res
        
    def __mul__(self,other):
        res = self.__class__()
        res._decval = self._decval*other
        return res
        
    def __div__(self,other):
        res = self.__class__()
        res._decval = self._decval/other
        return res
        
    def __truediv__(self,other):
        res = self.__class__()
        res._decval = self._decval//other
        return res
        
    def __pow__(self,other):
        res = self.__class__()
        res._decval = self._decval**other
        return res
        
    def __float__(self):
        return self.degrees
        
    def getDmsStr(self,secform='%05.2f',sep=(unichr(176),"'",'"'), sign=True, 
                  canonical=False, inclzero=True):
        """
        Generates the string representation of this AngularCoordinate as
        degrees, arcminutes, and arcseconds.
        
        :param secform: a formatter for the seconds
        :type secform: string
        :param sep:
            The seperator between components - defaults to degree sign, ' and "
            symbols.
        :type sep: string or 3-tuple of strings
        :param sign: Forces sign to be present before degree component.
        :type sign: boolean
        :param canonical: forces [+/-]dd:mm:ss.ss , overriding other arguments
        :param inclzero: 
            If True, a "0" is included whenever even if the degrees or minutes 
            are 0.  Otherise, the "0" and the corresponding separator are 
            omitted from the string. 
        :type inclzero: bool
        
        :returns: String representation of this object.
        """
        d,m,s = self.degminsec
        
        if canonical:
            secform = '%05.2f'
            sep = (':',':','')
            sign = True
        
        s = secform%s
        
        if int(float(s))>=60:
            s = secform%0
            m += 1
            if m==60:
                m = 0
                d += 1
        
        d,m=str(abs(d)),str(m)
        
        if isinstance(sep,basestring):
            if sep == 'dms':
                sep = ('d','m','s')
            sep = (sep,sep)
        
        tojoin = []
        
        if sign and self._decval  >= 0:
            tojoin.append('+')
        if self._decval<0:
            tojoin.append('-')
        
        if inclzero or d is not '0':
            tojoin.append(d)
            tojoin.append(sep[0])
                
        if inclzero or m is not '0':
            tojoin.append(m)
            tojoin.append(sep[1])
                
        tojoin.append(s)
        if len(sep)>2:
            tojoin.append(sep[2])
            
        return ''.join(tojoin)
        
    def getHmsStr(self,secform = None,sep = ('h','m','s'), canonical = False, 
                  inclzero=True):
        """
        gets the string representation of this AngularCoordinate as hours,
        minutes, and seconds
        
        secform is the formatter for the seconds component
        
        sep is the seperator between components - defaults to h, m, and s
        
        canonical forces [+/-]dd:mm:ss.ss , overriding other arguments
        
        Generates the string representation of this AngularCoordinate as hours,
        minutes, and seconds.
        
        :param secform: a formatter for the seconds component
        :type secform: string
        :param sep:
            The seperator between components - defaults to 'h', 'm', and 's'.
        :type sep: string or 3-tuple of strings
        :param canonical: forces [+/-]dd:mm:ss.ss , overriding other arguments
        :param inclzero: 
            If True, a "0" is included whenever even if the degrees or minutes 
            are 0.  Otherise, the "0" and the corresponding separator are 
            omitted from the string. 
        :type inclzero: bool
        
        :returns: String representation of this object.
        """
        
        h,m,s = self.hrsminsec
        
        if canonical:
            secform = '%05.2f'
            sep = (':',':','')
        
        if secform is None:
            s = str(s)
        else:
            s = secform%s
        
        if int(float(s))>=60:
            s = secform%0
            m += 1
            if m==60:
                m = 0
                h += 1
                if h==24:
                    h = 0
            
        h,m=str(h),str(m)
        
        if isinstance(sep,basestring):
            if sep == 'hms':
                sep = ('h','m','s')
            sep = (sep,sep)
        
        tojoin = []
        
        
        if inclzero or h is not '0':
            tojoin.append(h)
            tojoin.append(sep[0])
        
        if inclzero or m is not '0':
            tojoin.append(m)
            tojoin.append(sep[1])
                
        tojoin.append(s)
        if len(sep)>2:
            tojoin.append(sep[2])
            
        return ''.join(tojoin)

class AngularSeparation(AngularCoordinate):
    """
    This class represents a separation between two angular coordinates on the
    unit sphere.
    
    A constructor is available, but the most natural way to generate this object
    is to use the subtraction (-) operator on two :class:`AngularCoordinate`
    objects or two :class:`LatLongCoordinates` objects.
    """
    
    def __init__(self,*args):
        """
        Input arguments can be either:
        
        * AngularSeparation(:class:`AngularSeparation` object) 
            Generates a copy of the provided object.
        * AngularSeparation(sep) 
            Generates a separation of the provided distance with no starting point.
        * AngularSeparation(start,end) 
            Computes the separation from the start and end objects, which must
            be :class:`AngularCoordinate` objects.
          
        """
        if len(args) == 1:
            a = args[0]
            if a.__class__ == self.__class__:
                self._decval = args[0]._decval
                self._range = args[0]._range
                return
            
            sep = a._decval if hasattr(a,'_decval') else a
        elif len(args) == 2:
            a0,a1 = args
            a0 = a0._decval if hasattr(a0,'_decval') else a0
            a1 = a1._decval if hasattr(a1,'_decval') else a1
            sep = a1 - a0
        else:
            raise ValueError('improper number of inputs to AngularSeparation')
        
        AngularCoordinate.__init__(self,sep)
        
    def __add__(self,other):
        if isinstance(other,AngularCoordinate) and not self.__class__ == other.__class__:
            res = AngularCoordinate()
            res._decval = self._decval+other._decval
            return res
        else:
            return AngularCoordinate.__add__(self,other)
            
    #comparisons
    def __lt__(self,other):
        return self._decval < other._decval
    def __le__(self,other):
        return self._decval <= other._decval
    def __gt__(self,other):
        return self._decval > other._decval
    def __ge__(self,other):
        return self._decval >= other._decval
    def __eq__(self,other):
        return self._decval == other._decval
    def __ne__(self,other):
        return self._decval != other._decval
        
    def _getArcsec(self):
        return self.degrees*3600
    def _setArcsec(self,val):
        self.degrees = val/3600
    arcsec = property(_getArcsec,_setArcsec,doc=None)
    
    def _getArcmin(self):
        return self.degrees*60
    def _setArcmin(self,val):
        self.degrees = val/60
    arcmin = property(_getArcmin,_setArcmin,doc=None)
    
        
    def projectedSeparation(self,zord,usez=False,**kwargs):
        """
        Computes the physical projected separation assuming a given distance.
        
        kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
        
        :param zord: Redshift or distance
        :type zord: scalar number
        :param usez:
            If True, the input will be interpreted as a redshift, and kwargs
            will be passed into the distance calculation. The result will be in
            pc. Otherwise, `zord` will be interpreted as a distance.
        :type usez: boolean
        
        :returns: a float value for the separation (in pc if redshift is used) 
        """
        from .funcs import angular_to_physical_size
        
        return angular_to_physical_size(self.arcsec,zord,usez=usez,**kwargs)
    
    def separation3d(self,zord1,zord2,usez=False,**kwargs):
        """
        computes the 3d separation assuming the two points at the ends of this
        :class:`AngularSeparation` are at the distances `zord1` and `zord2`.  
        
        :param zord1: Redshift or distance for start point
        :type zord1: scalar number
        :param zord2: Redshift or distance for end point
        :type zord2: scalar number
        :param usez:
            If True, the inputs will be interpreted as a redshift, and kwargs
            will be passed into the distance calculation. The result will be in
            pc. Otherwise, `zord` will be interpreted as a distance.
        :type usez: boolean
        
        :returns: a float value for the separation (in pc if redshift is used) 
        """
        from math import sin,cos,sqrt
        
        if usez:
            d1 = cosmo_z_to_dist(zord1,disttype=2,**kwargs)*1e6 #pc
            d2 = cosmo_z_to_dist(zord1,disttype=2,**kwargs)*1e6 #pc
        else:
            if len(kwargs)>0:
                raise TypeError('if not using redshift, kwargs should not be provided')
            d1 = zord1
            d2 = zord2
    
        costerm = 2*d1*d2*cos(self._decval)
        return sqrt(d1*d1+d2*d2-costerm)
    
#<-----------------------------Coordinate systems------------------------------>

class _CoosysMeta(ABCMeta):
    """
    Metaclass for CoordinateSystem class and subclasses - needed to support 
    :class:`CoordinateSystem.registerTransform` decorator.
    """
    def __init__(cls,name,bases,dct):
        ABCMeta.__init__(cls,name,bases,dct)
        import inspect
        
        for k,v in inspect.getmembers(cls):
            if isinstance(v,_TransformerMethodDeco):
                for vfc,vtc in zip(v.fromclasses,v.toclasses):
                    fromclass = cls if vfc == 'self' else vfc
                    toclass = cls if vtc == 'self' else vtc
                    CoordinateSystem.registerTransform(fromclass,toclass,v.f,v.transtype)
                setattr(cls,k,staticmethod(v.f))
                
class _TransformerMethodDeco(object):
    """
    A class representing methods used for registering transforms  for the class
    the are in.
    """
    def __init__(self,f,fromclass,toclass,transtype=None):
        self.f = f
        self.fromclasses = [fromclass]
        self.toclasses = [toclass]
        self.transtype = transtype


#Default for optmizing convert functions - currently false because it's not smart enough
_convertoptimizedefault = False

class CoordinateSystem(object):
    """
    Base class of all coordinate systems. This class also holds the static
    methods that manage conversion between coordinate systems.
    
    *Subclassing*
    
    * Subclasses of :class:`CoordinateSystem` must override :meth:`__init__` to
      set initial values. 
    
    * :class:`CoordinateSystem` objects are intended to be quite small, so
      unless there is a reason to do otherwise, subclasses should have a
      :attr:`__slots__` class attribute (it should be a sequence of all the 
      valid attribute names for the object - see
      http://docs.python.org/reference/datamodel.html for an explanation of the
      `__slots__` mechanism).
      
    * The :attr:`transweight` class variable can be set to determine the
      weighting of this class when computing coordinate transformation pathways.
      Note that *smaller* weights are preferred paths (e.g. a larger weight is
      less likely to be visited).  See 
      :meth:`CoordinateSystem.getTransformGraph` for more details.
    
    """
    from collections import defaultdict as _defaultdict
    
    __metaclass__ = _CoosysMeta
    __slots__ = tuple()
    
    @abstractmethod
    def __init__(self):
        raise NotImplementedError
    
    _converters = _defaultdict(dict) #first index is from, second is to
    _transtypes = dict()
    
    
    @staticmethod
    def registerTransform(fromclass,toclass,func=None,transtype=None,
                          overwrite=True):
        """
        Register a function to transform coordinates from one system to another.
        
        The transformation function is called is func(fromobject) and should 
        return a new object of type `toclass`.  If called with no arguments, 
        the function should raise a :exc:`NotImplementedError` or behave in a 
        manner defined in subclasses (see e.g. :class:`LatLongCoordinates`).  
        If `transtype` is not None, the output of the transformation function is
        filered through the function applied for that type using
        :meth:`astropysics.CoordinateSystem.addTransType` .
        
        If the transformation function `func` is None, the function is taken to
        be a decorator, and if it is a method of a subclass of
        :class:`CoordinateSystem`, `fromclass` or `toclass` may be the string
        'self' . In this case, the function will use the class itself as the
        from or to class. The function will then be treated as a static method
        (e.g. it should not have `self` as the first argument).
        
        :param fromclass: The class to transform from.
        :type fromclass: subclass of :class:`CoordinateSystem` or 'self'
        :param toclass: The class to transform to.
        :type toclass: subclass of :class:`CoordinateSystem` or 'self'
        :param func: the function to perform the transform or None if decorator.
        :type func: a callable or None
        :param transtype: 
            A transformation type that will be used to determine how the
            transform is performed, or None for no transform type. (see
            :meth:`astropysics.CoordinateSystem.addTransType` for details).
        :type transtype: string or None
        :param overwrite: 
            If True, any already existing function will be silently overriden.
            Otherwise, a ValueError is raised.
        :type overwrite: boolean
        
        **Examples**::
        
            class MyCoordinates(CoordinateSystem):
                ...
            class YourCoordinates(CoordinateSystem):
                ...
            def transformer(mycooobj):
                ...
                return yourcoordobj
            CoordinateSystem.registerTransform(MyCoordinates,YourCoordinates,transformer)
            
            class TheirCoordinates(CoordinateSystem):
                @CoordinateSystem.registerTransform(MyCoordinates,'self')
                @classmethod
                def fromMy(cls,mycoordinates):
                    ...
                    return theircoordobj
        
        """
        if func is None:
            if fromclass != 'self':
                if not issubclass(fromclass,CoordinateSystem):
                    raise TypeError('from class for registerTransform must be a CoordinateSystem')
                
            if toclass != 'self':
                if not issubclass(toclass,CoordinateSystem):
                    raise TypeError('to class for registerTransform must be a CoordinateSystem')
            
            def make_or_extend_trans_meth_deco(f):
                if isinstance(f,_TransformerMethodDeco):
                    f.fromclasses.append(fromclass)
                    f.toclasses.append(toclass)
                elif callable(f):
                    return _TransformerMethodDeco(f,fromclass,toclass,transtype)
                else:
                    raise TypeError('Tried to apply registerTransform to a non-callable')
                
            return make_or_extend_trans_meth_deco
        
        else:
            if not issubclass(fromclass,CoordinateSystem) or not issubclass(toclass,CoordinateSystem):
               raise TypeError('to/from classes for registerTransform must be CoordinateSystems')
           
            if not overwrite and (toclass in CoordinateSystem._converters[fromclass]):
                #format requires 2.6
                #raise ValueError('function already exists to convert {0} to {1}'.format(fromclass,toclass))
                raise ValueError('function already exists to convert %s to %s'%(fromclass,toclass))
            
            if transtype is not None:
                try:
                    ttf = CoordinateSystem._transtypes[transtype]
                except KeyError:
                    raise KeyError('coordinate transformation type %s does not exist'%transtype)
                
                lfunc = lambda cobj:ttf(func(cobj),cobj,toclass)
                lfunc.basetrans = func
                lfunc.transtype = transtype
                CoordinateSystem._converters[fromclass][toclass] = lfunc
            else:
                func.transtype = None
                CoordinateSystem._converters[fromclass][toclass] = func
        
        CoordinateSystem._invalidateTransformCache()
        
        
    @staticmethod
    def getTransform(fromclass,toclass):
        """
        Returns the transformation function to go from `fromclass` to `toclass`.
        """
        return CoordinateSystem._converters[fromclass][toclass]
    
    @staticmethod
    def listAllTransforms():
        """
        Returns a list of 2-tuples (fromclass,toclass) of all the coordinate 
        system combinations that have registered transformation functions.
        """
        trlist = []
        for fr,l in CoordinateSystem._converters.iteritems():
            for li in l:
                trlist.append((fr,li))
        return trlist
    
    @staticmethod
    def listTransformsTo(toclass):
        """
        Returns a list of classes that can be transformed to the supplied class.
        """
        flist = []
        for fr,l in CoordinateSystem._converters.iteritems():
            for li in l:
                if li is toclass:
                    flist.append(fr)
        return flist
        
    @staticmethod
    def listTransformsFrom(fromclass):
        """
        Returns a list of classes that can be transformed from the supplied
        class.
        """
        if fromclass in CoordinateSystem._converters:
            return list(CoordinateSystem._converters[fromclass])
        else:
            return []
        
    @staticmethod
    def delTransform(fromclass,toclass):
        """
        Deletes the transformation function to go from `fromclass` to `toclass`.
        """
        del CoordinateSystem._converters[fromclass][toclass]
        CoordinateSystem._invalidateTransformCache()
        
    @staticmethod
    def addTransType(funcorname):
        """
        Registers a new transformation type. Transformation types are used to
        implement particular types of transformations without repeating similar
        code in each of the actual transformation functions. 
        
        The transformation type function (`funcorname`) will be called as
        transfunc(trans,coords,toclass), where trans is the output of the actual
        registered transformation function, coords is the coordinate object that
        is being converted, and toclass is the target class.
        
        :param funcorname: 
            The function to register as the transfromation type function. If a
            string, the string is taken to be the name to use for the
            transformation (intended for use as a function decorator).
            Otherwise, the transformation type name is taken from the name of
            the function with any intial _ removed.
        :type funcorname: callable or string
        
        :returns: 
            The function itself, to allow for use as a decorator. Note that this
            means that if used as a decorator inside a class, the function will
            be assigned as an *instance* method.  Alternative, if `funcorname` 
            is a string, a function to be called on the transformation function
            will be returned (see second example).
        
        **Examples**::
        
            @addTransType
            def trans1(trans,coord,tocls):
                return tocls(trans*coord.val)
            
            @addTransType('trans1')
            def _transfunc(trans,coord,tocls):
                return tocls(trans*coord.val)
        
        """
        def regtrans(func,typename):
            if typename in CoordinateSystem._transtypes:
                #go through and re-assign all existing transes to use the new one
                for k,v in CoordinateSystem._converters.iteritems():
                    for k2,v2 in v.items():
                        if v2.transtype == typename:
                            btfunc = v2.basetrans
                            coof = lambda cobj:func(btfunc(cobj),cobj,k2)
                            coof.transtype = typename
                            coof.basetrans = btfunc
                            CoordinateSystem._converters[k][k2] = coof
            CoordinateSystem._transtypes[typename] = func
            return func
            
        if isinstance(funcorname,basestring):
            typename = funcorname
            return lambda f:regtrans(f,typename)
        elif callable(funcorname):
            typename = funcorname.func_name
            if typename.startswith('_'):
                typename = typename[1:]
            return regtrans(funcorname,typename)
        else:
            raise TypeError('funcorname is neither a callable nor a string')
    
    @staticmethod    
    def getTransformPath(fromsys,tosys):
        """
        Determines the transformation path from one coordinate system to another
        for use with :meth:`convert`.
        
        :param fromsys: The starting coordinate system class
        :param tosys: The target coordinate system class
        :returns: 
            A list of coordinate classes with the shortest path from `fromsys`
            to `tosys` (*including* `fromsys` and `tosys`) or a callable with
            the transformation if a single-step direct transformation is
            available
        
        :except NotImplementedError: If no path can be found.
        """
        if tosys in CoordinateSystem._converters[fromsys]:
            return CoordinateSystem._converters[fromsys][tosys]
        else:
            failstr = 'cannot convert coordinate system %s to %s'%(fromsys.__name__,tosys.__name__)
            try:
                import networkx as nx
                
                g = CoordinateSystem.getTransformGraph()
                if nx.__version__>'1.4':
                    path = nx.shortest_path(g,fromsys,tosys,weight=True)
                else:
                    path = nx.shortest_path(g,fromsys,tosys,weighted=True)
                if not path:
                    raise NotImplementedError(failstr+'; no transform path could be found')
                return path
            except ImportError,e:
                if e.args[0] == 'No module named networkx':
                    raise NotImplementedError(failstr+'; networkx not installed')
                else:
                    raise
        
    _transgraph = None
    @staticmethod
    def getTransformGraph():
        """
        Returns a `networkx <http://networkx.lanl.gov/>` :class:`DiGraph` object
        representing a graph of the registered coordinate systems and the
        transformations between them.
        
        :except ImportError: If networkx is not installed.
        
        """
        import networkx as nx
        
        if CoordinateSystem._transgraph is None:
            CoordinateSystem._transgraph = g = nx.DiGraph()
            
            transes = []
            for a,b in CoordinateSystem.listAllTransforms():
                avgweight = (getattr(a,'transweight',1) +
                             getattr(b,'transweight',1))/2
                transes.append((a,b,dict(weight=avgweight)))
            g.add_edges_from(transes)
            
            
        return CoordinateSystem._transgraph.copy()
    
    
    _transformcache = _defaultdict(dict)
    @staticmethod
    def _invalidateTransformCache():
        """
        Called when transforms are changed to invalidate the caches
        """
        from collections import defaultdict
        CoordinateSystem._transformcache = defaultdict(dict)
        CoordinateSystem._transgraph = None

    def convert(self,tosys):
        """
        converts the coordinate system from it's current system to a new 
        :class:`CoordinateSystem` object.
        
        :param tosys: 
            The new coordinate system class. Should be a subclass of
            :class:`CoordinateSystem` .
        :returns: A new object of a class determined by `tosys`
        
        :except: raises :exc:`NotImplementedError` if conversion is not present
        """
        
        convpath = CoordinateSystem.getTransformPath(self.__class__,tosys)
        
        if callable(convpath):
            return convpath(self)
        else:
            currobj = self
            currsys = self.__class__
            for intersys in convpath[1:-1]:
                currobj = CoordinateSystem._converters[currsys][intersys](currobj)
                currsys = intersys
            return CoordinateSystem._converters[currsys][tosys](currobj)

class EpochalCoordinates(CoordinateSystem):
    """
    A base class for :class:`CoordinateSystem` classes that have *changeable*
    epochs associated with them. 
    
    *Subclassing*
    
    Subclasses must implement these methods:
    
        * :meth:`__init__` from :class:`CoordinateSystem`
        * :meth:`transformToEpoch` -- see the method's entry for details.
        
    Furthermore, subclasses should set :attr:`julianepoch` to determine if they
    are Julian or Besselian.
    
    
    """    
    #TODO:figure out if there's a way to put this back in to save space - right
    #now if two __slots__ classes are mixed together, this is thrown:
    #TypeError: Error when calling the metaclass bases multiple bases have instance lay-out conflict
    #__slots__ = ('_epoch',)
    
    julianepoch = True
    """
    If True, the epoch is Julian, otherwise, Besselian
    """
    
    def __getstate__(self):
        return {'_epoch':self._epoch}
    
    def __setstate__(self,d):
        self._epoch = d['_epoch']
    
    def _getEpoch(self):
        return self._epoch
    def _setEpoch(self,val):
        if val is None:
            self._epoch = None
        else:
            if val == 'now':
                from ..obstools import jd_to_epoch
                val = jd_to_epoch(None,self.julianepoch)
            if not hasattr(self,'_epoch') or self._epoch is None:
                self._epoch = float(val)
            else:
                self.transformToEpoch(float(val))
    epoch = property(_getEpoch,_setEpoch,doc="""
    Epoch for this coordinate as a float. 
    
    Setting with the string 'now' will set the epoch to the time at the moment
    the command is executed.
    
    If set, this coordinate will be transformed to the new epoch, unless the
    current Epoch is None, in which case the epoch will be set with no
    transformation. If transformation is *not* desired, first set the epoch to
    None, and then set to the new epoch.
    
    Julian vs. Besselian is determined by the :attr:`julianepoch` attribute.
    """)
    
    def _getEpochstr(self):
        #format requires 2.6
        #return '{0}{1}'.format('J' if self.julianepoch else 'B',self._epoch)
        if self._epoch is None:
            return ''
        else:
            return '%s%s'%('J' if self.julianepoch else 'B',self._epoch)
    def _setEpochstr(self,val):
        self.epoch = val
    epochstr = property(_getEpochstr,_setEpochstr,doc="""
    A string representation of the epoch of this object with a J or B prefixed
    for julian or besselian epochs.
    """)
    
    def _getJdepoch(self):
        from ..obstools import epoch_to_jd
        return epoch_to_jd(self._epoch,self.julianepoch)
    def _setJdepoch(self,val):
        from ..obstools import jd_to_epoch
        self._epoch = jd_to_epoch(val,self.julianepoch)
    jdepoch = property(_getJdepoch,_setJdepoch,doc="""
    Julian Date of the epoch for this object.
    """)
    
    def _getMjdepoch(self):
        from ..obstools import epoch_to_jd
        return epoch_to_jd(self._epoch,self.julianepoch,mjd=True)
    def _setMjdepoch(self,val):
        from ..obstools import jd_to_epoch
        self._epoch = jd_to_epoch(val,self.julianepoch,mjd=True)
    mjdepoch = property(_getMjdepoch,_setMjdepoch,doc="""
    Modified Julian Date of the epoch for this object.
    """)
    
    @abstractmethod
    def transformToEpoch(self,newepoch):
        """
        Subclasses should implement this method to transform their coordinates
        to a new epoch. At the end of this method after the necessary data
        transformations are performed, subclasses should call
        ``EpochalCoordinates.transformToEpoch(newepoch)``.
        """
        self._epoch = newepoch

class RectangularCoordinates(CoordinateSystem):
    """
    Rectangular/Cartesian Coordinates in three dimensions. Coordinates are
    accessed via the attributes :attr:`x`, :attr:`y`, and :attr:`z`.
    """
    
    __slots__ = ('x','y','z')
    
    def __init__(self,x,y,z):
        #: x cartesian coordinate value
        self.x = x
        #: y cartesian coordinate value
        self.y = y
        #: z cartesian coordinate value
        self.z = z
        
    def __getstate__(self):
        #TODO: watch if this creates probelms by not being a dict
        return dict(x=self.x,y=self.y,z=self.z)
    
    def __setstate__(self,d):
        self.x = d['x']
        self.y = d['y']
        self.z = d['z']
        
    def __str__(self):
        return '%s: x=%f,y=%f,z=%f'%(self.__class__.__name__,self.x,self.y,self.z)
        
    def __add__(self,other):
        from copy import deepcopy
        if hasattr(other,'x') and hasattr(other,'y') and hasattr(other,'z'):
            new = deepcopy(self)
            new.x = self.x+other.x
            new.y = self.y+other.y
            new.z = self.z+other.z
            return new
        raise TypeError('Object of type %s does not have x,y, and z for operand +'%other.__class__)
        
    def __sub__(self,other):
        from copy import deepcopy
        if hasattr(other,'x') and hasattr(other,'y') and hasattr(other,'z'):
            new = deepcopy(self)
            new.x = self.x-other.x
            new.y = self.y-other.y
            new.z = self.z-other.z
            return new
        raise TypeError('Object of type %s does not have x,y, and z for operand -'%other.__class__)
        
    def _getLength(self):
        from math import sqrt
        return sqrt(self.x**2+self.y**2+self.z**2)
    def _setLength(self,val):
        scaling = val/self._getLength()
        self.x *= scaling
        self.y *= scaling
        self.z *= scaling        
    length = property(_getLength,_setLength,doc="""
    The Length of this coordinate's vector e.g. distance from the origin. If
    set, the direction will be preserved but the vector will be scaled to the
    provided value.""")
    
    
CartesianCoordinates = RectangularCoordinates
    

class _LatLongMeta(_CoosysMeta):
    def __init__(cls,name,bases,dct):
        _CoosysMeta.__init__(cls,name,bases,dct)
        if cls._longlatnames_[0] is not None:
            setattr(cls,cls._longlatnames_[0],cls.long)
            setattr(cls,cls._longlatnames_[0]+'err',cls.longerr)
        if cls._longlatnames_[1] is not None:
            setattr(cls,cls._longlatnames_[1],cls.lat)
            setattr(cls,cls._longlatnames_[1]+'err',cls.laterr)

class LatLongCoordinates(CoordinateSystem):
    """
    This object represents an angular location on a sphere as represented in
    spherical coordinates with a latitude and longitude, and optionally a
    distance. Subclasses specify details such as transformations or epochs.
    
    A :class:`LatLongCoordinate` system is designed to use the transformation
    type (see :meth:`CoordinateSystem.addTransType`) 'smatrix'. Thus, if the
    transformation between two :class:`LatLongCoordinate` subclasses can be
    represented as a unitary matrix operating on position vectors on the unit
    sphere, the transformation function can be written as::
    
        @CoordinateSystem.registerTransform(InLLCoords,OutLLCoords,transtype='smatrix')
        def transform(incoord):
            ... compute the elements of a 3x3 transformation matrix...
            return np.mat([[a,b,c],[d,e,f],[g,h,i]])
        
    *Subclassing*
    
    Subclasses of :class:`LatLongCoordinates` can have the class attribute
    :attr:`_longlatnames_` as a 2-tuple of strings (longname,latname), with
    names for the two coordinates, e.g. ('ra','dec'). They can also include the
    :attr:`_longrange_` attribute, which specifies the range of valid values for
    the longitude (latitude is always -90 to 90 degrees), or None to place no
    restriction. See :class:`CoordinateSystem` for additional subclassing
    information.
    
    
    """
    __slots__ = ('_lat','_long','_laterr','_longerr','_dpc')
    __metaclass__ = _LatLongMeta
    _longlatnames_ = ('longitude','latitude')
    _longrange_ = None
    
    
    def __init__(self,long=0,lat=0,longerr=None,laterr=None,distancepc=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        
        self._lat = AngularCoordinate(range=(-90,90,360))
        self._long = AngularCoordinate(range=self._longrange_)
        self.distancepc = distancepc
        
        if hasattr(lat,'lat') and hasattr(lat,'long'):
            if long is 0 and laterr is None and longerr is None:
                self.lat = lat.lat
                self.long = lat.long
                self.laterr = lat.laterr
                self.longerr = lat.longerr
            else:
                raise ValueError("can't provide a LatLongCoordinates as a constructor and set other values simultaneously")
        else:
            self.lat = lat
            self.long = long
            self.laterr = laterr
            self.longerr = longerr
        
    def __getstate__(self):
        return dict([(k,getattr(self,k)) for k in LatLongCoordinates.__slots__])
    
    def __setstate__(self,d):
        for k in LatLongCoordinates.__slots__:
            setattr(self,k,d[k])
            
    def _getDistancepc(self):
        if callable(self._dpc):
            return self._dpc()
        else:
            return self._dpc
    def _setDistancepc(self,val):
        if val is None:
            self._dpc = None
        elif callable(val):
            self._dpc = val
        else:
            try:
                self._dpc = (float(val[0]),float(val[1]))
            except (TypeError,IndexError):
                self._dpc = (float(val),0)
    distancepc = property(_getDistancepc,_setDistancepc,doc="""
    Parallax distance to object in parsecs, or None to assume infinity. Set as
    either a float, a 2-tuple (distance,distance_error), or a no-argument
    callable that returns such a tuple. Getter always returns 2-tuple or None.
    """)
    
    def _getDistanceau(self):
        from ..constants import aupercm,cmperpc
        auperpc = cmperpc * aupercm
        if self._dpc is None:
            return None
        elif callable(self._dpc):
            res = self._dpc()
        else:
            res = self._dpc
        return (res[0]*auperpc,res[1]*auperpc)
    def _setDistanceau(self,val):
        from ..constants import cmperau,pcpercm
        pcperau = cmperau * pcpercm
            
        if val is None:
            self._dpc = None
        elif callable(val):
            self._dpc = lambda: tuple((v*pcperau for v in val()))
        else:
            try:
                self._dpc = (float(val[0])*pcperau,float(val[1])*pcperau)
            except (TypeError,IndexError):
                self._dpc = (float(val)*pcperau,0)
    distanceau = property(_getDistanceau,_setDistanceau,doc="""
    Parallax distance to object in AU, or None to assume infinity. Set as either
    a float, a 2-tuple (distance,distance_error), or a no-argument callable that
    returns such a tuple. Getter always returns 2-tuple or None.
    """)
        
    def _getLat(self):
        return self._lat
    def _setLat(self,val):
        if isinstance(val,AngularCoordinate):
            rads = val.radians%_twopi
        else:
            rads = AngularCoordinate(val).radians%_twopi
        #fix for radian range
        if rads > 3*pi/2:
            rads -= _twopi
        elif rads > pi/2:
            rads = pi - rads
        self._lat.radians = rads
    lat = property(_getLat,_setLat,doc="""
    Latitude of this object as a :class:`AngularCoordinate` object.  May be set
    using any valid input form for :class:`AngularCoordinate`.
    """)
    
    def _getLong(self):
        return self._long
    def _setLong(self,val):
        if isinstance(val,AngularCoordinate):
            self._long.radians = val.radians%_twopi
        else:
            self._long.radians = AngularCoordinate(val).radians%_twopi
    long = property(_getLong,_setLong,doc="""
    Longitude of this object as a :class:`AngularCoordinate` object.  May be set
    using any valid input form for :class:`AngularCoordinate`.
    """)
    
    def _getLaterr(self):
        return self._laterr
    def _setLaterr(self,val):
        if val is None:
            self._laterr = None
        elif isinstance(val,AngularSeparation):
            self._laterr = val
        else:
            self._laterr = AngularSeparation(val)
    laterr = property(_getLaterr,_setLaterr,doc="""
    Latitude error for this object as a :class:`AngularSeparation` object. May
    be set using any valid input form for :class:`AngularSeparation`.
    """)
    
    def _getLongerr(self):
        return self._longerr
    def _setLongerr(self,val):
        if val is None:
            self._longerr = None
        elif isinstance(val,AngularSeparation):
            self._longerr = val
        else:
            self._longerr = AngularSeparation(val)
    longerr = property(_getLongerr,_setLongerr,doc="""
    Longitude error for this object as a :class:`AngularSeparation` object. May
    be set using any valid input form for :class:`AngularSeparation`.
    """)
    
    def __str__(self):
        lat,long = self._lat.d,self._long.d
        #lat,long = self._lat.getDmsStr(),self._long.getDmsStr()
        #this requires 2.6 - switch later maybe
        #return '{0}: {1[0]}={2},{1[1]}={3}'.format(self.__class__.__name__,self._longlatnames_,long,lat)
        return '%s: %s=%f,%s=%f'%(self.__class__.__name__,self._longlatnames_[0],long,self._longlatnames_[1],lat)
    
    def getCoordinateString(self,sep=' ',labels=False,canonical=False,hmslong=False):
        coords = []
        if hmslong:
            coords.append(self._long.getHmsStr(canonical=canonical,sign=False))
        else:
            coords.append(self._long.getDmsStr(canonical=canonical,sign=False))
        coords[-1] = self._longlatnames_[0]+'='+coords[-1]
        coords.append(self._lat.getDmsStr(canonical=canonical))
        coords[-1] = self._longlatnames_[1]+'='+coords[-1]
        return sep.join(coords)
    
    def __eq__(self,other):
        if hasattr(other,'lat') and hasattr(other,'long'):
            return self._lat==other.lat and self._long==other.long
        else:
            return False
        
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __sub__(self,other):        
        if isinstance(other,LatLongCoordinates) or (hasattr(other,'lat') and hasattr(other,'long')):
            from math import cos,degrees,acos,asin,sin,sqrt
            
            b1 = self._lat.radians
            b2 = other.lat.radians
            db = abs(b2 - b1)
            dl = abs(other.long.radians - self._long.radians)
            
            #haversin(theta) = (1-cos(theta))/2 = sin^2(theta/2)
            #has better numerical accuracy if sin for theta ~ 0, cos ~ pi/2
            haversin = lambda t:(1-cos(t))/2 if pi/4 < (t%pi) < 3*pi/4 else sin(t/2)**2
            
            hdb = haversin(db)
            hdl = haversin(dl)
            
            havsep = hdb + cos(b1)*cos(b2)*hdl
            #archaversin
            sep = acos(1 - 2*havsep) if 0.25 < havsep <= 0.75 else 2*asin(havsep**0.5)
                
            #straightforward definition without the tweaks using haversin - this
            #is in principal faster, but in practice it ends up only about
            #10% faster due to the other overhead
            #sep = acos(sin(b1)*sin(b2)+cos(b1)*cos(b2)*cos(dl))
            
            return AngularSeparation(degrees(sep))
            
#            #small angle version
#            from math import cos,degrees,sqrt
#            clat = cos((self._lat.radians+other._lat.radians)/2)
#            dcorrlong = (self._long.radians - other._long.radians)*clat
#            dlat = self._lat.radians-other._lat.radians
#            sep = AngularSeparation(degrees(sqrt(dlat*dlat+dcorrlong*dcorrlong)))
#            return sep
        else:
            raise ValueError("unsupported operand type(s) for -: '%s' and '%s'"%(self.__class__,other.__class__))
        
    @staticmethod
    @CoordinateSystem.addTransType
    def _smatrix(m,coord,tocls):
        newcoord = tocls()
        newcoord.lat = coord._lat
        newcoord.laterr = coord._laterr
        newcoord.long = coord._long
        newcoord.longerr = coord._longerr
        newcoord.matrixRotate(m)
        return newcoord
        
    def matrixRotate(self,matrix,apply=True,fixrange=True,unitarycheck=False):
        """
        Applies the supplied  unitary rotation matrix to these coordinates. 
        
        :param matrix: the transformation matrix in cartesian coordinates
        :type matrix: a 3x3 :class:`numpy.matrix`
        :param apply: 
            If True, the transform will be applied inplace to the coordinates
            for this object
        :type apply: boolean
        :param fixrange: 
            If True the latitude is autmoatically fixed to be on (-pi/2,pi/2) 
            and the longitude is on (0,2pi).  Otherwise the raw coordinate is
            output.
        :type fixrange: boolean
        :param unitarycheck: 
            If True and the matrix is not unitary, a ValueError will be raised.
            Otherwise no check is performed.
        :type unitarycheck: boolean
        
        :returns: 
            (lat,long) as decimal radians after the transformation matrix is
            applied or (lat,long,laterr,longerr) if errors are nonzero
        """
        #for single values, math module is much faster than numpy 
        from math import sin,cos,atan2,sqrt
        
        m = np.asmatrix(matrix)
        
        if unitarycheck:
            mdagger = m.H
            rtol = 1e-5 if unitarycheck is True else unitarycheck
            if not np.allclose(mdagger*m,m*mdagger,rtol):
                raise ValueError('matrix not unitary')
        
        lat = self.lat.radians
        long = self.long.radians
        laterr = 0 if self.laterr is None else self.laterr.radians
        longerr = 0 if self.longerr is None else self.longerr.radians    
        
        sb = sin(lat)
        cb = cos(lat)
        sl = sin(long)
        cl = cos(long)
        
        #spherical w/ r=1 > cartesian
        x = cb*cl
        y = cb*sl
        z = sb
        
        #do transform
        v = np.matrix((x,y,z)).T
        xp,yp,zp = (m*v).A1
        
        #cartesian > spherical
        sp = sqrt(xp*xp+yp*yp) #cylindrical radius
        latp = atan2(zp,sp)
        longp = atan2(yp,xp)
        
        #propogate errors if they are present
        
        if laterr != 0 or longerr != 0:
            #TODO: check formulae
            #all of these are first order taylor expansions about the value
            dx = sqrt((laterr*sb*cl)**2+(longerr*cb*sl)**2)
            dy = sqrt((laterr*sb*sl)**2+(longerr*cb*cl)**2)
            dz = abs(laterr*cb)
            
            dv = np.matrix((dx,dy,dz))
            dxp,dyp,dzp = np.sqrt(np.power(m,2)*np.power(dv,2))
            
            #intermediate variables for dlatp - each of the partial derivatives
            chi = 1/(1+(zp/sp)**2) 
            #common factor chi not included below
            dbdx = x*z*sp**-3
            dbdy = y*z*sp**-3
            dbdz = 1/sp
            
            dlatp = chi*sqrt((dxp*dbdx)**2 + (dyp*dbdy)**2 + (dzp*dbdz)**2)
            dlongp = sqrt((dxp*yp*xp**-2)**2 + (dyp/xp)**2)/(1 + (yp/xp)**2) #indep of z
            
        else:
            laterr = None
            
        if fixrange:
            ao = (latp+_pio2)/_twopi
            latp = _twopi*abs((ao-np.floor(ao+0.5)))-_pio2
            longp = longp % _twopi
        
        if apply:
            self.lat.radians = latp
            self.long.radians = longp
            if laterr is not None:
                self.laterr.radians = dlatp
                self.longerr.radians = dlongp
        
        if laterr is None:
            return latp,longp
        else:
            return latp,longp,dlatp,dlongp
    
    def convert(self,tosys,optimize=_convertoptimizedefault):
        """
        Converts the coordinate system from it's current system to a new
        :class:`CoordinateSystem` object possibly with optimizations for
        matrix-based transformation of :class:`LatLongCoordinates` objects.
        
        .. warning::
            The transformation optimizations used if `optimize` is True are only
            correct if the conversion matricies are independent of the
            coordinate values themselves (e.g. are linear and
            epoch-independent). Examples that are good for optimization include
            FK4->FK5->ICRS (all are epoch-independent).
            
            In the future, this system will be smarter about
            knowing when it is safe to optimize, but for now, only use it if
            you're sure it will work correctly.
        
        :param tosys: 
            The new coordinate system class. Should be a subclass of
            :class:`CoordinateSystem` .
        :param bool optimize: 
            If True, speed up the transformation by composing matricies where
            possible. If False, the standard transformation is performed.
        :returns: A new object of a class determined by `tosys`
        
        
        :except: raises NotImplementedError if converters are not present
        """
        if tosys is self.__class__:
            return self
        
        if optimize:
            cache = CoordinateSystem._transformcache['smatrix']
            
            #add this object to the cache if its missing
            if self.__class__ not in cache:
                cache[self.__class__] = {}
                
            if tosys not in cache[self.__class__]:
                convs = CoordinateSystem.getTransformPath(self.__class__,tosys)
                
                if callable(convs): #direct transform
                    convs = [convs]
                else:
                    convclasses = convs
                    convfuncs = [CoordinateSystem._converters[c1][c2] for c1,c2 in zip(convclasses[:-1],convclasses[1:])]
                    convs = []
                    
                    #now we populate convs with converter functions that are 
                    #either multplied-together matricies if they are smatrix
                    #converters or the actual converter function otherwise
                    combinedmatrix = None
                    lastcls = None
                    for cls,cfunc in zip(convclasses[:-1],convfuncs):
                        #note that cls here is the *previous* conversion's end 
                        #class/current conversion's start class...
                        if cfunc.transtype=='smatrix':
                            mt = cfunc.basetrans(self)
                            
                            if hasattr(mt,'nocache') and mt.nocache:
                                cache = None
                            
                            if combinedmatrix is None:
                                combinedmatrix = mt
                            else:
                                combinedmatrix = mt * combinedmatrix
                        else:
                            if combinedmatrix is None:
                                convs.append(cfunc)
                            else:
                                convs.append(_OptimizerSmatrixer(combinedmatrix,cls))
                                convs.append(cfunc)
                                
                    if combinedmatrix is not None:
                        convs.append(_OptimizerSmatrixer(combinedmatrix,convclasses[-1]))
                
                #now cache this transform for future use unless it was banned above
                if cache is not None:
                    cache[self.__class__][tosys] = convs
                
            else:
                convs = cache[self.__class__][tosys]
            
            #now actually do the transforms
            coord = self
            for conv in convs:
                coord = conv(coord)
            return coord
        
        else:
            return CoordinateSystem.convert(self,tosys)
        
class _OptimizerSmatrixer(object):
    """
    Used internally to do the optimization of :meth`LatLongCoordinates.convert`
    """
    transtype = 'smatrix'
    def __init__(self,combinedmatrix,tocls):
        self.combinedmatrix = combinedmatrix
        self.tocls = tocls
    def __call__(self,coord):
        return LatLongCoordinates._smatrix(self.combinedmatrix,coord,self.tocls)
    def basetrans(self,coords):
        return self.combinedmatrix
        
    
class EpochalLatLongCoordinates(LatLongCoordinates,EpochalCoordinates):
    """
    A Coordinate system where the coordinates change as a function of time.
    The origin and orientation of some coordinate systems are tied to the motion
    of the Earth and Solar System and hence most be updated as time passes.  
    
    In general this only accounts for epoch-related changes in its own
    coordinate system. If (for example) one has a :class:`ITRSCoordinates`
    coordinate, changing the epoch only updates for polar motion. To properly
    update all epoch-related such as precession/nutation and earth rotation, the
    coordinate should be transformed to :class:`ICRSCoordinates` , update the
    epoch, and transform back to :class:`TIRSCoordinates` .
    
    """
    __slots__ = tuple()
    julianepoch = True 
    
    def __init__(self,long=0,lat=0,longerr=None,laterr=None,epoch=None,distancepc=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        LatLongCoordinates.__init__(self,long,lat,longerr,laterr,distancepc)
        self._epoch = epoch
        
    
    def __getstate__(self):
        d = LatLongCoordinates.__getstate__(self)
        d.update(EpochalCoordinates.__getstate__(self))
        return d
    
    def __setstate__(self,d):
        LatLongCoordinates.__setstate__(self,d)
        EpochalCoordinates.__setstate__(self,d)
    
    @add_docs(LatLongCoordinates.convert)
    def convert(self,tosys,optimize=_convertoptimizedefault):
        ''
        res = LatLongCoordinates.convert(self,tosys,optimize)
        if issubclass(tosys,EpochalLatLongCoordinates):
            res._epoch = self._epoch
        return res
    convert.__doc__ = LatLongCoordinates.convert.__doc__
        

class EquatorialCoordinatesBase(EpochalLatLongCoordinates):
    """
    This object represents an angular location on the unit sphere, specified in
    right ascension and declination.  Some of the subclasses are not strictly 
    speaking equatorial, but they are close, or are tied to the equatorial 
    position of a particular epoch.
        
    This is a superclass for all other Equatorial Coordinate systems -
    particular reference systems implement the :meth:`transformToEpoch` method.
    See the docstring for :class:`EpochalLatLongCoordinates` for subclassing
    suggestions.
    """
    
    __slots__ = tuple()
    _longlatnames_ = ('ra','dec')
    _longrange_ = (0,360)
    
        
    def __str__(self):
        rastr = self.ra.getHmsStr(canonical=True)
        decstr = self.dec.getDmsStr(canonical=True)
        #2.6 required for format
        #return '{3}: {0} {1} ({2})'.format(rastr,decstr,self.epoch,self.__class__.__name__)
        return  '%s: %s %s %s'%(self.__class__.__name__,rastr,decstr,self.epochstr)
    
    
    def __init__(self,*args,**kwargs):
        """
        Input for equatorial coordinates. Can follow any of the following forms:
        
        * EquatorialCoordinatesBase()
        * EquatorialCoordinatesBase(:class:`EquatorialCoordinatesBase`)
        * EquatorialCoordinatesBase('rastr decstr')
        * EquatorialCoordinatesBase((ra,dec))
        * EquatorialCoordinatesBase(ra,dec)
        * EquatorialCoordinatesBase(ra,fdec,raerr,decerr)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch,distancepc)
        
        Note that the default epoch is 2000 if not otherwise specified.  To 
        disable epoch tranformations, set the epoch to None.  If scalar values
        are provided, they are assummed to be degrees.
        
        """
        posargs = {}
        if len(args) == 0:
            pass
        if len(args) == 1:
            if isinstance(args[0],EquatorialCoordinatesBase):
                EpochalLatLongCoordinates.__init__(self, args[0].ra, args[0].dec,
                                                   args[0].raerr, args[0].decerr,
                                                   args[0].epoch, args[0].distancepc)
                return
            elif isinstance(args[0],basestring):
                sargs = args[0].split()
                posargs['ra'] = sargs[0]
                posargs['dec'] = sargs[1]
            else:
                posargs['ra'],posargs['dec'] = args[0]
        elif len(args) == 2:
            posargs['ra'] = args[0]
            posargs['dec'] = args[1]
        elif len(args) == 4:
            posargs['ra'] = args[0]
            posargs['dec'] = args[1]
            posargs['raerr'] = args[2]
            posargs['decerr'] = args[3]
        elif len(args) == 5:
            posargs['ra'] = args[0]
            posargs['dec'] = args[1]
            posargs['raerr'] = args[2]
            posargs['decerr'] = args[3]
            posargs['epoch'] = args[4]
        elif len(args) == 6:
            posargs['ra'] = args[0]
            posargs['dec'] = args[1]
            posargs['raerr'] = args[2]
            posargs['decerr'] = args[3]
            posargs['epoch'] = args[4]
            posargs['distancepc'] = args[5]
        
        for k,v in posargs.iteritems():
            if k in kwargs:
                raise ValueError('got multiple values for argument '+k)
            kwargs[k] = v
        
        kwargs.setdefault('ra',0)
        kwargs.setdefault('dec',0)
        kwargs.setdefault('raerr',None)
        kwargs.setdefault('decerr',None)
        kwargs.setdefault('epoch',2000)
        
        EpochalLatLongCoordinates.__init__(self,kwargs['ra'],kwargs['dec'],kwargs['raerr'],kwargs['decerr'])
        if 'epoch' in kwargs:
            self.epoch = kwargs['epoch']
        if 'distancepc' in kwargs:
            if 'distanceau' in kwargs:
                raise TypeError("can't specify distance in both pc and au")
            self.distancepc = kwargs['distancepc']
        elif 'distanceau' in kwargs:
            self.distanceau = kwargs['distanceau']
        else:
            self._dpc = None
    
    @add_docs(EpochalLatLongCoordinates.convert)
    def convert(self,tosys,optimize=_convertoptimizedefault):
        ''
        res = EpochalLatLongCoordinates.convert(self,tosys,optimize)
        if self._dpc is not None:
            res._dpc = self._dpc
        return res
    convert.__doc__ = EpochalLatLongCoordinates.convert.__doc__
    
class ICRSCoordinates(EquatorialCoordinatesBase):
    """
    Equatorial Coordinates tied to the International Celestial Reference System
    (ICRS). Strictly speaking this is not an Equatorial system, as it is an
    inertial frame that only aligns with Earth's equator at J2000, but it is
    nearly an equatorial system at J2000.
        
    .. note::
        Technically, this is actually the Barycentric Celestial Referense System
        (BCRS), distringuished from ICRS by having acounted for space motion. In
        astropysics, instead, space motion is already accounted for by using
        :class:`astropysics.coords.ephems.EphemerisObject` objects, which yield
        coordinates (often :class:`ICRSCoordinates`) at the epoch of
        observation.
        
    .. warning:: 
        Abberation of starlight is not yet implemented in transformations
        to/from ICRS.
        
    """
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        """
        ICRS is an inertial frame, so no transformation is necessary
        """
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    #Conversions to/from ICRS are in the other coordinate systems

    def __makeFrameBias():
        from ..utils import rotation_matrix
        
        #sing USNO circular 179 for frame bias -- all in milliarcsec
        da0 = -14.6
        xi0 = -16.6170
        eta0 = -6.8192
        #mas->degrees
        return rotation_matrix(-eta0/3600000,axis='x') *\
               rotation_matrix(xi0/3600000,axis='y') *\
               rotation_matrix(da0/3600000,axis='z')
    frameBiasJ2000 = __makeFrameBias()
    """
    Frame bias matrix such that vJ2000 = B*vICRS . 
    """
    #make it a static method just for clarity even though it isn't really visible
    __makeFrameBias = staticmethod(__makeFrameBias)
    
class RectangularICRSCoordinates(RectangularCoordinates,EpochalCoordinates):
    """
    Rectangular coordinates aligned to the ICRS with origin at the solar system
    barycenter. The positive z-axis points to the north celestial pole and the
    positive x-axis  is along with the (0,0) point of the equatorial ICRS. 
    
    .. note::
        Units for the coordinates are specified via the :attr:`unit` attribute.
        When converting *from* :class:`ICRSCoordinates`, distances default to AU
        if less than 1000 AU, otherwise, pc. If a distance is not present, the
        default distance is 1 (unspecified unit).
    
    """
    __slots__ = tuple()
    julianepoch = True
    
    def __init__(self,x,y,z,epoch=None,unit='pc'):
        RectangularCoordinates.__init__(self,x,y,z)
        self._epoch = epoch
        self._unit = None
        self.unit = unit
        
    def __getstate__(self):
        d = RectangularCoordinates.__getstate__(self)
        d.update(EpochalCoordinates.__getstate__(self))
        return d
    
    def __setstate__(self,d):
        RectangularCoordinates.__setstate__(self,d)
        EpochalCoordinates.__setstate__(self,d)
        
    def __str__(self):
        if self.epoch is None:
            epochstr = ''
        else:
            epochstr = ' ('+self.epochstr+')'
        return RectangularCoordinates.__str__(self) + epochstr
        
    def transformToEpoch(self,newepoch):
        EpochalCoordinates.transformToEpoch(self,newepoch)
        
    def _getUnit(self):
        return self._unit
    def _setUnit(self,val):
        from ..constants import auperpc
        
        if val is None:
            self._unit = None
        elif self._unit is None and (val in ('au','pc')):
            self._unit = val
        elif val=='au':
            if self._unit == 'pc':
                self.x *= auperpc
                self.y *= auperpc
                self.z *= auperpc
                self._unit = val
        elif val == 'pc':
            if self._unit == 'au':
                self.x /= auperpc
                self.y /= auperpc
                self.z /= auperpc
                self._unit = val
        else:
            raise ValueError("unit must be 'au' or 'pc' - got %s"%val)
        
    unit = property(_getUnit,_setUnit,doc="""The unit for these coordinates. 
        Must be 'au', 'pc', or None - setting to anything else will raise a
        :exc:`ValueError`. If not None, setting to a new unit will convert the
        values from AU to pc or vice versa.
        """)
    
    
    @CoordinateSystem.registerTransform('self',ICRSCoordinates)
    def _toICRS(ric):
        from math import asin,atan2,degrees
        
        x,y,z = ric.x,ric.y,ric.z        
        r = (x*x+y*y+z*z)**0.5
        dec = degrees(asin(z/r))
        ra = degrees(atan2(y,x))
        
        if ric.unit is None:
            return ICRSCoordinates(ra,dec,epoch=ric.epoch)
        elif ric.unit == 'pc':
            return ICRSCoordinates(ra,dec,distancepc=r,epoch=ric.epoch)
        elif ric.unit == 'au':
            return ICRSCoordinates(ra,dec,distanceau=r,epoch=ric.epoch)
        else:
            raise NotImplementedError('Unrecognized unit %s in RectICRS->ICRS'%ric.unit)
    
    @CoordinateSystem.registerTransform(ICRSCoordinates,'self')    
    def _fromICRS(ic):
        from math import sin,cos
        
        ra,dec = ic.ra.r,ic.dec.r
        if ic.distanceau is None:
            r = 1
            unit = None
        elif ic.distanceau>1000:
            r = ic.distancepc[0]
            unit = 'pc'
        else:
            r = ic.distanceau[0]
            unit = 'au'
        x = r*cos(ra)*cos(dec)
        y = r*sin(ra)*cos(dec)
        z = r*sin(dec)
        
        return RectangularICRSCoordinates(x,y,z,ic.epoch,unit)
    
    
class GCRSCoordinates(EquatorialCoordinatesBase):
    """
    Geocentric Celestial Reference System equatorial coordinates. The
    orientation of this coordinates is fixed to the ICRS orientation, but with
    origin at the earth geocenter.
        
    .. warning:: 
        Abberation of starlight not yet included in transforms.
        
    """
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        """
        Transforms from the current epoch to a new epoch by converting to ICRS
        and back again in the new epoch.  
        """
        c1 = self.convert(ICRSCoordinates)
        c1.epoch = newepoch
        c2 = c1.convert(GCRSCoordinates)
        
        self._lat = c2._lat
        self._long = c2._long
        self._laterr = c2._laterr
        self._longerr = c2._longerr
        self._dpc = c2._dpc
            
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
#TODO:implement direct spherical transformations if needed/wanted for precision
#    @CoordinateSystem.registerTransform('self',ICRSCoordinates,transtype='smatrix')
#    def _toICRS(gc):
#        return np.eye(3).view(np.matrix)
    
#    @CoordinateSystem.registerTransform(ICRSCoordinates,'self',transtype='smatrix')    
#    def _fromICRS(ic):
#        return np.eye(3).view(np.matrix)
    
class RectangularGCRSCoordinates(RectangularCoordinates,EpochalCoordinates):
    """
    Rectangular coordinates aligned to the GCRS with origin at the Earth
    geocenter. The positive z-axis points to the north celestial pole and the
    positive x-axis points down the (0,0) point of the equatorial GCRS (and
    thus, also ICRS). 
    
    The coordinates at this location are actually an "Astrometric place" - the
    actual location relative to the geocenter. This is disctinct from
    :class:`GCRSCoordinates` in that :class:`GCRSCoordinates` *includes*
    aberration and light deflection, while :class:`RectangularGCRSCoordinates`
    does not.
    
    .. note::
        Units for the coordinates are specified via the :attr:`unit` attribute.
        When converting *from* :class:`GCRSCoordinates`, distances default to AU
        if less than 1000 AU, otherwise, pc. If a distance is not present, the
        default distance is 1 (unspecified unit).
    """
    __slots__ = tuple()
    julianepoch = True
    
    def __init__(self,x,y,z,epoch=None,unit='pc'):
        RectangularCoordinates.__init__(self,x,y,z)
        self._epoch = epoch
        self._unit = None
        self.unit = unit
        
    def __getstate__(self):
        d = RectangularCoordinates.__getstate__(self)
        d.update(EpochalCoordinates.__getstate__(self))
        return d
    
    def __setstate__(self,d):
        RectangularCoordinates.__setstate__(self,d)
        EpochalCoordinates.__setstate__(self,d)
        
    def __str__(self):
        if self.epoch is None:
            epochstr = ''
        else:
            epochstr = ' ('+self.epochstr+')'
        return RectangularCoordinates.__str__(self) + epochstr
        
    def transformToEpoch(self,newepoch):
        EpochalCoordinates.transformToEpoch(newepoch)
        
    def _getUnit(self):
        return self._unit
    def _setUnit(self,val):
        from ..constants import auperpc
        
        if val is None:
            self._unit = None
        elif self._unit is None and (val in ('au','pc')):
            self._unit = val
        elif val=='au':
            if self._unit == 'pc':
                self.x *= auperpc
                self.y *= auperpc
                self.z *= auperpc
                self._unit = val
        elif val == 'pc':
            if self._unit == 'au':
                self.x /= auperpc
                self.y /= auperpc
                self.z /= auperpc
                self._unit = val
        else:
            raise ValueError("unit must be 'au' or 'pc' - got %s"%val)
        
    unit = property(_getUnit,_setUnit,doc="""The unit for these coordinates. 
        Must be 'au', 'pc', or None - setting to anything else will raise a
        :exc:`ValueError`. If not None, setting to a new unit will convert the
        values from AU to pc or vice versa.
        """)
    
    
    @CoordinateSystem.registerTransform('self',GCRSCoordinates)
    def _toGCRS(rgc):
        from math import asin,atan2,degrees
        #TODO:implement aberration and light deflection
        
        x,y,z = rgc.x,rgc.y,rgc.z        
        r = (x*x+y*y+z*z)**0.5
        dec = degrees(asin(z/r))
        ra = degrees(atan2(y,x))
        
        if rgc.unit is None:
            return GCRSCoordinates(ra,dec,epoch=rgc.epoch)
        elif rgc.unit == 'pc':
            return GCRSCoordinates(ra,dec,distancepc=r,epoch=rgc.epoch)
        elif rgc.unit == 'au':
            return GCRSCoordinates(ra,dec,distanceau=r,epoch=rgc.epoch)
        else:
            raise NotImplementedError('Unrecognized unit %s in RectGCRS->GCRS'%rgc.unit)
    
    @CoordinateSystem.registerTransform(GCRSCoordinates,'self')    
    def _fromGCRS(gc):
        from math import sin,cos
        #TODO:implement aberration and light deflection
        
        ra,dec = gc.ra.r,gc.dec.r
        if gc.distanceau is None:
            r = 1
            unit = None
        elif gc.distanceau>1000:
            r = gc.distancepc[0]
            unit = 'pc'
        else:
            r = gc.distanceau[0]
            unit = 'au'
        x = r*cos(ra)*cos(dec)
        y = r*sin(ra)*cos(dec)
        z = r*sin(dec)
        
        return RectangularGCRSCoordinates(x,y,z,gc.epoch,unit)
    
    @CoordinateSystem.registerTransform('self',RectangularICRSCoordinates)
    def _toRectICRS(rgc):
        from .ephems import earth_pos_vel
        from ..obstools import epoch_to_jd
        from ..constants import auperpc
        
        x = rgc.x
        y = rgc.y
        z = rgc.z
        unit = rgc.unit
        epoch = rgc.epoch
        if epoch is None:
            raise ValueError('cannot transform GCRS to ICRS without an epoch')
          
        if unit is None: #infitiely far, so no corrections
            return RectangularICRSCoordinates(x,y,z,epoch,unit=None)
        else: #do parallax correction
            xe,ye,ze = earth_pos_vel(epoch_to_jd(epoch),True)[0]
            
            if unit == 'au':
                xp = x - xe
                yp = y - ye
                zp = z - ze
            elif unit == 'pc':
                xp = x - xe/auperpc
                yp = y - ye/auperpc
                zp = z - ze/auperpc
            else:
                raise NotImplementedError('Unit %s not supported by GCRS->ICRS'%unit)
            
            return RectangularICRSCoordinates(xp,yp,zp,epoch,unit=unit)
    
    @CoordinateSystem.registerTransform(RectangularICRSCoordinates,'self')    
    def _fromRectICRS(ric):
        from .ephems import earth_pos_vel
        from ..obstools import epoch_to_jd
        from ..constants import auperpc
        
        x = ric.x
        y = ric.y
        z = ric.z
        unit = ric.unit
        epoch = ric.epoch
        
        if epoch is None:
            raise ValueError('cannot transform ICRS to GCRS without an epoch')
        
        if unit is None: #infitiely far, so no corrections
            return RectangularGCRSCoordinates(x,y,z,epoch,unit=None)
        else: #do parallax correction
            xe,ye,ze = earth_pos_vel(epoch_to_jd(epoch),True)[0]
            
            if unit == 'au':
                xp = x - xe
                yp = y - ye
                zp = z - ze
            elif unit == 'pc':
                xp = x - xe/auperpc
                yp = y - ye/auperpc
                zp = z - ze/auperpc
            else:
                raise NotImplementedError('Unit %s not supported by ICRS->GCRS'%unit)
            
            return RectangularGCRSCoordinates(xp,yp,zp,epoch,unit=unit)
    
    
def _precession_matrix_J2000_Capitaine(epoch):
        """
        Computes the precession matrix from J2000 to the given Julian Epoch.
        Expression from from Capitaine et al. 2003 as written in the USNO
        Circular 179.  This should match the IAU 2006 standard from SOFA 
        (although this has not yet been tested)
        """
        from ..utils import rotation_matrix
        
        T = (epoch-2000.0)/100.0
        #from USNO circular
        pzeta = (-0.0000003173,-0.000005971,0.01801828,0.2988499,2306.083227,2.650545)
        pz = (-0.0000002904,-0.000028596,0.01826837,1.0927348,2306.077181,-2.650545)
        ptheta = (-0.0000001274,-0.000007089,-0.04182264,-0.4294934,2004.191903,0)
        zeta = np.polyval(pzeta,T)/3600.0
        z = np.polyval(pz,T)/3600.0
        theta = np.polyval(ptheta,T)/3600.0
        
        return rotation_matrix(-z,'z') *\
               rotation_matrix(theta,'y') *\
               rotation_matrix(-zeta,'z')
               
               
def _load_nutation_data(datafn,seriestype):
    """
    Loads nutation series from saved data files.
    
    Seriestype can be 'lunisolar' or 'planetary'
    """
    from ..utils.io import get_package_data
    
    if seriestype == 'lunisolar':
        dtypes = [('nl',int),
                  ('nlp',int),
                  ('nF',int),
                  ('nD',int),
                  ('nOm',int),
                  ('ps',float),
                  ('pst',float),
                  ('pc',float),
                  ('ec',float),
                  ('ect',float),
                  ('es',float)]
    elif seriestype == 'planetary':
        dtypes = [('nl',int),
                  ('nF',int),
                  ('nD',int),
                  ('nOm',int),
                  ('nme',int),
                  ('nve',int),
                  ('nea',int),
                  ('nma',int),
                  ('nju',int),
                  ('nsa',int),
                  ('nur',int),
                  ('nne',int),
                  ('npa',int),  
                  ('sp',int),
                  ('cp',int),
                  ('se',int),
                  ('ce',int)]
    else:
        raise ValueError('requested invalid nutation series type')
    
    lines = [l for l in get_package_data(datafn).split('\n') if not l.startswith('#') if not l.strip()=='']
    
    lists = [[] for n in dtypes]
    for l in lines:
        for i,e in enumerate(l.split(' ')):
            lists[i].append(dtypes[i][1](e))
    return np.rec.fromarrays(lists,names=[e[0] for e in dtypes])

_nut_data_00a_ls = _load_nutation_data('iau00a_nutation_ls.tab','lunisolar')
_nut_data_00a_pl = _load_nutation_data('iau00a_nutation_pl.tab','planetary')
def _nutation_components20062000A(epoch):
    """
    :returns: eps,dpsi,deps in radians
    """
    from ..obstools import epoch_to_jd
    from .funcs import obliquity
    
    epsa = obliquity(epoch_to_jd(epoch),2006)
    
    raise NotImplementedError('2006/2000A nutation model not implemented')

    return epsa,dpsi,deps


    
_nut_data_00b = _load_nutation_data('iau00b_nutation.tab','lunisolar')
def _nutation_components2000B(intime,asepoch=True):
    """
    :param intime: time to compute the nutation components as a JD or epoch
    :type intime: scalar
    :param asepoch: if True, `intime` is interpreted as an epoch, otherwise JD
    :type asepoch: bool
    
    :returns: eps,dpsi,deps in radians
    """
    from ..constants import asecperrad
    from ..obstools import epoch_to_jd,jd2000
    from .funcs import obliquity
    
    if asepoch:
        jd = epoch_to_jd(intime)
    else:
        jd = intime
    epsa = np.radians(obliquity(jd,2000))
    t = (jd-jd2000)/36525
    
    #Fundamental (Delaunay) arguments from Simon et al. (1994) via SOFA
    #Mean anomaly of moon
    el = ((485868.249036 + 1717915923.2178*t)%1296000)/asecperrad
    #Mean anomaly of sun
    elp = ((1287104.79305 + 129596581.0481*t)%1296000)/asecperrad
    #Mean argument of the latitude of Moon
    F = ((335779.526232 + 1739527262.8478*t)%1296000)/asecperrad
    #Mean elongation of the Moon from Sun
    D = ((1072260.70369 + 1602961601.2090*t)%1296000)/asecperrad
    #Mean longitude of the ascending node of Moon
    Om = ((450160.398036 + -6962890.5431*t)%1296000)/asecperrad
    
    #compute nutation series using array loaded from data directory
    dat = _nut_data_00b
    arg = dat.nl*el + dat.nlp*elp + dat.nF*F + dat.nD*D + dat.nOm*Om
    sarg = np.sin(arg)
    carg = np.cos(arg)
    
    p1uasecperrad = asecperrad*1e7 #0.1 microasrcsecperrad
    dpsils = np.sum((dat.ps + dat.pst*t)*sarg + dat.pc*carg)/p1uasecperrad
    depsls = np.sum((dat.ec + dat.ect*t)*carg + dat.es*sarg)/p1uasecperrad
    #fixed offset in place of planetary tersm
    masecperrad = asecperrad*1e3 #milliarcsec per rad
    dpsipl = -0.135/masecperrad
    depspl =  0.388/masecperrad
    
    return epsa,dpsils+dpsipl,depsls+depspl #all in radians
               
def _nutation_matrix(epoch):
    """
    Nutation matrix generated from nutation components.
    
    Matrix converts from mean coordinate to true coordinate as
    r_true = M * r_mean
    """
    from ..utils import rotation_matrix
    
    #TODO: implement higher precision 2006/2000A model if requested/needed
    epsa,dpsi,deps = _nutation_components2000B(epoch) #all in radians
    
    return rotation_matrix(-(epsa + deps),'x',False) *\
           rotation_matrix(-dpsi,'z',False) *\
           rotation_matrix(epsa,'x',False)
           

def _load_CIO_locator_data(datafn):
    """
    Loads CIO locator series terms from saved data files.
    
    returns polycoeffs,termsarr (starting with 0th)
    """
    from ..utils.io import get_package_data
    
    lines = [l for l in get_package_data(datafn).split('\n') if not l.startswith('#') if not l.strip()=='']
    coeffs = []
    sincs = []
    coscs = []
    orders = []
    
    inorder = False
    for l in lines:
        if 'Polynomial coefficients:' in l:
            polys = l.replace('Polynomial coefficients:','').split(',')
            polys = np.array(polys,dtype=float)
        elif 'order' in l:
            if inorder:
                orders.append((np.array(coeffs,dtype=int),
                               np.array(sincs,dtype=float),
                               np.array(coscs,dtype=float)))
                coeffs = []
                sincs = []
                coscs = []
            inorder = True
        elif inorder:
            ls = l.split()
            coeffs.append(ls[:8])
            sincs.append(ls[8])
            coscs.append(ls[9])
    if inorder:
        orders.append((np.array(coeffs,dtype=int),
                       np.array(sincs,dtype=float),
                       np.array(coscs,dtype=float)))
        
    return polys,orders
_CIO_locator_data = _load_CIO_locator_data('iau00_cio_locator.tab')  


class CIRSCoordinates(EquatorialCoordinatesBase):
    """
    Represents an object as equatorial coordinates in the Celestial Intermediate
    Reference System. This is the post-2000 IAU system for equatorial
    coordinates that seperates the coordinate system from the dynamically
    complicated and somewhat imprecisely-defined ideas of the equinox and
    ecliptic. This system's fundamental plane is the equator of the Celestial
    Intermediate Pole (CIP) and the origin of RA is at the Celestial
    Intermediate Origin (CIO).
    
    Changes to the :attr:`epoch` will result in the coordinates being updated
    for precession nutation. Nutation currently uses the IAU 2000B model that
    should be good to ~1 mas. If aberration or annual parallax corrections are
    necessary, convert to :class:`ICRSCoordinates`, change the epoch, and then
    convert back to :class:`CIRSCoordinates`.
    
    To convert from these coordinates to :class:`HorizontalCoordinates`
    appropriate for observed coordinates, site information is necessary. Hence,
    the transformations from equatorial to horizontal coordinates are performed
    by the :class:`~astropysics.obstools.Site` class in the
    :mod:`~astropysics.obstools` module, and attempting to directly convert will
    raise an :exc:`TypeError`.
    """
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch using the
        IAU 2000 precessions from Capitaine, N. et al. 2003 as written in the
        USNO Circular 179.
        """
        if self.epoch is not None and newepoch is not None:
            M = self._CMatrix(self.epoch).T
            Mn = self._CMatrix(newepoch)
            
            self.matrixRotate(Mn*M)
            
        #this sets the epoch
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    @staticmethod    
    def _CMatrix(epoch):
        """
        The GCRS->CIRS transformation matrix
        """
        B = ICRSCoordinates.frameBiasJ2000
        if epoch is None:
            return B
        else:
            from math import sin,cos,atan,atan2,sqrt
            from ..utils import rotation_matrix
            
            P = _precession_matrix_J2000_Capitaine(epoch)
            N = _nutation_matrix(epoch)
            
            
            x,y,z = (N*P*B).A[2]
            xsq,ysq = x**2,y**2
            bz = 1/(1+z)
            s = CIRSCoordinates._CIOLocator(epoch)
                   
            #matrix components - see Circular 179 or IERS Conventions 2003
            a,b,c = 1-bz*xsq , -bz*x*y , -x
            d,e,f = -bz*x*y , 1 - bz*ysq , -y
            g,h,i = x , y , 1 - bz*(xsq+ysq)
            
            #return rotation_matrix(-s,'z',defrees=False)*np.mat([[a,b,c],
            #                                                     [d,e,f],
            #                                                     [g,h,i]]) 
            
            si = sin(s)
            co = cos(s)
            
            M = [[a*co - d*si,b*co - e*si,c*co - f*si],
                 [a*si + d*co,b*si + e*co,c*si + f*co],
                 [     g,          h,          i     ]]     
            return np.mat(M)
        
#            #SOFA implementation using spherical angles - numerically identical
#            r2 = x*x + y*y
#            e = atan2(y,x) if r2 != 0 else 0
#            d = atan(sqrt(r2/(1-r2)))
            
#            return rotation_matrix(-(e+s),'z',False) *\
#                   rotation_matrix(d,'y',False) *\
#                   rotation_matrix(e,'z',False)
        
    @staticmethod
    def _CIOLocator(epoch):
        """
        Returns the CIO locator s for the provided epoch. s is the difference in
        RA between the GCRS and CIP points for the ascending node of the CIP 
        equator.
        """
        #from ..obstools import jd2000,epoch_to_jd
        from ..constants import asecperrad
        
        from .ephems import _mean_anomaly_of_moon,_mean_anomaly_of_sun,\
                            _mean_long_of_moon_minus_ascnode,_long_earth,\
                            _mean_elongation_of_moon_from_sun,_long_venus,\
                            _mean_long_asc_node_moon,_long_prec
        
        #first need to find x and y for the CIP, as s+XY/2 is needed
        B = ICRSCoordinates.frameBiasJ2000
        P = _precession_matrix_J2000_Capitaine(epoch)
        N = _nutation_matrix(epoch)
        
        #N*P*B takes GCRS to true, so CIP is bottom row
        x,y,z = (N*P*B).A[2]
        
        #T = (epoch_to_jd(epoch) - jd2000)/36525
        T = (epoch-2000)/100
        
        fundargs = [] #fundamental arguments
        
        fundargs.append(_mean_anomaly_of_moon(T))
        fundargs.append(_mean_anomaly_of_sun(T))
        fundargs.append(_mean_long_of_moon_minus_ascnode(T))
        fundargs.append(_mean_elongation_of_moon_from_sun(T))
        fundargs.append(_mean_long_asc_node_moon(T))
        fundargs.append(_long_venus(T))
        fundargs.append(_long_earth(T))
        fundargs.append(_long_prec(T))
        
        fundargs = np.array(fundargs)
        
        polys,orders = _CIO_locator_data
        newpolys = polys.copy() #copy 0-values to add to
        
        for i,o in enumerate(orders):
            ns,sco,cco = o
            a = np.dot(ns,fundargs)
            newpolys[i] += np.sum(sco*np.sin(a) + cco*np.cos(a))
        
        return np.polyval(newpolys[::-1],T)/asecperrad - x*y/2.0
    
    @CoordinateSystem.registerTransform(GCRSCoordinates,'self',transtype='smatrix')
    def _fromGCRS(gcrsc):
        return CIRSCoordinates._CMatrix(gcrsc.epoch)
    @CoordinateSystem.registerTransform('self',GCRSCoordinates,transtype='smatrix')
    def _toGCRS(cirssys):
        return CIRSCoordinates._CMatrix(cirssys.epoch).T
            
class EquatorialCoordinatesEquinox(EquatorialCoordinatesBase):
    """
    Represents an object in *mean* geocentric apparent equatorial coordinates,
    using the pre-IAU2000 systems where the plane of the ecliptic is the
    fundamental plane and the origin is at the equinox of date (as set by
    :attr:`epoch`).
    
    Changes to the :attr:`epoch` will result in the coordinates being updated
    for precession, but not nutation, nor annual abberation. Neither are
    planned by the primary author of this package, as IAU 2000 recommends using
    only CIO-based systems, but if someone actually wants equinox-based
    nutation, feel free to implement it and pass it along.
    
    To convert from these coordinates to :class:`HorizontalCoordinates`
    appropriate for observed coordinates, site information is necessary. Hence,
    the transformations from equatorial to horizontal coordinates are performed
    by the :class:`~astropysics.obstools.Site` class in the
    :mod:`astropysics.obstools` module, and attempting to directly convert will
    raise a :exc:`TypeError`.
    """
    
    transweight = 1.1 #Equinox-based transforms are to be avoided in favor of CIRS
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch using the
        IAU 2000 precessions from Capitaine, N. et al. 2003 as written in the
        USNO Circular 179.
        """
        if self.epoch is not None and newepoch is not None:
            #convert from current to J2000
            B = _nutation_matrix(self.epoch) *\
                _precession_matrix_J2000_Capitaine(self.epoch).T 
                #transpose==inv; matrix is real unitary
                
            #convert to new epoch
            A = _nutation_matrix(newepoch) *\
                _precession_matrix_J2000_Capitaine(newepoch)
                
            self.matrixRotate(A*B)
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
    @CoordinateSystem.registerTransform(GCRSCoordinates,'self',transtype='smatrix')
    def _fromGCRS(gcrsc):
        B = ICRSCoordinates.frameBiasJ2000
        if gcrsc.epoch is None:
            return B
        else:
            P = _precession_matrix_J2000_Capitaine(gcrsc.epoch)
            N = _nutation_matrix(gcrsc.epoch)
            return N*P*B 
    @CoordinateSystem.registerTransform('self',GCRSCoordinates,transtype='smatrix')
    def _toGCRS(eqsys):
        return EquatorialCoordinatesEquinox._fromGCRS(eqsys).T
          
    @CoordinateSystem.registerTransform('self',CIRSCoordinates,transtype='smatrix')
    def _toCIRS(eqsys):
        if eqsys.epoch is None:
            return np.eye(3).view(np.matrix)
        else:
            from ..obstools import epoch_to_jd
            from .funcs import equation_of_the_origins
            from ..utils import rotation_matrix
            
            jd = epoch_to_jd(eqsys.epoch)
            eqo = equation_of_the_origins(jd)*15.  #hours>degrees
            return rotation_matrix(-eqo,'z',True)
    
    @CoordinateSystem.registerTransform(CIRSCoordinates,'self',transtype='smatrix')
    def _fromCIRS(cirssys):
        return EquatorialCoordinatesEquinox._toCIRS(cirssys).T
            
                
class ITRSCoordinates(EpochalLatLongCoordinates):
    """
    Coordinates based on the International Terrestrial Reference System. The
    particular implementation here assumes ITRS matches the WGS84 coordinates
    (used by GPS) - for astronomical purposes, this is a perfectly good
    assumption.
    
    Epoch transformations in this system only adjust for polar motion - to
    account for earth rotation, transform back to
    :class:`CIRSCoordinates` or :class:`EquatorialCoordinatesEquinox`,
    change the epoch, then transfrom back to :class:`ITRSCoordinates`.
    
    Because polar motion is not fully predictable, a number of methods are
    available for approximating it. To choose a method, set the
    :attr:`ITRSCoordinates.polarmotion` class attribute -- this will also affect
    all future transformations to :class:`ITRSCoordinates` from other coordinate
    systems. The following are valid values:
    
        * None
            Assumes the pole locations are fixed to the CIP at all times, aside
            from the tiny effect of s' (the TIO-CIO shift).
        * A 2-tuple of callables (xp,yp)
            They will be called as xp(epoch) and yp(epoch) and the result will
            be assumed to give the x and y coordinates of the poles in the CIP
            frame.
    
    .. note:: 
        The transformations from CIRS and Equinox systems to ITRS technically
        involve converting to TIRS (the Terrestrial Intermediate Reference
        System), distinguished from ITRS by no polar motion correction. While
        there is no class explicitly representing TIRS, ITRS with :attr:`epoch`
        set to None is equivalent to TIRS.
        
    """
    
    __slots__ = tuple('_dpc') 
    #_dpc is included for transformation to/from Equatorial-like systems
    _longrange_ = (-180,180)
    
    polarmotion = None
    """
    Technique of computing poles (see :class:`ITRSCoordinates` documentation)
    """
    
    @staticmethod
    def _TIOLocator(epoch):
        """
        s-prime, the offset between the0 of longitude for the CIO of CIRS and
        the TIO of TIRS (Terrestrial Intermediate Reference System) - TIRS and
        ITRS differ by polar motion. Return value in radians.
        
        This is really,really small, and isn't really necessary except for
        completeness
        """
        from ..constants import asecperrad
        T = (epoch-2000)/100
        return -47e-6*T/asecperrad
    
    @staticmethod
    def _WMatrix(epoch):
        from ..utils import rotation_matrix
        
        sp = ITRSCoordinates._TIOLocator(epoch)
        if ITRSCoordinates.polarmotion is None:
            xp = 0
            yp = 0
        else: #assume 2-sequence (xp,yp)
            xp,yp = ITRSCoordinates.polarmotion
            if callable(xp):
                xp = xp(epoch)
            if callable(yp):
                yp = yp(epoch)
        
        
        #TODO: test if the following, linear matrix is good enough to not 
        #bother with the "right" one:
        #[[1,-sp,-xp], 
        # [sp,1,yp],
        # [xp,-yp,1]] #can also do sp->0
        return rotation_matrix(-yp,'x') *\
               rotation_matrix(-xp,'y') *\
               rotation_matrix(sp,'z') 
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`ITRSCoordinates` to a new epoch, adjusting the 
        coordinate values for polar motion.
        """
        if self.epoch is not None and newepoch is not None:
            M = self._WMatrix(self.epoch).T
            Mn = self._WMatrix(newepoch)
            
            self.matrixRotate(Mn*M)
            
        #this sets the epoch
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
    
    @add_docs(EpochalLatLongCoordinates.convert)
    def convert(self,tosys,optimize=_convertoptimizedefault):
        ''
        res = EpochalLatLongCoordinates.convert(self,tosys,optimize)
        if self._dpc is not None:
            res._dpc = self._dpc
        return res
    
    @CoordinateSystem.registerTransform(CIRSCoordinates,'self',transtype='smatrix')
    def _fromEqC(eqc):
        from .funcs import earth_rotation_angle
        from ..obstools import epoch_to_jd
        from ..utils import rotation_matrix
        
        epoch = eqc.epoch
        if epoch is not None:
            jd = epoch_to_jd(eqc.epoch)
            
            era = earth_rotation_angle(jd,degrees=True)
            W = ITRSCoordinates._WMatrix(eqc.epoch)
            
            return W*rotation_matrix(era)
        else:
            return np.eye(3).view(np.matrix)
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self',transtype='smatrix')
    def _fromEqE(eqe):
        from .funcs import greenwich_sidereal_time
        from ..utils import rotation_matrix
        from ..obstools import epoch_to_jd
        
        epoch = eqe.epoch
        if epoch is not None:
            jd = epoch_to_jd(eqe.epoch)
            try:
                gst = greenwich_sidereal_time(jd,True)*15. #hours -> degrees
            except Exception,e:
                from warnings import warn
                warn('temporarily bypassing problem with greenwich_sidereal_time:%s'%e)
                gst = greenwich_sidereal_time(jd,'simple')*15. #hours -> degrees
            W = ITRSCoordinates._WMatrix(eqe.epoch)
            
            return W*rotation_matrix(gst)
        else:
            return np.eye(3).view(np.matrix)  
    
    @CoordinateSystem.registerTransform('self',CIRSCoordinates,transtype='smatrix')
    def _toEqC(itrsc):
        #really we want inverse, but rotations are unitary -> inv==transpose
        #we provide itrsc in the call because the epoch is needed
        return ITRSCoordinates._fromEqC(itrsc).T 
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox,transtype='smatrix')
    def _toEqE(itrsc):
        #really we want inverse, but rotations are unitary -> inv==transpose
        #we provide itrsc in the call because the epoch is needed
        return ITRSCoordinates._fromEqE(itrsc).T 
            
class FK5Coordinates(EquatorialCoordinatesEquinox):
    """
    Equatorial Coordinates fixed to the FK5 reference system. 
    """
    
    __slots__ = tuple()
    
    @staticmethod
    def _precessionMatrixJ(epoch1,epoch2):
        """
        Computes the precession matrix from one Julian epoch to another
        """
        from ..utils import rotation_matrix
        
        T = (epoch1 - 2000)/100
        dt = (epoch2 - epoch1)/100
        
        pzeta = (0.017998,0.000344,0.30188,-0.000139,1.39656,2306.2181)
        temp = pzeta[5] + T*(pzeta[4]+T*pzeta[3])
        zeta = dt*(temp + dt*((pzeta[2]+pzeta[1]*T) + dt*pzeta[0]))/3600
        
        pz = (0.018203,-0.000066,1.09468)
        z = dt*(temp + dt*((pz[2]+pz[1]*T) + dt*pz[0]))/3600
        
        ptheta = (-0.041833,-0.000217,-0.42665,-0.000217,-0.85330,2004.3109)
        temp = ptheta[5] + T*(ptheta[4]+T*ptheta[3])
        theta = dt*(temp + dt*((ptheta[2]+ptheta[1]*T) + dt*ptheta[0]))/3600
        
        return rotation_matrix(-z,'z') *\
               rotation_matrix(theta,'y') *\
               rotation_matrix(-zeta,'z')
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch. Uses the
        algorithm resolved by the IAU in 1976 as written in Meeus, as well as
        Lieske, J.H. 1979.
        
        According to SOFA, this for becomes less valid at the following levels:
        
        * 1960 CE  to 2040 CE
            < 0.1"
        * 1640 CE to 2360 CE
            < 1"
        * 500 BCE to 3000 CE
            < 3"
        * 1200 BCE to 3900 CE
            > 10"
        * 4200 BCE to 5600 CE
            > 100"
        * 6800 BCE to 8200 CE
            > 1000"
            
        """
        if self.epoch is not None and newepoch is not None:
            self.matrixRotate(self._precessionMatrixJ(self.epoch,newepoch))
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
    @CoordinateSystem.registerTransform(ICRSCoordinates,'self',transtype='smatrix')
    def _fromICRS(icrsc):
        """
        B-matrix from USNO circular 179 
        """
        from ..utils import rotation_matrix
        
        eta0 = -19.9/3600000
        xi0 = 9.1/3600000
        da0 = -22.9/3600000
        B = rotation_matrix(-eta0,'x') *\
            rotation_matrix(xi0,'y') *\
            rotation_matrix(da0,'z')
            
        epoch = icrsc.epoch
        if icrsc.epoch is None:
            return B
        else:
            return FK5Coordinates._precessionMatrixJ(2000,icrsc.epoch)*B
    
    @CoordinateSystem.registerTransform('self',ICRSCoordinates,transtype='smatrix')
    def _toICRS(fk5c):
        return FK5Coordinates._fromICRS(fk5c).T
    
class FK4Coordinates(EquatorialCoordinatesEquinox):
    """
    Equatorial Coordinates fixed to the FK4 reference system.  Note that this 
    implementation does *not* correct for the elliptic terms of aberration
    as of yet.
    
    Epoch is Besselian.
    """
    
    __slots__ = tuple()
    
    julianepoch = False
    
    def __init__(self,*args,**kwargs):
        """
        Input for FK4 coordinates. Can follow any of the following forms:
        
        * EquatorialCoordinatesBase()
        * EquatorialCoordinatesBase(:class:`EquatorialCoordinatesBase`)
        * EquatorialCoordinatesBase('rastr decstr')
        * EquatorialCoordinatesBase((ra,dec))
        * EquatorialCoordinatesBase(ra,dec)
        * EquatorialCoordinatesBase(ra,fdec,raerr,decerr)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch,distancepc)
        
        The epoch of FK4 coordinates defaults to B1950.
        """
        args = list(args)
        args.insert(0,self)
        EquatorialCoordinatesEquinox.__init__(*args,**kwargs)
        if self._epoch==2000.:
            self._epoch = 1950.
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch. Uses the
        method of Newcomb (pre-IAU1976) to compute precession.
        """
        if self.epoch is not None and newepoch is not None:
            self.matrixRotate(self._precessionMatrixB(self.epoch,newepoch))
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    @staticmethod
    def _precessionMatrixB(epoch1,epoch2):
        """
        computes the precession matrix from one Besselian epoch to another using
        Newcomb's method.
        """
        from ..utils import rotation_matrix
        
        #tropical years
        t1 = (epoch1-1850.0)/1000.0    
        t2 = (epoch2-1850.0)/1000.0
        dt = t2 - t1
        
        zeta1 = 23035.545 + t1*139.720+0.060*t1*t1
        zeta2 = 30.240 - 0.27*t1
        zeta3 = 17.995
        pzeta = (zeta3,zeta2,zeta1,0)
        zeta = np.polyval(pzeta,dt)/3600
        
        z1 = 23035.545 + t1*139.720 + 0.060*t1*t1
        z2 = 109.480 + 0.39*t1
        z3 = 18.325
        pz = (z3,z2,z1,0)
        z = np.polyval(pz,dt)/3600
        
        theta1 = 20051.12 - 85.29*t1 - 0.37*t1*t1
        theta2 = -42.65 - 0.37*t1
        theta3 = -41.8
        ptheta = (theta3,theta2,theta1,0)
        theta = np.polyval(ptheta,dt)/3600
        
        
        return rotation_matrix(-z,'z') *\
               rotation_matrix(theta,'y') *\
               rotation_matrix(-zeta,'z')
        
               
    @CoordinateSystem.registerTransform('self',FK5Coordinates,transtype='smatrix')
    def _toFK5(fk4c):
        from ..obstools import epoch_to_jd,jd_to_epoch
        
        
        #B1950->J2000 matrix from Murray 1989 A&A 218,325
        B = np.mat([[0.9999256794956877,-0.0111814832204662,-0.0048590038153592],
                    [0.0111814832391717,0.9999374848933135,-0.0000271625947142],
                    [0.0048590037723143,-0.0000271702937440,0.9999881946023742]])
        
        if fk4c.epoch is not None and fk4c.epoch != 1950:
            jd = epoch_to_jd(fk4c.epoch,False)
            jepoch = jd_to_epoch(jd)
            T = (jepoch - 1950)/100
            
            #now add in correction terms for FK4 rotating system
            B[0,0] += -2.6455262e-9*T
            B[0,1] += -1.1539918689e-6*T
            B[0,2] += 2.1111346190e-6*T
            B[1,0] += 1.1540628161e-6*T
            B[1,1] += -1.29042997e-8*T
            B[1,2] += 2.36021478e-8*T
            B[2,0] += -2.1112979048e-6*T
            B[2,1] += -5.6024448e-9*T
            B[2,2] += 1.02587734e-8*T
            
            PB = FK4Coordinates._precessionMatrixB(fk4c.epoch,1950)
            
            return B*PB
        else:
            return B
    
    @CoordinateSystem.registerTransform(FK5Coordinates,'self',transtype='smatrix')
    def _fromFK5(fk5c):
        #need inverse because Murray's matrix is *not* a true rotation matrix
        return FK4Coordinates._toFK5(fk5c).I
        
class EclipticCoordinatesCIRS(EpochalLatLongCoordinates):
    """
    Ecliptic Coordinates (beta, lambda) such that the fundamental plane passes
    through the ecliptic at the current epoch.
    
    Note that because the concept of the ecliptic can be complicated or even
    ill-defined, ecliptic coordinates is astropysics are simply defined as
    tied to a particular set of equatorial coordinates with a given obliquity
    model.  For :class:`EclipticCoordinatesCIRS`, the equatorial 
    coordinates are :class:`CIRSCoordinates` with obliquity
    given by the IAU 2006 obliquity model (see :func:`obliquity`)
    """
    
    __slots__ = ()
    _longlatnames_ = ('lamb','beta')
    _longrange_ = (0,360)
    
    obliqyear = 2006
    
    def __init__(self,lamb=0,beta=0,lamberr=None,betaerr=None,epoch=2000,
                      distanceau=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        EpochalLatLongCoordinates.__init__(self,lamb,beta,lamberr,betaerr,epoch)
        self.distanceau = distanceau
        
    @CoordinateSystem.registerTransform('self',CIRSCoordinates,transtype='smatrix')
    def _toEq(eclsc):
        from .funcs import obliquity
        from ..utils import rotation_matrix
        
        return rotation_matrix(-obliquity(eclsc.jdepoch,EclipticCoordinatesCIRS.obliqyear),'x')
        
    @CoordinateSystem.registerTransform(CIRSCoordinates,'self',transtype='smatrix')
    def _fromEq(eqc):
        from .funcs import obliquity
        from ..utils import rotation_matrix
        
        return rotation_matrix(obliquity(eqc.jdepoch,EclipticCoordinatesCIRS.obliqyear),'x')
    
    def transformToEpoch(self,newepoch):
        if self.epoch is not None and newepoch is not None:
            eqc = self.convert(CIRSCoordinates)
            eqc.epoch = newepoch
            newval = eqc.convert(self.__class__)
            self._lat._decval = newval._lat._decval
            self._long._decval = newval._long._decval
            self._laterr._decval = newval._lat._decval
            self._longerr._decval = newval._longerr._decval
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
        
class EclipticCoordinatesEquinox(EpochalLatLongCoordinates):
    """
    Ecliptic Coordinates (beta, lambda) such that the fundamental plane passes
    through the ecliptic at the current epoch.
    
    Note that because the concept of the ecliptic can be complicated or even
    ill-defined, ecliptic coordinates is astropysics are simply defined as tied
    to a particular set of equatorial coordinates with a given obliquity model.
    For :class:`EclipticCoordinatesEquinox`, the equatorial coordinates are
    :class:`EquatorialCoordinatesEquinox` with obliquity given by the IAU 1980
    obliquity model (see :func:`~astropysics.coords.funcs.obliquity`)
    """
    
    __slots__ = ()
    _longlatnames_ = ('lamb','beta')
    _longrange_ = (0,360)
    
    obliqyear = 1980
    
    def __init__(self,lamb=0,beta=0,lamberr=None,betaerr=None,epoch=2000,
                      distanceau=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        EpochalLatLongCoordinates.__init__(self,lamb,beta,lamberr,betaerr,epoch)
        self.distanceau = distanceau
        
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox,transtype='smatrix')
    def _toEq(eclsc):
        from .funcs import obliquity
        from ..utils import rotation_matrix
        
        return rotation_matrix(-obliquity(eclsc.jdepoch,EclipticCoordinatesEquinox.obliqyear),'x')
        
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self',transtype='smatrix')
    def _fromEq(eqc):
        from .funcs import obliquity
        from ..utils import rotation_matrix
        
        return rotation_matrix(obliquity(eqc.jdepoch,EclipticCoordinatesEquinox.obliqyear),'x')
        
    def transformToEpoch(self,newepoch):
        if self.epoch is not None and newepoch is not None:
            eqc = self.convert(EquatorialCoordinatesEquinox)
            eqc.epoch = newepoch
            newval = eqc.convert(self.__class__)
            self._lat._decval = newval._lat._decval
            self._long._decval = newval._long._decval
            self._laterr._decval = newval._lat._decval
            self._longerr._decval = newval._longerr._decval
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
class RectangularGeocentricEclipticCoordinates(RectangularCoordinates,EpochalCoordinates):
    """
    Rectangular coordinates oriented so that the x-y plane lies in the plane of
    the ecliptic at the specified epoch. Distances are in AU. Origin is at the
    center of mass of the Earth.
    
    Note that the epoch should not be set directly - if precession is desired
    desired, convert to an Ecliptic coordinate system, do the precession, and
    convert back.
    """
    __slots__ = tuple()
    julianepoch = True
    
    def __init__(self,x,y,z,epoch=None):
        RectangularCoordinates.__init__(self,x,y,z)
        self._epoch = epoch
        
    def __getstate__(self):
        d = RectangularCoordinates.__getstate__(self)
        d.update(EpochalCoordinates.__getstate__(self))
        return d
    
    def __setstate__(self,d):
        RectangularCoordinates.__setstate__(self,d)
        EpochalCoordinates.__setstate__(self,d)
        
    def __str__(self):
        if self.epoch is None:
            epochstr = ''
        else:
            epochstr = ' ('+self.epochstr+')'
        return RectangularCoordinates.__str__(self) + epochstr
        
    def transformToEpoch(self,newepoch):
        EpochalCoordinates.transformToEpoch(newepoch)
    
    @CoordinateSystem.registerTransform('self',EclipticCoordinatesCIRS)
    def _toEcC(rec):
        from math import asin,atan2,degrees
        
        x,y,z = rec.x,rec.y,rec.z        
        r = (x*x+y*y+z*z)**0.5
        beta = degrees(asin(z/r))
        lamb = degrees(atan2(y,x))
        
        return EclipticCoordinatesCIRS(lamb,beta,distanceau=r,epoch=rec.epoch)
    
    @CoordinateSystem.registerTransform('self',EclipticCoordinatesEquinox)
    def _toEcQ(rec):
        from math import asin,atan2,degrees
        
        x,y,z = rec.x,rec.y,rec.z        
        r = (x*x+y*y+z*z)**0.5
        beta = degrees(asin(z/r))
        lamb = degrees(atan2(y,x))
        
        return EclipticCoordinatesEquinox(lamb,beta,distanceau=r,epoch=rec.epoch)
    
    @CoordinateSystem.registerTransform(EclipticCoordinatesCIRS,'self')    
    def _fromEcC(ec):
        from math import sin,cos,degrees
        
        l,b = ec.lamb.r,ec.beta.r
        if ec.distanceau is None:
            r = 1
        else:
            r = ec.distanceau[0]
            
        x = r*cos(l)*cos(b)
        y = r*sin(l)*cos(b)
        z = r*sin(b)
        
        return RectangularGeocentricEclipticCoordinates(x,y,z,ec.epoch)
    
    @CoordinateSystem.registerTransform(EclipticCoordinatesEquinox,'self')    
    def _fromEcQ(ec):
        from math import sin,cos,degrees
        
        l,b = ec.lamb.r,ec.beta.r
        if ec.distanceau is None:
            r = 1
        else:
            r = ec.distanceau[0]
            
        x = r*cos(l)*cos(b)
        y = r*sin(l)*cos(b)
        z = r*sin(b)
            
        return RectangularGeocentricEclipticCoordinates(x,y,z,ec.epoch)
    
class GalacticCoordinates(EpochalLatLongCoordinates):
    __slots__ = tuple()
    _longlatnames_ = ('l','b')
    _longrange_ = (0,360)
    
    _ngp_J2000 = FK5Coordinates(192.859508, 27.128336,epoch=2000)
    _long0_J2000 = AngularCoordinate(122.932)
    _ngp_B1950 = FK4Coordinates(192.25, 27.4,epoch=1950)
    _long0_B1950 = AngularCoordinate(123)
    
    def transformToEpoch(self,newepoch):
        """
        Galactic coordinates are nominally inertial, although the definition is
        a bit unclear in that regard.
        """
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    @CoordinateSystem.registerTransform(FK5Coordinates,'self',transtype='smatrix')
    def _fromFK5(fk5coords):
        from ..utils import rotation_matrix
        
        epoch = 2000 if fk5coords.epoch is None else fk5coords.epoch
        
        mat = rotation_matrix(180 - GalacticCoordinates._long0_J2000.d,'z') *\
              rotation_matrix(90 - GalacticCoordinates._ngp_J2000.dec.d,'y') *\
              rotation_matrix(GalacticCoordinates._ngp_J2000.ra.d,'z') *\
              FK5Coordinates._precessionMatrixJ(epoch,2000)
        mat.nocache = True #can't cache because of the need to get the epoch
        return mat
    
    @CoordinateSystem.registerTransform('self',FK5Coordinates,transtype='smatrix')
    def _toFK5(galcoords):
        return GalacticCoordinates._fromFK5(galcoords).T
    
    @CoordinateSystem.registerTransform(FK4Coordinates,'self',transtype='smatrix')
    def _fromFK4(fk4coords):
        from ..utils import rotation_matrix
        
        epoch = 1950 if fk4coords.epoch is None else fk4coords.epoch
        
        mat = rotation_matrix(180 - GalacticCoordinates._long0_B1950.d,'z') *\
              rotation_matrix(90 - GalacticCoordinates._ngp_B1950.dec.d,'y') *\
              rotation_matrix(GalacticCoordinates._ngp_B1950.ra.d,'z') *\
              FK4Coordinates._precessionMatrixB(epoch,1950)
        mat.nocache = True #can't cache because of the need to get the epoch
        return mat
    
    @CoordinateSystem.registerTransform('self',FK4Coordinates,transtype='smatrix')
    def _toFK4(galcoords):
        return GalacticCoordinates._fromFK4(galcoords).T
        
class SupergalacticCoordinates(EpochalLatLongCoordinates):   
    __slots__ = tuple()
    _longlatnames_ = ('sgl','sgb')
    _longrange_ = (0,360)
    
    _nsgp_gal = GalacticCoordinates(47.37,6.32) #glactc is 47.47=WRONG
    _sg0_gal = GalacticCoordinates(137.37,0)
    
    def transformToEpoch(self,newepoch):
        """
        Supergalactic coordinates are nominally inertial, although the
        definition is a bit unclear in that regard.
        """
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    @CoordinateSystem.registerTransform('self',GalacticCoordinates,transtype='smatrix')
    def _toGal(sgalcoords):
        return SupergalacticCoordinates._fromGal(sgalcoords).T
    
    @CoordinateSystem.registerTransform(GalacticCoordinates,'self',transtype='smatrix')
    def _fromGal(galcoords):
        from ..utils import rotation_matrix
        
        z1r = rotation_matrix(SupergalacticCoordinates._nsgp_gal.l.d,'z')
        yr = rotation_matrix(90 - SupergalacticCoordinates._nsgp_gal.b.d,'y')
        z2r = rotation_matrix(180 - SupergalacticCoordinates._sg0_gal.l.d +\
                              SupergalacticCoordinates._nsgp_gal.l.d,'z')
        return z2r*yr*z1r
    
class HorizontalCoordinates(LatLongCoordinates):
    """
    This object represents an angular location on the unit sphere, with the 
    north pole of the coordinate position fixed to the local zenith
    
    To convert from other :class:`Coordinate` types to horizontal positions, see 
    :class:`astropysics.obstools.Site`, as site information is required for
    these corrections
    """  
    __slots__ = tuple()
    _longlatnames_ = ('az','alt')
    _longrange_ = (0,360)
    
    def __init__(self,alt=0,az=0,alterr=None,azerr=None,distancepc=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        LatLongCoordinates.__init__(self,az,alt,azerr,alterr,distancepc)
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self')
    @CoordinateSystem.registerTransform(CIRSCoordinates,'self')
    def _toHoriz(incoosys=None):
        raise TypeError('use astropysics.obstools.Site methods to transform celestial to terrestrial coordinates')
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox)
    @CoordinateSystem.registerTransform('self',CIRSCoordinates)
    def _fromHoriz(incoosys=None):
        raise TypeError('use astropysics.obstools.Site methods to transform terrestrial to celestial coordinates')
    
    
#Now that all the coordinate systems have been made, add the diagram to the docs
#That shows the graph of the built-in transforms

postbuiltin = """
A similar diagram can be generated after the user has created and registered
custom coordinates and transforms::
    
    from networkx import to_agraph,relabel_nodes,draw_networkx
    from astropysics.coords import CoordinateSystem
    graph = CoordinateSystem.getTransformGraph()
    dotgraph = to_agraph(relabel_nodes(graph,lambda n:n.__name__))
    dotgraph.graph_attr.update(dict(size='12.0, 12.0',fontsize=12))
    dotgraph.write('mygraph.dot')
    draw_networkx(graph)
    
This will save a graphviz dot file and displaying the graph with matplotlib,
showing both builtin and custom-added coordinates and transforms.
    
"""

try:
    from networkx import to_agraph,relabel_nodes
    graph = to_agraph(relabel_nodes(CoordinateSystem.getTransformGraph(),lambda n:n.__name__))
    graph.graph_attr.update(dict(size=r'12.0, 12.0',fontsize=12))
    transstr="""
Built-in Transforms
^^^^^^^^^^^^^^^^^^^

A number of coordinate systems are provided built into astropysics. Most of
these have pre-defined standard transformations. The built-in coordinate classes
with defined transformation are shown in the diagram below.

.. graphviz::

    """+graph.string().replace('\n','\n    ')+postbuiltin
    __doc__ = __doc__.replace('{transformdiagram}',transstr)
    del to_agraph,relabel_nodes,graph
except ImportError:
    #if networkx or pygraphviz isn't present, drop the diagram but add a warning that it's missing
    warningstr = """
Builtin Coordinate System Transforms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
    A diagram showing the relationships between the pre-defined transformations
    should be here, but this copy of the documentation was built without
    `networkx <http://networkx.lanl.gov/>` and `pygraphviz
    <http://networkx.lanl.gov/pygraphviz/>` available to build the diagram.
    Please re-build this file after those packages are installed to see the
    diagram.
    """+postbuiltin
    __doc__ = __doc__.replace('{transformdiagram}',warningstr)
    del warningstr
    
    
#<--------------------------Convinience Functions------------------------------>

    
def angular_string_to_dec(instr,hms=True,degrees=True):
    """
    Convinience function to convert a angular coordinate string to a decimal
    value.
    
    :param hms: 
        If True, the coordinate will be assumed to be h:m:s, otherwise d:m:s.
        This will be ignored if the coordinates are specified as ##h##m##s or
        ##d##m##s, or if the input is not in sexigesimal form.
    :type hms: boolean
    :param degrees: 
        If True, the output will be decimal degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: Decimal value in radians or degrees
    """
    ac = AngularCoordinate(instr)
    if degrees:
        return ac.degrees
    else:
        return ac.radians


def objects_to_coordinate_arrays(posobjs,coords='auto',degrees=True):
    """
    converts a sequence of position objects into an array of coordinates.  
    
    `coords` determines the order of the output coordinates - it can be a 
    comma-seperated list of coordinate names or if 'auto', it will be 'lat,long'
    for all coordinate systems except for Equatorial, which will use 'ra,dec'
    
    if `degrees` is True, returned arrays are in degrees, otherwise radians
    """
    if coords=='auto':
        coordnames = None
    else:
        coordnames = coords.split(',')
        
    coords = []
    if degrees:
        for o in posobjs:
            if coordnames is None:
                if isinstance(o,EquatorialCoordinates):
                    coords.append((o.ra.d,o.dec.d))
                else:
                    coords.append((o.lat.d,o.long.d))
            else:
                coords.append([getattr(o,c).d for c in coordnames])
    else:
        for o in posobjs:
            if coordnames is None:
                if isinstance(o,EquatorialCoordinates):
                    coords.append((o.ra.r,o.dec.r))
                else:
                    coords.append((o.lat.r,o.long.r))
            else:
                coords.append([getattr(o,c).r for c in coordnames])

    return np.array(coords).T

