#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""

======
coords
======

The :mod:`coords` module contains classes and functions specifying locations of
objects on the sky as well as coordinate converstions and distance computations,
including cosmological distance and redshift calculations.

Some of the calculations involved make use of the currently selected cosmology
(see :mod:`astropysics.constants`) and hence may not function as expected if a
particularly strange cosmology is in use.

.. seealso::

    `Kapteyn libraries <http://www.astro.rug.nl/software/kapteyn/index.html>`_
        A set of python libraries with excellent coordinate transform and other
        capabilities.
        
    `Meeus, Jean H. "Astronomical Algorithms" ISBN 0943396352 <http://www.willbell.com/MATH/mc1.htm>`_ 
        An authoritative reference on coordinates, ephemerides, and related
        transforms in astronomy.
   
    `Standards Of Fundamental Astronomy (SOFA) <http://www.iausofa.org/>`_ 
        The IAU reference implementations for coordinates and earth rotation.
   
    `USNO Circular 179 <http://aa.usno.navy.mil/publications/docs/Circular_179.pdf>`_ 
        An excellent description of the IAU 2000 resolutions and related
        background for defining ICRS, CIO, and related standards.
   
.. note::
    Timekeeping/conversion functions are kept in :mod:`astropysics.obstools`,
    although some are crucial to coordinate systems.
   
.. warning:: 
    While the framework for the coordinate transformations is done, the
    actual implementation is still a work in progress, so be aware that some of 
    the transformations will be throwing :exc:`NotImplementedError` until 
    everything is in place.

.. todo:: Tutorials


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.coords
   :parts: 1

Module API
----------

"""

#TODO: WCSlib or similar support a la Kapteyn?

from __future__ import division,with_statement

from .constants import pi,asecperrad
from .utils import rotation_matrix
from .io import _get_package_data
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

#<----------------coordinate classes and related functions------------------>




class AngularCoordinate(object):
    """
    The value of an angular coordinate on the unit sphere.
    
    Arithmetic operators can be applied to the coordinate, and will be applied 
    directly to the numerical value in radians.  For + and -, two angular 
    coordinates may be used, although for -, an AngularSeperation object will
    be returned.
    
    """
    import re as _re
    __slots__=('_decval','_range')
    
    __acregex = _re.compile(r'(?:([+-])?(\d+(?:[.]\d*)?)(hours|h|degrees|d|radians|rads|rad|r| |:(?=\d+:\d+[.]?\d*$)))?(?:(\d+(?:[.]\d*)?)(m|\'|[:]| ))?(?:(\d+(?:[.]\d*)?)(s|"|$))?$')
    __decregex = _re.compile(r'[+-]\d+([.]\d*)') 
    
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
        
        .. testsetup::
        
            from astropysics.coords import AngularCoordinate
            from math import pi
        
        .. doctest::
        
            >>> ac = AngularCoordinate(2.5)
            >>> print ac
            +2d30'00.00"
            >>> print AngularCoordinate(ac)
            +2d30'00.00"
            >>> print AngularCoordinate(pi,radians=True)
            +180d00.00"
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
                    raise ValueError('invalid or ambiguous string input for AngularCoordinate')
            else:
                decm = self.__decregex.match(sinpt)
                if decm:
                    if radians:
                        self.r = float(decm.group(0))
                    else:
                        self.d = float(decm.group(0))
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
        self.degrees=abs(dms[0])+abs(dms[1])/60.+abs(dms[2])/3600.
        if dms[0]<0:
            self._decval*=-1
    def _getDegminsec(self):
        fulldeg=abs(self.degrees)
        deg=int(fulldeg)
        fracpart=fulldeg-deg
        min=int(fracpart*60.)
        sec=fracpart*3600.-min*60.
        return -deg if self.degrees < 0 else deg,min,sec
    degminsec = property(_getDegminsec,_setDegminsec,doc="""
    The value of this :class:`AngularCoordinate` as an (degrees,minutes,seconds)
    tuple, with degrees and minutes as integers and seconds as a float.    
    """)
    dms = degminsec
    
    def _setHrsminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set hrsminsec as a length-3 iterator')
        self.degrees=15*(dms[0]+dms[1]/60.+dms[2]/3600.)
    def _getHrsminsec(self):
        factorized=self.degrees/15.
        hrs=int(factorized)
        mspart=factorized-hrs
        min=int(mspart*60.)
        sec=mspart*3600.-min*60.
        return hrs,min,sec
    hrsminsec = property(_getHrsminsec,_setHrsminsec,doc="""
    The value of this :class:`AngularCoordinate` as an (hours,minutes,seconds)
    tuple, with hours and minutes as integers and seconds as a float.    
    """)
    hms = hrsminsec
    
    def _setDecdeg(self,deg):
        rads = deg*pi/180.
        if self.range is not None:
            self._checkRange(rads)
        self._decval = rads
    def _getDecdeg(self):
        return self._decval*180/pi
    degrees = property(_getDecdeg,_setDecdeg,doc="""
    The value of this :class:`AngularCoordinate` in decimal degrees.   
    """)
    d = degrees    
    
    def _setRad(self,rads):
        if self.range is not None:
            self._checkRange(rads)
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
            self._checkRange(rads)
        self._decval = rads
    def _getDechr(self):
        return self._decval*12/pi
    hours = property(_getDechr,_setDechr,doc="""
    The value of this :class:`AngularCoordinate` in decimal hours.   
    """)
    h = hours
    
    def _checkRange(self,rads):
        if self._range is not None:
            low,up = self._range
            if not low <= rads <= up:
                raise ValueError('Attempted to set angular coordinate outside range')
    def _setRange(self,newrng):
        oldrange = self._range        
        try:
            if newrng is None:
                self._range = None
            else:
                from math import radians
                newrng = tuple(newrng)
                if len(newrng) != 2:
                    raise ValueError('range is not a 2-sequence')
                elif newrng[0] > newrng[1]:
                    raise ValueError('lower edge of range is not <= upper')
                newrng = (radians(newrng[0]),radians(newrng[1]))
            self._range = newrng
            self._checkRange(self._decval)
        except ValueError:
            self._range = oldrange
            raise ValueError('Attempted to set range when value is out of range')
    def _getRange(self):
        if self._range is None:
            return None
        else:
            from math import degrees
            return degrees(self._range[0]),degrees(self._range[1])
    range=property(_getRange,_setRange,doc="""
    The acceptable range of angles for this :class:`AngularCoordinate` as a
    2-tuple (lower,upper) with both angles in degrees. Will raise a
    :exc:`ValueError` if the current value is outside the range.
    """)
    
    ############################################################################
    #                                                                          #
    #                                                                          #
    #  and then a monkey showed up                                             #
    #  and gave you a kiss                                                     #
    #  and you were confused                                                   #
    #                                                                          #
    #                                             comment poem by Rie O:-)     #
    ############################################################################
   
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
        if type(other) == AngularCoordinate:
            from math import degrees
            res = AngularSeperation()
            return AngularSeperation(degrees(other._decval),degrees(self._decval))
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
        
    def getDmsStr(self,secform='%05.2f',sep=(unichr(176),"'",'"'), sign=True, canonical=False):
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
        
        :returns: String representation of this object.
        """
        d,m,s = self.degminsec
        
        if canonical:
            sgn = '' if self._decval < 0 else '+'
            return '%s%02.i:%02.i:%05.2f'%(sgn,d,m,s)
        
        d,m=str(d),str(m)
        
        s = secform%s
        
        if isinstance(sep,basestring):
            if sep == 'dms':
                sep = ('d','m','s')
            sep = (sep,sep)
        
        tojoin = []
        
        if sign and self._decval  >= 0:
            tojoin.append('+')
        
        if d is not '0':
            tojoin.append(d)
            tojoin.append(sep[0])
                
        if m is not '0':
            tojoin.append(m)
            tojoin.append(sep[1])
                
        tojoin.append(s)
        if len(sep)>2:
            tojoin.append(sep[2])
            
        return ''.join(tojoin)
        
    def getHmsStr(self,secform = None,sep = ('h','m','s'), canonical = False):
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
        
        :returns: String representation of this object.
        """
        
        h,m,s = self.hrsminsec
        
        if canonical:
            return '%02.i:%02.i:%06.3f'%(h,m,s)
        
        h,m=str(h),str(m)
        if secform is None:
            s = str(s)
        else:
            s = secform%s
        
        if isinstance(sep,basestring):
            if sep == 'hms':
                sep = ('h','m','s')
            sep = (sep,sep)
        
        tojoin = []
        
        tojoin.append(h)
        tojoin.append(sep[0])
            
        tojoin.append(m)
        tojoin.append(sep[1])
                
        tojoin.append(s)
        if len(sep)>2:
            tojoin.append(sep[2])
            
        return ''.join(tojoin)
    
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

class AngularSeperation(AngularCoordinate):
    """
    This class represents a seperation between two angular coordinates on the
    unit sphere.
    
    A constructor is available, but the most natural way to generate this object
    is to use the subtraction (-) operator on two :class:`AngularCoordinate`
    objects or two :class:`LatLongCoordinates` objects.
    """
    
    def __init__(self,*args):
        """
        Input arguments can be either:
        
        * AngularSeperation(:class:`AngularSeperation` object) 
            Generates a copy of the provided object.
        * AngularSeperation(sep) 
            Generates a seperation of the provided distance with no starting point.
        * AngularSeperation(start,end) 
            Computes the seperation from the start and end objects, which must
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
            raise ValueError('improper number of inputs to AngularSeperation')
        
        super(AngularSeperation,self).__init__(sep)
        
    def __add__(self,other):
        if isinstance(other,AngularCoordinate) and not self.__class__ == other.__class__:
            res = AngularCoordinate()
            res._decval = self._decval+other._decval
            return res
        else:
            return super(AngularSeperation,self).__add__(other)
        
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
    
        
    def projectedSeperation(self,zord,usez=False,**kwargs):
        """
        Computes the physical projected seperation assuming a given distance.
        
        kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
        
        :param zord: Redshift or distance
        :type zord: scalar number
        :param usez:
            If True, the input will be interpreted as a redshift, and kwargs
            will be passed into the distance calculation. The result will be in
            pc. Otherwise, `zord` will be interpreted as a distance.
        :type usez: boolean
        
        :returns: a float value for the seperation (in pc if redshift is used) 
        """
        return angular_to_physical_size(self.arcsec,zord,usez=usez,**kwargs)
    
    def seperation3d(self,zord1,zord2,usez=False,**kwargs):
        """
        computes the 3d seperation assuming the two points at the ends of this
        :class:`AngularSeperation` are at the distances `zord1` and `zord2`.  
        
        :param zord1: Redshift or distance for start point
        :type zord1: scalar number
        :param zord2: Redshift or distance for end point
        :type zord2: scalar number
        :param usez:
            If True, the inputs will be interpreted as a redshift, and kwargs
            will be passed into the distance calculation. The result will be in
            pc. Otherwise, `zord` will be interpreted as a distance.
        :type usez: boolean
        
        :returns: a float value for the seperation (in pc if redshift is used) 
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
                        
class CoordinateSystem(object):
    """
    Base class of all coordinate systems.  Contains machinery for managing 
    conversion between coordinate systems. Subclasses must override
    :meth:`__init__`
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
                raise ValueError('function already exists to convert {0} to {1}'.format(fromclass,toclass))
            
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
                func.transtype = transtype
                CoordinateSystem._converters[fromclass][toclass] = func
        
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
        
        **Examples**
        
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
        
    def convert(self,tosys):
        """
        converts the coordinate system from it's current system to a new 
        :class:`CoordinateSystem` object.
        
        :param tosys: The new coordinate system 
        :type tosys: A subclass of :class:`CoordinateSystem`
        
        :except: raises :exc:`NotImplementedError` if conversion is not present
        """
        try:
            return CoordinateSystem._converters[self.__class__][tosys](self)
        except KeyError:
            strf = 'cannot convert coordinate system {0} to {1}'
            raise NotImplementedError(strf.format(self.__class__.__name__,tosys))
    

class RectangularCoordinates(CoordinateSystem):
    """
    Rectangular/Cartesian Coordinates in three dimensions. Coordinates are
    accessed via the attributes :attr:`x`, :attr:`y`, and :attr:`z`.
    """
    
    __slots__ = ('x','y','z')
    
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        #: description of z
        self.z = z
        
    def __getstate__(self):
        #TODO: watch if this creates probelms by not being a dict
        return dict(x=self.x,y=self.y,z=self.z)
    
    def __setstate__(self,d):
        self.x = d['x']
        self.y = d['y']
        self.z = d['z']
    
CartesianCoordinates = RectangularCoordinates

class RectangularEclipticCoordinates(RectangularCoordinates):
    """
    Rectangular coordinates oriented so that the x-y plane lies in the plane of
    the ecliptic at the specified epoch. Distances are in AU. Can be either
    geocentric or heliocentric.
    
    Note that the epoch cannot be set directly - if precession is desired
    desired, convert to some other coordinate system, set the epoch, and convert
    back.
    """
    __slots__ = ('_geo','_epoch')
    
    def __init__(self,x,y,z,geocentric,epoch=None):
        RectangularCoordinates.__init__(x,y,z)
        self.geocentric = geocentric
        self._epoch = epoch
        
    def _getGeocentric(self):
        return self._geo
    def _setGeocentric(self,val):
        self._geo = bool(val)
    geocentric = property(_getGeocentric,_setGeocentric,doc="""
    True if the coordinate system origin is the Earth 
    """)
    
    def _getHeliocentric(self):
        return not self._geo
    def _setHeliocentric(self,val):
        self._geo = not bool(val)
    heliocentric = property(_getHeliocentric,_setHeliocentric,doc="""
    True if the coordinate system origin is the Sun 
    """)
    
    @property
    def epoch(self):
        """
        Epoch of observation
        """
        return self._epoch
    

class _LatLongMeta(_CoosysMeta):
    def __init__(cls,name,bases,dct):
        _CoosysMeta.__init__(cls,name,bases,dct)
        if cls._latlongnames_[0] is not None:
            setattr(cls,cls._latlongnames_[0],cls.lat)
            setattr(cls,cls._latlongnames_[0]+'err',cls.laterr)
        if cls._latlongnames_[1] is not None:
            setattr(cls,cls._latlongnames_[1],cls.long)
            setattr(cls,cls._latlongnames_[1]+'err',cls.longerr)

class LatLongCoordinates(CoordinateSystem):
    """
    This object represents an angular location on a unit sphere as represented
    in spherical coordinates with a latitude and longitude.  Subclasses specify 
    details such as transformations or epochs.
    
    Note that converter functions (described in
    :meth:`CoordinateSystem.registerTranform`) for LatLongCoordinates 
    subclasses are epected to have the signature f(inobj,matrix=False). They
    should return a rotation matrix if called with the second argument True. The
    matrix will be applied to cartesian coordinates on the unit sphere to
    perform transformations.
    
    """
    __slots__ = ('_lat','_long','_laterr','_longerr')
    __metaclass__ = _LatLongMeta
    _latlongnames_ = (None,None)
    _longrange_ = None
    
    hubcoosys = None #this is set to EquatorialCoordinatesCIRS below
    """
    The coordinate system that :class:`LatLongCoordinates` try to transform 
    to/from if no direct transformation is available.  If None, it defaults to
    :class:`EquatorialCoordinates`.
    """
    
    def __init__(self,lat=0,long=0,laterr=None,longerr=None):
        """
        See the associated attribute docstrings for the meaning of the inputs.  
        """
        
        self._lat = AngularCoordinate(range=(-90,90))
        self._long = AngularCoordinate(range=self._longrange_)
        
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
        return dict([(k,getattr(k)) for k in LatLongCoordinates.__slots__])
    
    def __setstate__(self,d):
        for k in LatLongCoordinates.__slots__:
            setattr(k,d[k])
        
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
        elif isinstance(val,AngularSeperation):
            self._laterr = val
        else:
            self._laterr = AngularSeperation(val)
    laterr = property(_getLaterr,_setLaterr,doc="""
    Latitude error for this object as a :class:`AngularSeperation` object. May
    be set using any valid input form for :class:`AngularSeperation`.
    """)
    
    def _getLongerr(self):
        return self._longerr
    def _setLongerr(self,val):
        if val is None:
            self._longerr = None
        elif isinstance(val,AngularSeperation):
            self._longerr = val
        else:
            self._longerr = AngularSeperation(val)
    longerr = property(_getLongerr,_setLongerr,doc="""
    Longitude error for this object as a :class:`AngularSeperation` object. May
    be set using any valid input form for :class:`AngularSeperation`.
    """)
    
    def __str__(self):
        lat,long = self._lat.d,self._long.d
        #lat,long = self._lat.getDmsStr(),self._long.getDmsStr()
        return '{0}: {1[0]}={2},{1[1]}={3}'.format(self.__class__.__name__,self._latlongnames_,lat,long)
    
    def getCoordinateString(self,sep=' ',labels=False,canonical=False,hmslong=False):
        coords = []
        if hmslong:
            coords.append(self._long.getHmsStr(canonical=canonical,sign=False))
        else:
            coords.append(self._long.getDmsStr(canonical=canonical,sign=False))
        coords[-1] = self._latlongnames_[1]+'='+coords[-1]
        coords.append(self._lat.getDmsStr(canonical=canonical))
        coords[-1] = self._latlongnames_[0]+'='+coords[-1]
        return sep.join(coords)
    
    def __eq__(self,other):
        if hasattr(other,'_lat') and hasattr(other,'_long'):
            return self._lat==other._lat and self._long==other._long
        else:
            return False
        
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __sub__(self,other):        
        if isinstance(other,self.__class__):
            from math import cos,degrees,acos,asin,sin,sqrt
            
            b1 = self._lat.radians
            b2 = other._lat.radians
            db = abs(b2 - b1)
            dl = abs(other._long.radians - self._long.radians)
            
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
            
            return AngularSeperation(degrees(sep))
            
#            #small angle version
#            from math import cos,degrees,sqrt
#            clat = cos((self._lat.radians+other._lat.radians)/2)
#            dcorrlong = (self._long.radians - other._long.radians)*clat
#            dlat = self._lat.radians-other._lat.radians
#            sep = AngularSeperation(degrees(sqrt(dlat*dlat+dcorrlong*dcorrlong)))
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
        
    def getConversionMatrix(self,tosys):
        """
        Generates and returns the rotation matrix to convert from this object's
        coordinate system into the requested coordinate system class.
        
        :param tosys: The target coordinate system
        :type tosys: A subclass of :class:`LatLongCoordinates`
        
        :returns: 
            A 3x3 rotation matrix or None if the conversion is not matrix-based
        
        :except NotImplementedError: 
            If no conversion exists between these coordinates.
            
        """
        if tosys is self.__class__:
            return np.eye(3).view(np.matrix)
        else:
            try:
                conv = CoordinateSystem._converters[self.__class__][tosys]
                if conv.transtype != 'smatrix':
                    return None
                else:
                    return conv(self)
            except (KeyError,TypeError),e:
                strf = 'cannot generate matrix to transform from {0} to {1}'
                raise NotImplementedError(strf.format(self.__class__.__name__,tosys))
        
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
    
    def convert(self,tosys):
        """
        Converts the coordinate system from it's current system to a new
        :class:`CoordinateSystem` object. For :class:`LatLongCoordinate`
        objects, if a converter is not available, a conversion to
        :attr:`LatLongCoordinates.hubcoosys` (which by default is
        :class:`EquatorialCoordinatesCIRS`) and then to the target system will be
        applied.
        
        :param tosys: The new coordinate system 
        :type tosys: A subclass of :class:`CoordinateSystem`
        
        :except: raises NotImplementedError if converters are not present
        """
        convs = CoordinateSystem._converters[self.__class__]
        hubcoosys = self.hubcoosys
#        if hubcoosys is None:
#            hubcoosys = EquatorialCoordinatesCIRS #this can be reassigned by user
        
        if not issubclass(tosys,CoordinateSystem):
            raise TypeError('Coordinate system to convert to must be a subclass of CoordinateSystem ')
        elif tosys in convs:
            return CoordinateSystem.convert(self,tosys)
        elif hubcoosys in convs:
            try:
                m1 = convs[hubcoosys](self) 
                #calling hubcoosys->target conv w/self technically incorrect, but often gets the right result
                m2 = CoordinateSystem._converters[hubcoosys][tosys](self) 
                return LatLongCoordinates._smatrix(m2*m1,self,tosys)
            except Exception,e:
                niestr = 'could not convert from {0} to {1} by way of {2}'
                niestr = niestr.format(self.__class__,tosys,self.hubcoosys)
                raise NotImplementedError(nies,e)
        else:
            niestr = 'could not convert from {0} to {1}'
            raise NotImplementedError(niestr.format(self.__class__,tosys))
        
    
class EpochalLatLongCoordinates(LatLongCoordinates):
    """
    The origin and orientation of some coordinate systems are tied to the motion
    of the Earth. Hence, a particular location of the equinox as specified by
    the epoch.
    
    Subclasses must implement the :meth:`transformToEpoch` method.
    """
    __slots__ = ('_epoch',)
    
    #: If True, this coordinate system uses Julian Epochs.  Otherwise, Besselian
    julianepoch = True 
    
    def __init__(self,lat=0,long=0,laterr=None,longerr=None,epoch=None):
        LatLongCoordinates.__init__(self,lat,long,laterr,longerr)
        self._epoch = epoch
        
    
    def __getstate__(self):
        d = LatLongCoordinates.__getstate__(self)
        d['_epoch'] = self._epoch
        return d
    
    def __setstate__(self,d):
        LatLongCoordinates.__setstate__(self,d)
        self._epoch = d['_epoch']
    
    def _getEpoch(self):
        return self._epoch
    def _setEpoch(self,val):
        if val is None:
            self._epoch = None
        else:
            if val == 'now':
                from obstools import jd_to_epoch
                val = jd_to_epoch(None)
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
        return '{0}{1}'.format('J' if self.julianepoch else 'B',self._epoch)
    def _setEpochstr(self,val):
        self.epoch = val
    epochstr = property(_getEpochstr,_setEpochstr,doc="""
    A string representation of the current Epoch with a J or B prefixed for
    julian or besselian epochs.
    """)
    
    def _getJdepoch(self):
        from .obstools import epoch_to_jd
        return epoch_to_jd(self._epoch,self.julianepoch)
    def _setJdepoch(self,val):
        from .obstools import jd_to_epoch
        self._epoch = jd_to_epoch(val,self.julianepoch)
    jdepoch = property(_getJdepoch,_setJdepoch,doc="""
    Julian Date of the current epoch.  If set, it is assumed to be a Julian
    Epoch.
    """)
    
    def convert(self,tosys):
        res = LatLongCoordinates.convert(self,tosys)
        if issubclass(tosys,EpochalLatLongCoordinates):
            res._epoch = self._epoch
        return res
    convert.__doc__ = LatLongCoordinates.convert.__doc__
    
    @abstractmethod
    def transformToEpoch(self,newepoch):
        """
        Subclasses should implement this method to transform their coordinates
        to a new epoch.  Subclasses should always call this method *after* the
        transformation is performed
        """
        self._epoch = newepoch
        

class EquatorialCoordinatesBase(EpochalLatLongCoordinates):
    """
    This object represents an angular location on the unit sphere, specified in
    right ascension and declination.  Some of the subclasses are not strictly 
    speaking equatorial, but they are close, or are tied to the equatorial 
    position of a particular epoch.
    
    This is a superclass for all other Equatorial Coordinate systems - particular
    reference systems implement the :meth:`transformToEpoch` method.
    """
    
    __slots__ = ('_dpc',)
    _latlongnames_ = ('dec','ra')
    _longrange_ = (0,360)
    
        
    def __str__(self):
        rastr = self.ra.getHmsStr(canonical=True)
        decstr = self.dec.getDmsStr(canonical=True)
        return '{3}: {0} {1} ({2})'.format(rastr,decstr,self.epoch,self.__class__.__name__)
    
    _dpc = None #default distance should be infinity
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
            except TypeError:
                self._dpc = (float(val),0)
    distancepc = property(_getDistancepc,_setDistancepc,doc="""
    Parallax distance to object in parsecs, or None to assume infinity. Set as
    either a float, a 2-tuple (distance,distance_error), or a no-argument
    callable that returns such a tuple. Getter always returns 2-tuple or None.
    """)
    
    def _getDistanceau(self):
        from .constants import aupercm,cmperpc
        auperpc = cmperpc * aupercm
        if self._dpcs is None:
            return None
        elif callable(self._dpc):
            return self._dpc() * auperpc
        else:
            return self._dpc * auperpc
    def _setDistanceau(self,val):
        from .constants import cmperau,pcpercm
        pcperau = cmperau * pcpercm
            
        if val is None:
            self._dpc = None
        elif callable(val):
            self._dpc = lambda: val()*pcperau
        else:
            try:
                self._dpc = (float(val[0])*pcperau,float(val[1])*pcperau)
            except TypeError:
                self._dpc = (float(val)*pcperau,0)
    distanceau = property(_getDistanceau,_setDistanceau,doc="""
    Parallax distance to object in AU, or None to assume infinity. Set as either
    a float, a 2-tuple (distance,distance_error), or a no-argument callable that
    returns such a tuple. Getter always returns 2-tuple or None.
    """)
    
    def __getstate__(self):
        d = EpochalLatLongCoordinates.__getstate__(self)
        d['_dpc'] = self._dpc
        return d
    
    def __setstate__(self,d):
        EpochalLatLongCoordinates.__setstate__(self,d)
        self._dpc = d['_dpc']
    
    def __init__(self,*args,**kwargs):
        """
        Input for equatorial coordinates. Can follow any of the following forms:
        
        * EquatorialCoordinatesBase()
        * EquatorialCoordinatesBase(:class:`EquatorialCoordinatesBase`)
        * EquatorialCoordinatesBase('rastr decstr')
        * EquatorialCoordinatesBase(ra,dec)
        * EquatorialCoordinatesBase(ra,fdec,raerr,decerr)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch)
        * EquatorialCoordinatesBase(ra,dec,raerr,decerr,epoch,distancepc)
        
        Note that the default epoch is 2000 if not otherwise specified.  To 
        disable epoch tranformations, set the epoch to None.
        
        """
        posargs = {}
        if len(args) == 0:
            pass
        if len(args) == 1:
            if isinstance(args[0],EquatorialCoordinatesBase):
                super(EquatorialCoordinates,self).__init__(args[0])
                self.epoch = args[0].epoch
                return
            elif isinstance(args[0],basestring):
                sargs = args[0].split()
                posargs['ra'] = sargs[0]
                posargs['dec'] = sargs[1]
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
        
        EpochalLatLongCoordinates.__init__(self,kwargs['dec'],kwargs['ra'],kwargs['decerr'],kwargs['raerr'])
        if 'epoch' in kwargs:
            self.epoch = kwargs['epoch']
        if 'distancepc' in kwargs:
            self.distancepc = kwargs['distancepc']
            
    def convert(self,tosys):
        res = EpochalLatLongCoordinates.convert(self,tosys)
        if self._dpc is not None:
            res._dpc = self._dpc
        return res
    convert.__doc__ = EpochalLatLongCoordinates.convert.__doc__
            
            
            
def _precession_matrix_J2000_Capitaine(epoch):
        """
        Computes the precession matrix from J2000 to the given Julian Epoch.
        Expression from from Capitaine et al. 2003 as written in the USNO
        Circular 179.  This should match the IAU 2006 standard from SOFA 
        (although this has not yet been tested)
        """
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
    
    lines = [l for l in _get_package_data(datafn).split('\n') if not l.startswith('#') if not l.strip()=='']
    
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
    from obstools import epoch_to_jd
    
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
    from obstools import epoch_to_jd
    if asepoch:
        jd = epoch_to_jd(intime)
    else:
        jd = intime
    epsa = np.radians(obliquity(jd,2000))
    t = (jd-epoch_to_jd(2000))/36525
    
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
    #TODO: implement higher precsion 2006/2000A model if requested/needed
    epsa,dpsi,deps = _nutation_components2000B(epoch) #all in radians
    
    return rotation_matrix(-(epsa + deps),'x',False) *\
           rotation_matrix(-dpsi,'z',False) *\
           rotation_matrix(epsa,'x',False)
           

def _load_CIO_locator_data(datafn):
    """
    Loads CIO locator series terms from saved data files.
    
    returns polycoeffs,termsarr (starting with 0th)
    """
    
    lines = [l for l in _get_package_data(datafn).split('\n') if not l.startswith('#') if not l.strip()=='']
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


class EquatorialCoordinatesCIRS(EquatorialCoordinatesBase):
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
    should be good to ~1 mas. If abberration or annual parallax corrections are
    necessary, convert to :class:`ICRSCoordinates`, change the epoch, and then
    convert back to :class:`EquatorialCoordinatesCIRS`.
    
    To convert from these coordinates to :class:`HorizontalCoordinates`
    appropriate for observed coordinates, site information is necessary. Hence,
    the transformations from equatorial to horizontal coordinates are performed
    by the :class:`~astropysics.obstools.Site` class in the
    :mod:`~astropysics.obstools` module, and attempting to directly convert will
    raise an :exc:`TypeError`.
    """
    
    #if True, nutation will be used, if False, it will be ignored
    _nuton = True
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch using the
        IAU 2000 precessions from Capitaine, N. et al. 2003 as written in the
        USNO Circular 179.
        """
        if self.epoch is not None and newepoch is not None:
            if EquatorialCoordinatesCIRS._nuton:
                #convert from current to J2000
                B = _precession_matrix_J2000_Capitaine(self.epoch).T #T==inv; real unitary
                Bn = _nutation_matrix(self.epoch).T 
                #convert to new epoch
                A = _precession_matrix_J2000_Capitaine(newepoch)
                An = _nutation_matrix(newepoch) 
                M = An*A*B*Bn
            else:
                #convert from current to J2000
                B = _precession_matrix_J2000_Capitaine(self.epoch).T #T==inv; real unitary
                #convert to new epoch
                A = _precession_matrix_J2000_Capitaine(newepoch)
                M = A*B
            self.matrixRotate(M)
            
        #this sets the epoch
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
        
    @staticmethod
    def _CIOLocator(epoch):
        """
        Returns the CIO locator s for the provided epoch. s is the difference in
        RA between the GCRS and CIP points for the ascending node of the CIP 
        equator.
        """
        from .obstools import jd2000,epoch_to_jd
        from .constants import asecperrad
        
        from .ephems import _mean_anomaly_of_moon,_mean_anomaly_of_sun,\
                            _mean_long_of_moon_minus_ascnode,_long_earth,\
                            _mean_elongation_of_moon_from_sun,_long_venus,\
                            _mean_long_asc_node_moon,_long_prec
        
        #first need to find x and y for the CIP, as s+XY/2 is needed
        B = ICRSCoordinates.frameBiasJ2000
        P = _precession_matrix_J2000_Capitaine(epoch)
        N = _nutation_matrix(epoch)
        
        #B*P*N takes GCRS to true, so CIP is bottom row
        x,y,z = (B*P*N).A[2]
        
        T = (epoch_to_jd(epoch) - jd2000)/36525
        
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
        
#set default hub coordinate system to EquatorialCoordinatesCIRS
LatLongCoordinates.hubcoosys = EquatorialCoordinatesCIRS
            
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
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        """
        Transforms these :class:`EquatorialCoordinates` to a new epoch using the
        IAU 2000 precessions from Capitaine, N. et al. 2003 as written in the
        USNO Circular 179.
        """
        if self.epoch is not None and newepoch is not None:
            #convert from current to J2000
            B = _precession_matrix_J2000_Capitaine(self.epoch).T #T==inv; real unitary
            #convert to new epoch
            A = _precession_matrix_J2000_Capitaine(newepoch)
            self.matrixRotate(A*B)
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
          
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesCIRS,transtype='smatrix')
    def _toCIRS(eqsys):
        raise NotImplementedError
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesCIRS,'self',transtype='smatrix')
    def _fromCIRS(cirssys):
        raise NotImplementedError
            

class ICRSCoordinates(EquatorialCoordinatesBase):
    """
    Equatorial Coordinates tied to the International Celestial Reference System
    (ICRS). Strictly speaking this is not an Equatorial system, as it is an
    inertial frame that only aligns with Earth's equator at J2000, but it is
    nearly an equatorial system at J2000.
    
    .. note::
        Strictly speaking this class is actually the Barycentric Celestial
        Reference System (BCRS), as any space/proper motions are computed before 
        this object (generally in :class:`EphemerisObject`). But BCRS's orientation 
        is tied to ICRS (it can be considered , so it is called ICRS here for
        clarity.
        
    .. warning:: 
        Abberation and annual parallax are not yet included in transformations!
        Also, all general relativistic effects are currently ignored.
        
    """
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        """
        ICRS is an inertial frame, so no transformation is necessary
        """
        EpochalLatLongCoordinates.transformToEpoch(self,newepoch)
    
    
    #conversion methods and attributes are below
    
    
    def __makeFrameBias():
        #see USNO circular 179 for frame bias -- all in milliarcsec
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
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox,transtype='smatrix')
    def _toEqE(icrsc):
        B = ICRSCoordinates.frameBiasJ2000
        if icrsc.epoch is None:
            return B
        else:
            P = _precession_matrix_J2000_Capitaine(icrsc.epoch)
            return P*B #no nutation implemented for equinox-based
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesCIRS,transtype='smatrix')
    def _toEqC(icrsc):
        B = ICRSCoordinates.frameBiasJ2000
        if icrsc.epoch is None:
            return B
        else:
            P = _precession_matrix_J2000_Capitaine(icrsc.epoch)
            if EquatorialCoordinatesCIRS._nuton:
                N = _nutation_matrix(icrsc.epoch)
                return N*P*B
            else:
                return P*B
                
class ITRSCoordinates(EquatorialCoordinatesBase):
    """
    ???
        
    """
    
    __slots__ = tuple()
    
    def transformToEpoch(self,newepoch):
        #implement polar motion
        raise NotImplementedError
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesCIRS,'self',transtype='smatrix')
    def _fromEqC(eqc):
        raise NotImplementedError
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self',transtype='smatrix')
    def _fromEqE(eqe):
        from .obstools import greenwich_sidereal_time
        raise NotImplementedError    
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesCIRS,transtype='smatrix')
    def _toEqC(itrsc):
        #really we want inverse, but rotations are unitary -> inv==transpose
        #we provide eqc in the call because the epoch is needed
        return ICRSCoordinates._fromEqC(itrsc).T 
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox,transtype='smatrix')
    def _toEqE(itrsc):
        #really we want inverse, but rotations are unitary -> inv==transpose
        #we provide eqc in the call because the epoch is needed
        return ICRSCoordinates._fromEqE(itrsc).T 
            
class FK5Coordinates(EquatorialCoordinatesEquinox):
    """
    Equatorial Coordinates fixed to the FK5 reference system. 
    """
    
    __slots__ = tuple()
    
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
    
    @staticmethod
    def _precessionMatrixJ(epoch1,epoch2):
        """
        Computes the precession matrix from one Julian epoch to another
        """
        from .obstools import epoch_to_jd,jd2000
        
        jd1 = epoch_to_jd(epoch1)
        jd2 = epoch_to_jd(epoch2)
        
        T = (jd1 - jd2000)/36525
        dt = (jd2 - jd1)/36525
        
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
    
class FK4Coordinates(EquatorialCoordinatesEquinox):
    """
    Equatorial Coordinates fixed to the FK4 reference system.  Note that this 
    implementation does *not* correct for the elliptic terms of aberration
    as of yet.
    
    Epoch is Besselian.
    """
    
    __slots__ = tuple()
    
    julianepoch = False
    
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
        #tropical years
        t1 = (epoch1-1850.0)/1000.0    
        t2 = (epoch2-1850.0)/1000.0
        dt = t2 - t1
        
        zeta1 = 23035.545 + t1*139.720+0.060*t1*t1
        zeta2 = 30.240 - 0.27*t1
        zeta3 = 17.995
        pzeta = (zeta3,zeta2,zeta1,0)
        zeta = np.polyval(pzeta,dt)
        
        z1 = 23035.545 + t1*139.720 + 0.060*t1*t1
        z2 = 109.480 + 0.39*t1
        z3 = 18.325
        pz = (z3,z2,z1,0)
        z = np.polyval(pz,dt)
        
        theta1 = 20051.12 - 85.29*t1 - 0.37*t1*t1
        theta2 = -42.65 - 0.37*t1
        theta3 = -41.8
        ptheta = (theta3,theta2,theta1,0)
        theta = np.polyval(ptheta,dt)
        
        return rotation_matrix(-z,'z') *\
               rotation_matrix(theta,'y') *\
               rotation_matrix(-zeta,'z')
        
        
class EclipticCoordinatesCIRS(EpochalLatLongCoordinates):
    """
    Ecliptic Coordinates (beta, lambda) such that the fundamental plane passes
    through the ecliptic at the current epoch.
    
    Note that because the concept of the ecliptic can be complicated or even
    ill-defined, ecliptic coordinates is astropysics are simply defined as
    tied to a particular set of equatorial coordinates with a given obliquity
    model.  For :class:`EclipticCoordinatesCIRS`, the equatorial 
    coordinates are :class:`EquatorialCoordinatesCIRS` with obliquity
    given by the IAU 2006 obliquity model (see :func:`obliquity`)
    """
    
    __slots__ = ()
    _latlongnames_ = ('beta','lamb')
    _longrange_ = (0,360)
    
    obliqyear = 2006
    
    def __init__(self,beta=0,lamb=0,betaerr=None,lamberr=None,epoch=2000):
        EpochalLatLongCoordinates(beta,lamb,betaerr,lamberr,epoch)
        
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesCIRS,transtype='smatrix')
    def _toEq(eclsc):
        return rotation_matrix(-obliquity(eclsc.jdepoch,EclipticCoordinatesCIRS.obliqyear),'x')
        
    @CoordinateSystem.registerTransform(EquatorialCoordinatesCIRS,'self',transtype='smatrix')
    def _fromEq(eqc):
        return rotation_matrix(obliquity(eclsc.jdepoch,EclipticCoordinatesCIRS.obliqyear),'x')
    
    def transformToEpoch(self,newepoch):
        if self.epoch is not None and newepoch is not None:
            eqc = self.convert(EquatorialCoordinatesCIRS)
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
    ill-defined, ecliptic coordinates is astropysics are simply defined as
    tied to a particular set of equatorial coordinates with a given obliquity
    model.  For :class:`EclipticCoordinatesEquinox`, the equatorial 
    coordinates are :class:`EquatorialCoordinatesEquinox` with obliquity
    given by the IAU 1980 obliquity model (see :func:`obliquity`)
    """
    
    __slots__ = ()
    _latlongnames_ = ('beta','lamb')
    _longrange_ = (0,360)
    
    obliqyear = 1980
    
    def __init__(self,beta=0,lamb=0,betaerr=None,lamberr=None,epoch=2000):
        EpochalLatLongCoordinates(beta,lamb,betaerr,lamberr,epoch)
        
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox,transtype='smatrix')
    def _toEq(eclsc):
        return rotation_matrix(-obliquity(eclsc.jdepoch,EclipticCoordinatesEquinox.obliqyear),'x')
        
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self',transtype='smatrix')
    def _fromEq(eqc):
        return rotation_matrix(obliquity(eclsc.jdepoch,EclipticCoordinatesEquinox.obliqyear),'x')
        
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
    
class GalacticCoordinates(LatLongCoordinates):
    __slots__ = tuple()
    _latlongnames_ = ('b','l')
    _longrange_ = (0,360)
    
    _ngp_J2000 = FK5Coordinates(192.859508, 27.128336,epoch=2000)
    _long0_J2000 = AngularCoordinate(122.932)
    _ngp_B1950 = FK4Coordinates(192.25, 27.4,epoch=1950)
    _long0_B1950 = AngularCoordinate(123)
    
    @CoordinateSystem.registerTransform(FK5Coordinates,'self',transtype='smatrix')
    def _fromFK5(fk5coords):
        return rotation_matrix(180 - GalacticCoordinates._long0_J2000.d,'z') *\
               rotation_matrix(90 - GalacticCoordinates._ngp_J2000.dec.d,'y') *\
               rotation_matrix(GalacticCoordinates._ngp_J2000.ra.d,'z')
    
    @CoordinateSystem.registerTransform('self',FK5Coordinates,transtype='smatrix')
    def _toFK5(galcoords):
        return GalacticCoordinates._fromFK5(galcoords).T
    
    @CoordinateSystem.registerTransform(FK4Coordinates,'self',transtype='smatrix')
    def _fromFK4(fk4coords):
        return rotation_matrix(180 - GalacticCoordinates._long0_B1950.d,'z') *\
               rotation_matrix(90 - GalacticCoordinates._ngp_B1950.dec.d,'y') *\
               rotation_matrix(GalacticCoordinates._ngp_B1950.ra.d,'z')
    
    @CoordinateSystem.registerTransform('self',FK4Coordinates,transtype='smatrix')
    def _toFK4(galcoords):
        return GalacticCoordinates._fromFK4(galcoords).T
        
class SupergalacticCoordinates(LatLongCoordinates):   
    __slots__ = tuple()
    _latlongnames_ = ('sgb','sgl')
    _longrange_ = (0,360)
    
    _nsgp_gal = FK4Coordinates(6.32,47.37,epoch=1950)
    _sglong0_gal = 137.37
    _nsgp_J2000 = FK5Coordinates(283.75420420,15.70894043,epoch=2000)
    _sglong0_J2000 = 42.30997710
    
    @CoordinateSystem.registerTransform('self',GalacticCoordinates,transtype='smatrix')
    def _toGal(sgalcoords):
        return SupergalacticCoordinates._fromGal(sgalcoords).T
    
    @CoordinateSystem.registerTransform(GalacticCoordinates,'self',transtype='smatrix')
    def _fromGal(galcoords):
        latang = SupergalacticCoordinates._nsgp_gal.lat.d
        longang = SupergalacticCoordinates._nsgp_gal.long.d
        return rotation_matrix(90-latang,'x')*rotation_matrix(90+longang,'z')
    
class HorizontalCoordinates(LatLongCoordinates):
    """
    This object represents an angular location on the unit sphere, with the 
    north pole of the coordinate position fixed to the local zenith
    
    To convert from other :class:`Coordinate` types to horizontal positions, see 
    :class:`astropysics.obstools.Site`, as site information is required for
    these corrections
    """  
    __slots__ = tuple()
    _latlongnames_ = ('alt','az')
    _longrange_ = (0,360)
    
    @CoordinateSystem.registerTransform(EquatorialCoordinatesEquinox,'self')
    @CoordinateSystem.registerTransform(EquatorialCoordinatesCIRS,'self')
    def _toHoriz(incoosys=None):
        raise TypeError('use astropysics.obstools.Site methods to transform celestial to terrestrial coordinates')
    
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesEquinox)
    @CoordinateSystem.registerTransform('self',EquatorialCoordinatesCIRS)
    def _fromHoriz(incoosys=None):
        raise TypeError('use astropysics.obstools.Site methods to transform terrestrial to celestial coordinates')


#rotation_matrix function should really be defined here, but it is above 
#because of it's use for generating class/module variables


    
def obliquity(jd,algorithm=2006):
    """
    Computes the obliquity of the Earth at the requested Julian Date. 
    
    :param jd: julian date at which to compute obliquity
    :type jd: scalar or array-like
    :param algorithm: 
        Year of algorithm based on IAU adoption. Can be 2006, 2000 or 1980. The
        2006 algorithm is mentioned in Circular 179, but the canonical reference
        for the IAU adoption is apparently Hilton et al. 06 is composed of the
        1980 algorithm with a precession-rate correction due to the 2000
        precession models, and a description of the 1980 algorithm can be found
        in the Explanatory Supplement to the Astronomical Almanac.
    
    :type algorithm: int
    
    :returns: mean obliquity in degrees (or array of obliquities)
    
    .. seealso::
        
        * Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351. 2000
        * USNO Circular 179
        * Explanatory Supplement to the Astronomical Almanac: P. Kenneth
          Seidelmann (ed), University Science Books (1992).
    """
    from .obstools import jd2000
    
    T = (jd-jd2000)/36525.0
    
    if algorithm==2006:
        p = (-0.0000000434,-0.000000576,0.00200340,-0.0001831,-46.836769,84381.406)
        corr = 0
    elif algorithm==2000:
        p = (0.001813,-0.00059,-46.8150,84381.448)
        corr = -0.02524*T
    elif algorithm==1980:
        p = (0.001813,-0.00059,-46.8150,84381.448)
        corr = 0
    else:
        raise ValueError('invalid algorithm year for computing obliquity')
        
    return (np.polyval(p,T)+corr)/3600.

def earth_rotation_angle(jd,degrees=True):
    """
    Earth Rotation Angle (ERA) for a given Julian Date.
    
    :param jd: The Julian Date or a sequence of JDs
    :type jd: scalar or array-like
    :param degrees: 
        If True, the ERA is returned in degrees, if None, 1=full rotation.  
        Otherwise, radians.
    :type degrees: bool or None
    
    :returns: ERA or an array of angles (if `jd` is an array) 
    
    """
    d = jd - 2451545.0 #days since 2000
    res = (0.7790572732640 + 0.00273781191135448*d + (d%1.0))%1.0
    
    if degrees is None:
        return res
    elif degrees:    
        return res*360
    else:
        return res*2*pi
        
def greenwich_sidereal_time(jd,apparent=True):
    """
    Computes the Greenwich Sidereal Time for a given Julian Date.
    
    :param jd: The Julian Date or a sequence of JDs
    :type jd: scalar or array-like
    :param apparent: 
        If True, the Greenwich Apparent Sidereal Time (GAST) is returned,
        computed from the IAU 2000B nutation model. In the special case that
        'simple' is given, a faster (but much lower precision) nutation model
        will be used. If False, the Greenwich Mean Sidereal Time (GMST) is
        returned, instead.
    :type apparent: 
    
    :returns: GMST or GAST in hours or an array of times (if `jd` is an array) 
        
    .. seealso:: 
        :func:`equation_of_the_equinoxes`
        USNO Circular 179 and http://aa.usno.navy.mil/faq/docs/GAST.php
    
    
    """
    from .constants import asecperrad
    
    era = earth_rotation_angle(jd,False) #in radians
    
    t = (jd - 2451545.0)/36525
    
    gmst = era + (0.014506 + 4612.156534*t + 1.3915817*t**2 - 0.00000044*t**3 -\
            0.000029956*t**4 - 0.0000000368*t**5)/asecperrad
            
    if apparent:
        if apparent == 'simple':
            eps =  np.radians(23.4393 - 0.0000004*d) #obliquity
            L = np.radians(280.47 + 0.98565*d) #mean longitude of the sun
            omega = np.radians(125.04 - 0.052954*d) #longitude of ascending node of moon
            dpsi = -0.000319*np.sin(omega) - 0.000024*np.sin(2*L) #nutation longitude
            coor = 0
        else:
            from .coords import _nutation_components2000B
            eps,dpsi,deps = _nutation_components2000B(jd,False)
            dpsi = dpsi
            raise  NotImplementedError('need to implement complementary terms for equation of the equinoxes from iauEect00 0 use "simple" for now')
            coor = 0
        return ((gmst + dpsi*np.cos(eps))*12/pi + coor)%24
    else:
        return (gmst*12/pi)%24
    
#    #previous algorithm described on USNO web site http://aa.usno.navy.mil/faq/docs/GAST.php
#    jd0 = np.round(jd-.5)+.5
#    h = (jd - jd0) * 24.0
#    d = jd - 2451545.0
#    d0 = jd0 - 2451545.0
#    t = d/36525
    
#    #mean sidereal time @ greenwich
#    gmst = 6.697374558 + 0.06570982441908*d0 + 0.000026*t**2 + 1.00273790935*h
#           #- 1.72e-9*t**3 #left off as precision to t^3 is unneeded
   
#    if apparent:
#        eps =  np.radians(23.4393 - 0.0000004*d) #obliquity
#        L = np.radians(280.47 + 0.98565*d) #mean longitude of the sun
#        omega = np.radians(125.04 - 0.052954*d) #longitude of ascending node of moon
#        dpsi = -0.000319*np.sin(omega) - 0.000024*np.sin(2*L) #nutation longitude
#        return (gmst + dpsi*np.cos(eps))%24.0
#    else:
#        return gmst%24.0 

def equation_of_the_equinoxes(jd):
    """
    Computes equation of the equinoxes GAST-GMST.
    
    :param jd: The Julian Date or a sequence of JDs.
    :type jd: scalar or array-like
    
    :returns: the equation of the equinoxes for the provided date in hours.
    
    """
    return greenwich_sidereal_time(jd,True) - greenwich_sidereal_time(jd,False)

def equation_of_the_origins(jd):
    """
    Computes the equation of the origins ERA - GAST
    (ERA = Earth Rotation Angle, GAST = Greenwich Apparent Sidereal Time) 
    
    :param jd: The Julian Date or a sequence of JDs.
    :type jd: scalar or array-like
    
    :returns: the equation of the origins for the provided date in hours.
    
    """
    return earth_rotation_angle(jd,None)*24. - greenwich_sidereal_time(jd,True)


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

    

#<--------------------Functional coordinate transforms------------------------->
def cartesian_to_polar(x,y,degrees=False):
    """
    Converts arrays in 2D rectangular Cartesian coordinates to polar
    coordinates.
    
    :param x: First cartesian coordinate
    :type x: :class:`numpy.ndarray`
    :param y: Second cartesian coordinate
    :type y: :class:`numpy.ndarray`
    :param degrees: 
        If True, the output theta angle will be in degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: 
        (r,theta) where theta is measured from the +x axis increasing towards
        the +y axis
    """
    r = (x*x+y*y)**0.5
    t = np.arctan2(y,x)
    if degrees:
        t = np.degrees(t)
    
    return r,t

def polar_to_cartesian(r,t,degrees=False):
    """
    Converts arrays in 2D polar coordinates to rectangular cartesian
    coordinates.
    
    Note that the spherical coordinates are in *physicist* convention such that
    (1,0,pi/2) is x-axis.
    
    :param r: Radial coordinate
    :type r: :class:`numpy.ndarray`
    :param t: Azimuthal angle from +x-axis increasing towards +y-axis
    :type t: :class:`numpy.ndarray`
    :param degrees: 
        If True, the input angles will be in degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: arrays (x,y)
    """
    if degrees:
        t=np.radians(t)
        
    return r*np.cos(t),r*np.sin(t)

def cartesian_to_spherical(x,y,z,degrees=False):
    """
    Converts three arrays in 3D rectangular cartesian coordinates to
    spherical polar coordinates.
    
    Note that the spherical coordinates are in *physicist* convention such that
    (1,0,pi/2) is x-axis.
    
    :param x: First cartesian coordinate
    :type x: :class:`numpy.ndarray`
    :param y: Second cartesian coordinate
    :type y: :class:`numpy.ndarray`
    :param z: Third cartesian coordinate
    :type z: :class:`numpy.ndarray`
    :param degrees: 
        If True, the output theta angle will be in degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: arrays (r,theta,phi) 
    """
    xsq,ysq,zsq=x*x,y*y,z*z
    r=(xsq+ysq+zsq)**0.5
    #t=np.arccos(z,r) #TODO:check to make even more efficient
    t=np.arctan2((xsq+ysq)**0.5,z)
    p=np.arctan2(y,x)
    if degrees:
        t,p=np.degrees(t),np.degrees(p)
    return r,t,p

def spherical_to_cartesian(r,t,p,degrees=False):
    """
    Converts arrays in 3D spherical polar coordinates to rectangular cartesian
    coordinates.
    
    Note that the spherical coordinates are in *physicist* convention such that
    (1,0,pi/2) is x-axis.
    
    :param r: Radial coordinate
    :type r: :class:`numpy.ndarray`
    :param t: Colatitude (angle from z-axis)
    :type t: :class:`numpy.ndarray`
    :param p: Azimuthal angle from +x-axis increasing towards +y-axis
    :type p: :class:`numpy.ndarray`
    :param degrees: 
        If True, the input angles will be in degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: arrays (x,y,z)
    """
    if degrees:
        t,p=np.radians(t),np.radians(p)
    x=r*np.sin(t)*np.cos(p)
    y=r*np.sin(t)*np.sin(p)
    z=r*np.cos(t)
    
    return x,y,z

def latitude_to_colatitude(lat,degrees=False):
    """
    converts from latitude  (i.e. 0 at the equator) to colatitude/inclination 
    (i.e. "theta" in physicist convention).
    """
    if degrees:
        return 90 - lat
    else:
        return pi/2 - lat

def colatitude_to_latitude(theta,degrees=False):
    """
    Converts from colatitude/inclination (i.e. "theta" in physicist convention) 
    to latitude (i.e. 0 at the equator).
    
    :param theta: input colatitude
    :type theta: float or array-like
    :param degrees: 
        If True, the input is interpreted as degrees, otherwise radians.
    :type degrees: bool
    
    :returns: latitude
    
    """
    if degrees:
        return 90 - theta
    else:
        return pi/2 - theta

def cartesian_to_cylindrical(x,y,z,degrees=False):
    """
    Converts three arrays in 3D rectangular Cartesian coordinates to cylindrical
    polar coordinates.
    
    :param x: x cartesian coordinate
    :type x: float or array-like
    :param y: y cartesian coordinate
    :type y: float or array-like
    :param z: z cartesian coordinate
    :type z: float or array-like
    :param degrees: 
        If True, the output angles will be in degrees, otherwise radians.
    :type degrees: bool
    
    :returns: 
        Cylindrical coordinates as a (rho,theta,z) tuple (theta increasing from
        +x to +y, 0 at x-axis).
    """
    s,t = cartesian_to_polar(x,y)
    return s,t,z
    
def cylindrical_to_cartesian(s,t,z,degrees=False):
    """
    Converts three arrays in cylindrical polar coordinates to 3D rectangular
    Cartesian coordinates.
    
    :param s: radial polar coordinate
    :type s: float or array-like
    :param t: polar angle (increasing from +x to +y, 0 at x-axis)
    :type t: float or array-like
    :param z: z coordinate
    :type z: float or array-like
    :param degrees: 
        If True, the output angles will be in degrees, otherwise radians.
    :type degrees: bool
    
    :returns: Cartesian coordinates as an (x,y,z) tuple.
    """
    x,y = polar_to_cartesian(s,t,degrees)
    return x,y,z
    
def offset_proj_sep(rx,ty,pz,offset,spherical=False):
    """
    computes the projected seperation for a list of points in galacto-centric
    coordinates as seen from a point offset (an [[x,y,z]] 2-sequence)
    
    spherical determines if the inputs are spherical coords or cartesian.  If it
    is 'degrees', spherical coordinates will be used, converting from degrees to
    radians
    """
    if spherical is 'degrees':
        x,y,z=spherical_to_cartesian(rx,ty,pz,True)
    elif spherical:
        x,y,z=spherical_to_cartesian(rx,ty,pz,False)
    else:
        x,y,z=rx,ty,pz
    
    offset=np.array(offset)
    if offset.shape[1]!=3 or len(offset.shape)!=2:
        raise ValueError('offset not a sequnce of 3-sequence')
    
    ohat=(offset.T*np.sum(offset*offset,1)**-0.5)
    return np.array(np.matrix(np.c_[x,y,z])*np.matrix(ohat))


def sky_sep_to_3d_sep(pos1,pos2,d1,d2):
    """
    Compute the full 3D seperation between two objects at distances `d1` and
    `d2` and angular positions `pos1` and `pos2` (:class:`LatLongCoordinates`
    objects, or an argument that will be used to generate a
    :class:`EquatorialCoordinates` object)
    
    :param pos1: on-sky position of first object
    :type pos1: :class:`LatLongCoordinates` or initializer
    :param pos2: on-sky position of second object
    :type pos2: :class:`LatLongCoordinates` or initializer
    :param d1: distance to first object
    :type d1: scalar
    :param d2: distance to second object
    :type d2: scalar
    
    .. testsetup::
    
        from astropysics.coords import sky_sep_to_3d_sep
    
    .. doctest::
    
        >>> p1 = LatLongCoordinates(0,0)
        >>> p2 = LatLongCoordinates(0,10)
        >>> '%.10f'%sky_sep_to_3d_sep(p1,p2,20,25)
        '6.3397355613'
        >>> '%.10f'%sky_sep_to_3d_sep('0h0m0s +0:0:0','10:20:30 +0:0:0',1,2)
        '2.9375007333'
        
    """    
    if not isinstance(pos1,LatLongCoordinates):
        pos1 = EquatorialCoordinates(pos1)
    if not isinstance(pos2,LatLongCoordinates):
        pos2 = EquatorialCoordinates(pos2)
        
    return (pos1-pos2).seperation3d(d1,d2)

def radec_str_to_decimal(ra,dec):
    if isinstance(ra,basestring):
        if not isinstance(dec,basestring):
            raise ValueError('either both ra and dec must be a strings or neither')
        
        ra = AngularCoordinate(ra,sghms=True).d
        dec = AngularCoordinate(dec,sghms=False).d
    else:
        if isinstance(dec,basestring):
            raise ValueError('either both ra and dec must be a strings or neither')   
        if len(ra) != len(dec):
            raise ValueError("length of ra and dec don't match")
        
        ras,decs=[],[]
        for r,d in zip(ra,dec):
            ras.append(AngularCoordinate(r,sghms=True).d)
            decs.append(AngularCoordinate(d,sghms=False).d)
        ra,dec = ras,decs
    return ra,dec

def match_coords(a1,b1,a2,b2,eps=1,multi=False):
    """
    Match one coordinate array to another within a specified tolerance. Distance
    is determined by the cartesian distance between the two arrays added in 
    quadrature.
    
    :param a1: the first coordinate for the first set of coordinates
    :type a1: 1D :class:`numpy.ndarray`
    :param a2: the second coordinate for the first set of coordinates
    :type a2: 1D :class:`numpy.ndarray`
    :param a1: the first coordinate for the second set of coordinates
    :type a1: 1D :class:`numpy.ndarray`
    :param a2: the second coordinate for the second set of coordinates
    :type a2: 1D :class:`numpy.ndarray`
    :param multi:
        Determines behavior if more than one coordinate pair matches.  Can be:
        
        * True: raise an exception if more than one match is found
        * 'warn': a warning will be issued if more than one match is found
        * 'print': a statement will be printed  if more than one match is found
        * 'full': the 2D array with matches as booleans along the axes will be returned
        * 'count': a count of matches will be returned instead of a mask
        * 'index': a list of match indecies will be returned instead of a mask
        * False: do nothing - just return if something matched
    
    :returns: 
        Tuple (mask of matches for array 1, mask of matches for array 2) or as
        described in the corresponding `multi` parameter value.
    
    **Examples**
    
    .. testsetup::
    
        from astropysics.coords import match_coords
        from numpy import array
    
    .. doctest::
    
        >>> ra1 = array([1,2,3,4])
        >>> dec1 = array([0,0,0,0])
        >>> ra2 = array([7,6,5,4])
        >>> dec2 = array([.5,.5,.5,.5])
        >>> match_coords(ra1,dec1,ra2,dec2,1)
        (array([False, False, False,  True], dtype=bool), array([False, False, 
        False,  True], dtype=bool))
    """
    def find_sep(A,B):
        At = np.tile(A,(len(B),1))
        Bt = np.tile(B,(len(A),1))
        return At.T-Bt
    sep1=find_sep(a1,a2)
    sep2=find_sep(b1,b2)
    
    matches = (sep1*sep1+sep2*sep2)**0.5 < eps
    if multi:
        if multi == 'full':
            return matches.T
        elif multi == 'count':
            return np.sum(np.any(matches,axis=1)),np.sum(np.any(matches,axis=0)) 
        elif multi == 'index':
            return np.where(matches)
        elif multi == 'warn':
            s1,s2 = np.sum(matches,axis=1),np.sum(matches,axis=0) 
            from warnings import warn
            
            for i in np.where(s1>1)[0]:
                warn('1st index %i has %i matches!'%(i,s1[i]))
            for j in np.where(s2>1)[0]:
                warn('2nd index %i has %i matches!'%(j,s2[j]))
            return s1>0,s2>0
        elif multi == 'print':
            s1,s2 = np.sum(matches,axis=1),np.sum(matches,axis=0) 
            
            for i in np.where(s1>1)[0]:
                print '1st index',i,'has',s1[i],'matches!'
            for j in np.where(s2>1)[0]:
                print '2nd index',j,'has',s2[j],'matches!'
            return s1>0,s2>0
        else:
            raise ValueError('unrecognized multi mode')
    else:
        return np.any(matches,axis=1),np.any(matches,axis=0)
    
def seperation_matrix(v,w=None,tri=False):
    """
    This function takes a n(x?x?x?) array and produces an array given by
    A[i][j] = v[i]-v[j]. if w is not None, it produces A[i][j] = v[i]-w[j]
    
    If the input has more than 1 dimension, the first is assumed to be the 
    one to expand
    
    If tri is True, the lower triangular part of the matrix is set to 0
    (this is really only useful if w is None)
    """
    if w is None:
        w = v
    
    shape1 = list(v.shape)
    shape1.insert(1,1)
    shape2 = list(w.shape)
    shape2.insert(0,1)
    
    A = v.reshape(shape1)-w.reshape(shape2)
    
    if tri:
        return np.triu(A)
    else:
        return A
    

#<--------------------Cosmological distances and conversions------------------->
def cosmo_z_to_dist(z,zerr=None,disttype=0,inttol=1e-6,normed=False,intkwargs={}):
    """
    Calculates the cosmolgical distance to some object given a redshift. Note
    that this uses H0,omegaM,omegaL, and omegaR from the current
    :class:`astropyscs.constants.Cosmology` -- if any of those do not exist in
    the current cosmology this will fail.
    
    The distance type can be one of the following:
    
    * 'comoving'(0) : comoving distance (in Mpc)
    * 'luminosity'(1) : luminosity distance (in Mpc)
    * 'angular'(2) : angular diameter distance (in Mpc)
    * 'lookback'(3) : lookback time (in Gyr)
    * 'distmod'(4) : distance modulus
    
    :param z: 
        The redshift at which to compute the distance, or None to compute the
        maximum value for this distance (for luminosity and distmod this is
        infinite)
    :type z: array, scalar, or None
    :param zerr: Symmetric error in redshift
    :type zerr: array, scalar, or None
    :param disttype:
        The type of distance to compute -- may be any of the types described
        above.
    :type disttype: A string or int
    :param inttol: fractional precision of the output (used in integrals)
    :type inttol: A float<1
    :param normed: 
        If True, normalize output by result for `z` == None.  If a scalar, 
        normalize by the distance at that redshift. If False, no normalization.
    :type normed: boolean
    :param intkwargs: keywords for integrals (see :mod:`scipy.integrate`)
    :type intkwargs: a dictionary   
    
    
    :returns: 
        Distance of type selected by `disttype` in above units or normalized as
        controlled by `normed` parameter. If `zerr` is not None, the output is
        (z,zupper,zlower), otherwise just z.
        
    **Examples**
    
    In these examples we are assuming the WMAP7 BAOH0 cosmological parameters.   
     
    .. testsetup::
    
        from astropysics.constants import choose_cosmology
        choose_cosmology('wmap7baoh0')
        from astropysics.coords import cosmo_z_to_dist
    
    .. doctest::
    
        >>> '%.6f'%cosmo_z_to_dist(0.03)
        '126.964723'
        >>> '%.6f'%cosmo_z_to_dist(0.2)
        '815.469170'
        >>> '%.6f'%cosmo_z_to_dist(0.2,disttype=1)
        '978.563004'
        >>> '%.6f'%cosmo_z_to_dist(0.2,disttype='luminosity')
        '978.563004'
        >>> '%.6f'%cosmo_z_to_dist(0.2,disttype='angular')
        '679.557642'
        >>> '%.3f'%cosmo_z_to_dist(1,disttype='lookback')
        '7.789'
        >>> '%.2f'%cosmo_z_to_dist(0.5,disttype='distmod')
        '42.27'
        >>> '%.6f'%cosmo_z_to_dist(0.2,disttype='angular',normed=True)
        '0.382326'
        >>> '%.6f'%cosmo_z_to_dist(0.8,disttype='angular',normed=True)
        '0.879027'
        >>> '%.6f'%cosmo_z_to_dist(1.64,disttype='angular',normed=True)
        '1.000000'
        >>> '%.6f'%cosmo_z_to_dist(2.5,disttype='angular',normed=True)
        '0.956971'
        
    """
    from operator import isSequenceType
    from scipy.integrate import quad as integrate
    from numpy import array,vectorize,abs,isscalar
    
    from .constants import H0,omegaM,omegaL,omegaR,c
    
    c=c/1e5 #convert to km/s
    if type(disttype) == str:
        disttypemap={'comoving':0,'luminosity':1,'angular':2,'lookback':3,'distmod':4}
        try:
            disttype=disttypemap[disttype]
        except KeyError,e:
            e.message='invalid disttype string'
            raise
    
    flipsign = disttype < 0
    disttype = abs(disttype)
    
    if z is None:
        if normed:
            return 1.0
        if disttype == 2:
            #find maximum value for angular diam dist
            from scipy.optimize import fminbound
            res = upper = 5
            while abs(res-upper) < inttol:
                #-2 flips sign so that we get a minimum instead of a maximum
                res = fminbound(cosmo_z_to_dist,0,upper,(None,-2,inttol,normed,intkwargs),inttol,full_output=1)
                res = -res[1] #this is the actual value -- res[0] is the redshift at which it occurs
            return res
        else:
            #iterate towards large numbers until convergence achieved
            iterz = 1e6
            currval = cosmo_z_to_dist(iterz,None,disttype,inttol,False,intkwargs)
            lastval = currval + 2*inttol
            while(abs(lastval-currval)>inttol):
                lastval = currval
                iterz *= 10
                currval = cosmo_z_to_dist(iterz,None,disttype,inttol,False,intkwargs)
            return currval
        
    z = array(z,copy=False)
    a0 = 1/(z+1)
    omegaK = 1 - omegaM - omegaL - omegaR
    
    if disttype != 3:
        #comoving distance out to scale factor a0: integral(da'/(a'^2 H(a')),a0,1)
        #H^2 a^4=omegaR +omegaM a^1 + omegaE a^4 + omegaK a^2
        def integrand(a,H0,R,M,L,K): #1/(a^2 H)
            return (R + M*a + L*a**4 + K*a**2)**-0.5/H0
    else:
        #lookback time
        def integrand(a,H0,R,M,L,K): #1/(a^2 H)
            return a*(R + M*a + L*a**4 + K*a**2)**-0.5/H0
        
    if isSequenceType(a0):
        integratevec = vectorize(lambda x:integrate(integrand,x,1,args=(H0,omegaR,
                                             omegaM,omegaL,omegaK),**intkwargs))
        res=integratevec(a0)
        intres,interr = res[0],res[1]        
        try:
            if any(interr/intres > inttol):
                raise Exception('Integral fractional error for one of the integrals is beyond tolerance')
        except ZeroDivisionError:
            pass
        
    else:
        res=integrate(integrand,a0,1,args=(H0,omegaR,omegaM,omegaL,omegaK),**intkwargs)
        intres,interr=res[0],res[1]
        
        try:
            if interr/intres > inttol:
                raise Exception('Integral fractional error is '+str(interr/intres)+', beyond tolerance'+str(inttol))
        except ZeroDivisionError:
            pass
    
    if disttype == 3: #lookback integrand
        d = c*intres*3.26163626e-3
        #d = c*intres*3.08568025e19/24/3600/365.25e9
    else: 
        dc = c*intres #comoving distance 
        
        if disttype == 0:
            d = dc
        elif disttype == 1:
            d = dc/a0
        elif disttype == 2:
            if omegaK == 0:
                d = dc*a0
            else:
                angfactor = H0*complex(-omegaK)**0.5
                d = c*(np.sin(angfactor*intres)/angfactor).real*a0
        elif disttype == 4:
            from .phot import distance_modulus
            d = distance_modulus(c*intres/a0*1e6,autocosmo=False)
        else:
            raise KeyError('unknown disttype')
        
    if normed:
        nrm = 1/cosmo_z_to_dist(None if normed is True else normed,None,disttype,inttol,intkwargs)
    else:
        nrm = 1
        
    if flipsign:
        nrm *= -1
        
    if zerr is None:
        return nrm*d
    else:
        if not isscalar(zerr):
            zerr = array(zerr,copy=False) 
        upper=cosmo_z_to_dist(z+zerr,None,disttype,inttol,intkwargs)
        lower=cosmo_z_to_dist(z-zerr,None,disttype,inttol,intkwargs)
        return nrm*d,nrm*(upper-d),nrm*(d-lower)
    
def cosmo_dist_to_z(d,derr=None,disttype=0,inttol=1e-6,normed=False,intkwargs={}):
    """
    Convert a distance to a redshift. See :func:`cosmo_z_to_dist` for meaning of
    parameters. Note that if `d` is None, the maximum distance will be returned.
    """
    from scipy.optimize import brenth
    maxz=10000.0
    
    if derr is not None:
        raise NotImplementedError
    
    if d is None:
        if disttype==2:
            #find maximum value for angular diam dist
            from scipy.optimize import fminbound
            res = upper = 5
            while abs(res-upper) < inttol:
                #-2 flips sign so that we get a minimum instead of a maximum
                res = fminbound(cosmo_z_to_dist,0,upper,(None,-2,inttol,normed,intkwargs),inttol,full_output=1)
                res = res[0] #this is the redshift, -res[1] is the distance value
            return res
        else:
            d = cosmo_z_to_dist(None,None,disttype,inttol,normed,intkwargs)
    
    f=lambda z,dmin:dmin-cosmo_z_to_dist(z,None,disttype,inttol,normed,intkwargs)
    try:
        while f(maxz,d) > 0:
            maxz=maxz**2
    except OverflowError:
        raise ValueError('input distance %g impossible'%float(d))
        
    zval = brenth(f,0,maxz,(d,),xtol=inttol)
    
    return zval
    
    
def cosmo_z_to_H(z,zerr=None):
    """
    Calculates the hubble constant as a function of redshift for the current
    :class:`astropysics.constant.Cosmology` .  
    
    :param z: redshift
    :type z: scalar or array-like
    :param zerr: uncertainty in redshift 
    :type zerr: scalar, array-like, or None
    
    :returns: 
        Hubble constant for the given redshift, or (H,upper_error,lower_error)
        if `zerr` is not None
    """
    from .constants import get_cosmology
    c = get_cosmology()
    if zerr is None:
        return c.H(z)
    else:
        H=c.H(z)
        upper=c.H(z+zerr)
        lower=c.H(z-zerr)
        return H,upper-H,lower-H

def angular_to_physical_size(angsize,zord,usez=True,**kwargs):
    """
    Converts an observed angular size (in arcsec or as an AngularSeperation 
    object) to a physical size.
    
    :param angsize: Angular size in arcsecond.
    :type angsize: float or an :class:`AngularSeperation` object
    :param zord: Redshift or distance
    :type zord: scalar number
    :param usez:
        If True, the input will be interpreted as a redshift, and kwargs
        will be passed into the distance calculation. The result will be in
        pc. Otherwise, `zord` will be interpreted as a distance.
    :type usez: boolean
    
    kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
    
    :returns: a scalar value for the physical size (in pc if redshift is used) 
    """
    if usez:
        d = cosmo_z_to_dist(zord,disttype=2,**kwargs)*1e6 #pc
    else:
        if len(kwargs)>0:
            raise TypeError('if not using redshift, kwargs should not be provided')
        d = zord
    
    if hasattr(angsize,'arcsec'):
        angsize = angsize.arcsec
    sintheta = np.sin(angsize/asecperrad)
    return d*(1/sintheta/sintheta-1)**-0.5
    #return angsize*d/asecperrad

def physical_to_angular_size(physize,zord,usez=True,objout=False,**kwargs):
    """
    Converts a physical size (in pc) to an observed angular size (in arcsec or 
    as an AngularSeperation object if objout is True)
    
    if usez is True, zord is interpreted as a redshift, and cosmo_z_to_dist 
    is used to determine the distance, with kwargs passed into cosmo_z_to_dist 
    otherwise, zord is taken directly as a angular diameter distance (in pc) 
    and kwargs should be absent
    
    :param physize: Physical size in pc
    :type physize: float
    :param zord: Redshift or distance
    :type zord: scalar number
    :param usez:
        If True, the input will be interpreted as a redshift, and kwargs
        will be passed into the distance calculation. The result will be in
        pc. Otherwise, `zord` will be interpreted as a distance.
    :type usez: boolean
    :param objout: 
        If True, return value is an :class:`AngularSeperation` object,
        otherwise, angular size in arcsec.
    :type: bool
    
    kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
    
    :returns: 
        The angular size in acsec, or an :class:`AngularSeperation` object if
        `objout` is True.
        
    """
    if usez:
        d = cosmo_z_to_dist(zord,disttype=2,**kwargs)*1e6 #pc
    else:
        if len(kwargs)>0:
            raise TypeError('if not using redshift, kwargs should not be provided')
        d = zord
        
    r=physize
    res = asecperrad*np.arcsin(r*(d*d+r*r)**-0.5)
    
    if objout:
        return AngularSeperation(res/3600)
    else:
        return res
    
    
#<---------------------DEPRECATED transforms----------------------------------->

#galactic coordate reference positions from IAU 1959 and wikipedia
_galngpJ2000=EquatorialCoordinatesEquinox('12h51m26.282s','+27d07m42.01s')
_galngpB1950=EquatorialCoordinatesEquinox('12h49m0s','27d24m0s')
_gall0J2000=122.932
_gall0B1950=123

def celestial_transforms(ai,bi,transtype=1,epoch='J2000',degin=True,degout=True):
    """
    :deprecated:
    
    Use this to transform between Galactic,Equatorial, and Ecliptic coordinates
    
    transtype can be a number from the table below, or 'ge','eg','gq','qg','gc',
    'cg','cq','qc'
    
    transtype   From           To       |  transtype    From         To
        1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec    
        2     Galactic       RA-DEC     |     5       Ecliptic      Galactic  
        3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
        
    adapted from IDL procedure EULER 
    (http://astro.uni-tuebingen.de/software/idl/astrolib/astro/euler.html)
    """
    #   J2000 coordinate conversions are based on the following constants
    #   (see the Hipparcos explanatory supplement).
    #  eps = 23.4392911111d              Obliquity of the ecliptic
    #  alphaG = 192.85948d               Right Ascension of Galactic North Pole
    #  deltaG = 27.12825d                Declination of Galactic North Pole
    #  lomega = 32.93192d                Galactic longitude of celestial equator  
    #  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
    #  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
    #  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator   
    
    from warnings import warn
    warn('celestial_transforms function is deprecated - use general coordinate transform framework',DeprecationWarning)           

    if epoch == 'B1950':
            psi   = ( 0.57595865315, 4.9261918136,0, 0,0.11129056012, 4.7005372834)     
            stheta =( 0.88781538514,-0.88781538514, 0.39788119938,-0.39788119938, 0.86766174755,-0.86766174755)    
            ctheta =( 0.46019978478, 0.46019978478,0.91743694670, 0.91743694670, 0.49715499774, 0.49715499774)    
            phi =   ( 4.9261918136,  0.57595865315,  0, 0,  4.7005372834, 0.11129056012)

    elif epoch == 'J2000':
            psi   = ( 0.57477043300,4.9368292465,0,0,0.11142137093, 4.71279419371)     
            stheta =( 0.88998808748,-0.88998808748,  0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025)   
            ctheta =( 0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207,  0.49714719172, 0.49714719172)    
            phi  =  ( 4.9368292465,  0.57477043300,  0, 0,        4.71279419371, 0.11142137093)
    else:
            raise ValueError('unknown epoch')
            
    from math import pi
    from numpy import array,sin,cos,arcsin,arctan2
    twopi   =   2.0*pi
    fourpi  =   4.0*pi
    deg_to_rad = 180.0/pi
    
    
    if degin:
        ai,bi=array(ai),array(bi)
    else:
        ai,bi=np.degrees(ai),np.degrees(bi)
    
    if type(transtype) == int:
        i = transtype - 1
    else:
        transd={'ge':1,'eg':0,'gq':1,'qg':0,'gc':5,'cg':4,'cq':3,'qc':2}
        i  = transd[transtype]
    a  = ai/deg_to_rad - phi[i]
    b = bi/deg_to_rad
    sb = sin(b) 
    cb = cos(b)
    cbsa = cb * sin(a)
    b  = -stheta[i] * cbsa + ctheta[i] * sb
    try:
            b[b>1.0]=1.0
    except TypeError: #scalar
            if b > 1:
                    b=array(1.0)
    bo = arcsin(b)*deg_to_rad
    a =  arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
    ao = ( (a+psi[i]+fourpi) % twopi) * deg_to_rad
    
    if not degout:
        ao,bo = np.radians(ao),np.radians(bo)
    
    return ao,bo

_B1950toJ2000xyz=np.matrix([[0.999926,  -0.011179,  -0.004859],
                            [0.011179,   0.999938,  -0.000027],
                            [0.004859,   0.000027,   0.999988]])

def epoch_transform(ra,dec,inepoch='B1950',outepoch='J2000',degrees=True):
    """
    :deprecated:
    """
    from warnings import warn
    warn('epoch_transform function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
    if inepoch != 'B1950' and inepoch != 'J2000':
        raise ValueError('unrecognized epoch '+inepoch)
    if outepoch != 'B1950' and outepoch != 'J2000':
        raise ValueError('unrecognized epoch '+outepoch)
    if degrees:
        ra,dec=np.radians(ra),np.radians(dec)
    else:
        ra,dec=np.array(ra),np.array(dec)
    
    if inepoch == outepoch:
        trans=np.matrix(np.eye(3))
    elif inepoch == 'B1950' and outepoch == 'J2000':
        trans=_B1950toJ2000xyz
    elif inepoch == 'J2000' and outepoch == 'B1950':
        trans=_B1950toJ2000xyz.I
    else:
        raise ('unrecognized epochs')
    
    x=np.cos(ra)*np.cos(dec)
    y=np.sin(ra)*np.cos(dec)
    z=np.sin(dec)
    
    v=np.matrix((x,y,z))
    xp,yp,zp=trans*v
    
    rap=np.arctan2(yp,xp)
    decp=np.arcsin(zp)
    
    return rap,decp

def galactic_to_equatorial(l,b,epoch='J2000',strout=None):
    """
    :deprecated:
    
    convinience function for celestial_transforms
    if strout is None, will automatically decide based on inputs
    """
    from warnings import warn
    warn('galactic_to_equatorial function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
    from operator import isSequenceType
    
    if type(l) == str:
        l=AngularCoordinate(l).degrees
        if strout is None:
            strout=True
    if type(b) == str:
        b=AngularCoordinate(b).degrees
        if strout is None:
            strout=True
    ra,dec = celestial_transforms(l,b,transtype='ge',epoch=epoch)
    if strout:
        if not isSequenceType(ra):
            ra=[ra]
        if not isSequenceType(dec):
            dec=[dec]
        rao,deco=[],[]
        for rai in ra:
            rao.append(AngularCoordinate(rai).getHmsStr())
        for deci in dec:
            deco.append(AngularCoordinate(deci).getDmsStr())
        return rao,deco
    else:
        return ra,dec
    
    
def equatorial_to_galactic(ra,dec,epoch='J2000',strout=None):
    """
    :deprecated:
    
    convinience function for celestial_transforms
    if strout is None, will automatically decide based on inputs
    """
    from warnings import warn
    warn('equatorial_to_galactic function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
    from operator import isSequenceType
    
    if type(ra) == str:
        ra=AngularCoordinate(ra).degrees
        if strout is None:
            strout=True
    if type(dec) == str:
        dec=AngularCoordinate(dec).degrees
        if strout is None:
            strout=True
    l,b = celestial_transforms(ra,dec,transtype='eg',epoch=epoch)
    if strout:
        if not isSequenceType(l):
            l=[l]
        if not isSequenceType(b):
            b=[b]
        lo,bo=[],[]
        for li in l:
            lo.append(AngularCoordinate(li).getDmsStr())
        for bi in b:
            bo.append(AngularCoordinate(bi).getDmsStr())
        return lo,bo
    else:
        return l,b
