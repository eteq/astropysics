#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains objects and functions for specifying locations as well 
as calculating distances and similar tools.

Some of the calculations involved make use of the currently selected cosmology 
(see astropysics.constants) and hence may not function as expected if a 
particularly strange cosmology is in use.

*see also*
Kapteyn libraries: http://www.astro.rug.nl/software/kapteyn/index.html
Pyephem: http://rhodesmill.org/pyephem/
"""

#TODO: WCSlib or similar support - Kapteyn?


from __future__ import division,with_statement
from .constants import pi
import numpy as np

#<----------------coordinate classes and related functions------------------>


class AngularCoordinate(object):
    """
    An angular coordinate that can be initialized in various formats.
    
    Arithmetic operators can be applied to the coordinate, and will be applied 
    directly to the numerical value in radians.  For + and -, two angular 
    coordinates may be used, although for -, an AngularSeperation object will
    be returned
    """
    import re as _re
    __slots__=('_decval','_range')
    
    def __setdegminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set as a length-3 iterator')
        self.degrees=abs(dms[0])+abs(dms[1])/60.+abs(dms[2])/3600.
        if dms[0]<0:
            self._decval*=-1
    def __getdegminsec(self):
        fulldeg=abs(self.degrees)
        deg=int(fulldeg)
        fracpart=fulldeg-deg
        min=int(fracpart*60.)
        sec=fracpart*3600.-min*60.
        return -deg if self.degrees < 0 else deg,min,sec
    
    def __sethrsminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set as a length-3 iterator')
        self.degrees=15*(dms[0]+dms[1]/60.+dms[2]/3600.)
    def __gethrsminsec(self):
        factorized=self.degrees/15.
        hrs=int(factorized)
        mspart=factorized-hrs
        min=int(mspart*60.)
        sec=mspart*3600.-min*60.
        return hrs,min,sec
    
    def __setdegdec(self,deg):
        rads = deg*pi/180.
        if self.range is not None:
            self.__checkRange(rads)
        self._decval = rads
    def __getdegdec(self):
        return self._decval*180/pi
    
    def __setrad(self,rads):
        if self.range is not None:
            self.__checkRange(rads)
        self._decval = rads
    def __getrad(self):
        return self._decval
    
    def __sethrdec(self,hr):
        rads = hr*pi/12
        if self.range is not None:
            self.__checkRange(rads)
        self._decval = rads
    def __gethrdec(self):
        return self._decval*12/pi
    
    def __checkRange(self,rads):
        if self._range is not None:
            low,up = self._range
            if not low <= rads <= up:
                raise ValueError('Attempted to set angular coordinate outside range')
    def __setrange(self,newrng):
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
            self.__checkRange(self._decval)
        except ValueError:
            self._range = oldrange
            raise ValueError('Attempted to set range when value is out of range')
    def __getrange(self):
        if self._range is None:
            return None
        else:
            from math import degrees
            return degrees(self._range[0]),degrees(self._range[1])
    
    def __str__(self):
        return self.getDmsStr()
    
    degminsec=property(fget=__getdegminsec,fset=__setdegminsec)
    dms=degminsec
    hrsminsec=property(fget=__gethrsminsec,fset=__sethrsminsec)
    hms=hrsminsec
    degrees=property(fget=__getdegdec,fset=__setdegdec)
    d=degrees
    radians=property(fget=__getrad,fset=__setrad)
    r=radians
    hours=property(fget=__gethrdec,fset=__sethrdec)
    h=hours
    range=property(fget=__getrange,fset=__setrange)
    
    ############################################################################
    #                                                                          #
    #                                                                          #
    #  and then a monkey showed up                                             #
    #  and gave you a kiss                                                     #
    #  and you were confused                                                   #
    #                                                                          #
    #                                             comment poem by Rie O:-)     #
    ############################################################################
   
    __purehre=_re.compile(r'.*?(\d+(?:\.?\d+))(?:h|hr).*')
    __hmre=_re.compile('.*?(\d{1,2})(?:h|hr)\s*(\d{1,2})(?:m|min).*')
    __hmsre=_re.compile(r'.*?(\d{1,2})(?:h|hr)\s*(\d{1,2})(?:m|min)\s*(\d+(?:\.?\d*))(?:s|sec).*')
    __puredre=_re.compile(r'.*?([+-]?\s*\d+(?:\.?\d+))(?:d|deg).*')
    __dmre=_re.compile(r'.*?([+-])?(\d{1,2})(?:d|deg)\s*(\d{1,2})(?:m|min).*')
    __dmsre=_re.compile(r'.*?([+-])?(\d{1,2})(?:d|deg)\s*(\d{1,2})(?:m|min)\s*(\d+(?:\.?\d*))(?:s|sec).*')
    __sexre=_re.compile(r'.*?(\+|\-)?(\d{1,3})[: ](\d{1,2})[: ](\d+(?:.\d+)?).*')
    __radsre=_re.compile(r'.*?(\d+(?:\.?\d+))(?:r|rad).*')
    def __init__(self,inpt=None,sghms=None,range=None,radians=False):
        """
        If an undecorated 3-element iterator, `inpt` is taken to be deg,min,sec, 
        othewise, input is cast to a float and treated as decimal degrees
        
        if sexigesimal, input will be interpreted as h:m:s if `sghms`
        is True, or d:m:s if `sghms` is False.  If None, a +/- will indicate
        d:m:s and nothing indicates r:m:s
        
        `range` sets the valid range of coordinates either any value (if None)
        or a 2-sequence (lowerdegrees,upperdegrees)
        
        if `radians` is True and the input is a float or ambiguous string, the
        value will be taken to be in radians, otherwise degrees
        """
        self._range = None
        
        if isinstance(inpt,AngularCoordinate):
            self._decval=inpt._decval
        elif isinstance(inpt,basestring):
            sexig=self.__sexre.match(inpt)
            hm=self.__purehre.match(inpt)
            hmm=self.__hmre.match(inpt)
            hmsm=self.__hmsre.match(inpt)
            dm=self.__puredre.match(inpt)
            dmm=self.__dmre.match(inpt)
            dmsm=self.__dmsre.match(inpt)
            radsm=self.__radsre.match(inpt)
            if sexig:
                t=sexig.group(2,3,4)
                if sghms is None:
                    if sexig.group(1) is None:
                        self.hrsminsec=int(t[0]),int(t[1]),float(t[2])
                    else:
                        sgn = 1 if sexig.group(1) == '+' else -1
                        self.degminsec=int(t[0]),int(t[1]),float(t[2])
                        self._decval *= sgn
                else:
                    sgn = -1 if sexig.group(1) == '-' else 1
                    if sghms:
                        self.hrsminsec=int(t[0]),int(t[1]),float(t[2])
                    else:
                        self.degminsec=int(t[0]),int(t[1]),float(t[2]) 
                    self._decval *= sgn
            elif hmsm:
                t=hmsm.group(1,2,3)
                self.hrsminsec=int(t[0]),int(t[1]),float(t[2])
            elif hmm:
                t=hmm.group(1,2)
                self.hrsminsec=int(t[0]),int(t[1]),0
            elif hm:
                self.hours=float(hm.group(1))
            elif dmsm:
                sgn = -1 if dmsm.group(1) =='-' else 1
                t=dmsm.group(2,3,4)
                self.degminsec=int(t[0]),int(t[1]),float(t[2])
                self._decval *= sgn
            elif radsm:
                self.radians=float(hm.group(1))
            elif dmm:
                t=dmm.group(1,2)
                self.degminsec=int(t[0]),int(t[1]),0
            elif dm:
                self.degrees=float(dm.group(1))
            else:
                try:
                    if radians:
                        self.radians = float(inpt)
                    else:
                        self.degrees = float(inpt)
                except:
                    raise ValueError('Unrecognized string format '+inpt)
            
        elif hasattr(inpt,'__iter__') and len(inpt)==3:
            self.degminsec=inpt
        elif inpt is None:
            self._decval=0
        elif radians:
            self._decval=float(inpt)
        else:
            self._decval=float(inpt)*pi/180.
        
        self.range = range
            
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
        
    #unichr(176) for deg symbol
    def getDmsStr(self,secform='%05.2f',sep=('d',"'",'"'), sign=True, canonical=False):
        """
        gets the string representation of this AngularCoordinate as degrees,
        minutes, and seconds
        
        `secform` should be a string to use as a formatter for the seconds
        
        `sep` is the seperator between components - defaults to d,
        ' and " symbols, can be a single string or a 3-tuple of strings
        
        `sign` forces sign to be present before degree component
        
        `canonical` forces [+/-]dd:mm:ss.ss , overriding other arguments
        """
        d,m,s = self.degminsec
        
        if canonical:
            sgn = '-' if self._decval < 0 else '+'
            return '%s%02.i:%02.i:%05.2f'%(sgn,d,m,s)
        
        d,m=str(d),str(m)
        
        s = secform%s
        
        if isinstance(sep,basestring):
            if sep == 'dms':
                sep = ('d','m','s')
            sep = (sep,sep)
        
        tojoin = []
        
        if sign:
            d='+'+d if d >= 0 else d
        
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
    
    def getSizeAtDistance(self,distance):
        """
        computes the size of an angular seperation given by this coordinate
        """
        from math import cos
        return distance*(2-2*cos(self.radians))**0.5
    
def angular_string_to_dec(instr,hms=True):
    """
    convinience function to convert a sexigesimal angular position to a decimal
    value.
    
    if hms is True, the coordinate will be assumed to be h:m:s, otherwise d:m:s
    """
    ac = AngularCoordinate(instr)
    return ac.degrees

class AngularSeperation(AngularCoordinate):
    __slots__ = ('start',)
    def __init__(self,*args,**kwargs):
        """
        inputs can be either AngularSeperation(sep), or
        AngularSeperation(start,end).  kwargs can be any of the kwargs to
        the initializer of AngularCoordinate.
        """
        if len(args) == 1:
            a = args[0]
            if a.__class__ == self.__class__:
                self._decval = args[0]._decval
                self._range = args[0]._range
                self.start = args[0].start
                return
            
            sep = a._decval if hasattr(a,'_decval') else a
            start = None
        elif len(args) == 2:
            a0,a1 = args
            a0 = a0._decval if hasattr(a0,'_decval') else a0
            a1 = a1._decval if hasattr(a1,'_decval') else a1
            sep = a1 - a0
            start = a0
        else:
            raise ValueError('inproper number of inputs to AngularSeperation')
        
        super(AngularSeperation,self).__init__(sep,**kwargs)
        self.start = start
        
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
        computes the physical projected seperation assuming a given distance.
        This implicitly assumes small-angle approximation.
        
        if `usez` is True, the input will be interpreted as a redshift, and
        kwargs will be passed into the distance calculation.
        """
        return angular_to_physical_size(self.arcsec,zord,usez=True,**kwargs)

class _LatLongMeta(type):
    def __init__(cls,name,bases,dct):
        super(_LatLongMeta,cls).__init__(name,bases,dct)
        if cls._latlongnames_[0] is not None:
            setattr(cls,cls._latlongnames_[0],cls.lat)
            setattr(cls,cls._latlongnames_[0]+'err',cls.laterr)
        if cls._latlongnames_[1] is not None:
            setattr(cls,cls._latlongnames_[1],cls.long)
            setattr(cls,cls._latlongnames_[1]+'err',cls.longerr)
            
#    def __call__(cls,*args,**kwargs):
#        obj = super(_LatLongMeta,_LatLongMeta).__call__(cls,**objkwargs) #object __init__ is called here
#        return obj

class LatLongPosition(object):
    """
    This object represents an angular location on a unit sphere as represented
    in spherical coordinates with a latitude and longitude.  Subclasses specify 
    details such as transformations or epochs.
    """
    __slots__ = ('_lat','_long','_laterr','_longerr')
    __metaclass__ = _LatLongMeta
    _latlongnames_ = (None,None)
    _longrange_ = None
    
    def __init__(self,lat=0,long=0,laterr=None,longerr=None):
        """
        `lat` and `long` are latitude and longitude for the position, respectively,
        and may be any valid input to AngularCoordinate (if a number, default
        is in degrees)        
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
                raise ValueError("can't provide a LatLongPosition as a constructor and set other values simultaneously")
        else:
            self.lat = lat
            self.long = long
            self.laterr = laterr
            self.longerr = longerr
        
    def __getstate__(self):
        return dict([(k,getattr(k)) for k in LatLongPosition.__slots__])
    
    def __setstate__(self,d):
        for k in LatLongPosition.__slots__:
            setattr(k,d[k])
        
    def _getLat(self):
        return self._lat
    def _setLat(self,val):
        if isinstance(val,AngularCoordinate):
            self._lat.radians = val.radians
        else:
            self._lat.radians = AngularCoordinate(val).radians
    lat = property(_getLat,_setLat,doc=None)
    
    def _getLong(self):
        return self._long
    def _setLong(self,val):
        if isinstance(val,AngularCoordinate):
            self._long.radians = val.radians
        else:
            self._long.radians = AngularCoordinate(val).radians
    long = property(_getLong,_setLong,doc=None)
    
    def _getLaterr(self):
        return self._laterr
    def _setLaterr(self,val):
        if val is None:
            self._laterr = None
        elif isinstance(val,AngularSeperation):
            self._laterr = val
        else:
            self._laterr = AngularSeperation(val)
    laterr = property(_getLaterr,_setLaterr,doc=None)
    
    def _getLongerr(self):
        return self._longerr
    def _setLongerr(self,val):
        if val is None:
            self._longerr = None
        elif isinstance(val,AngularSeperation):
            self._longerr = val
        else:
            self._longerr = AngularSeperation(val)
    longerr = property(_getLongerr,_setLongerr,doc=None)
    
    def __str__(self):
        return '{0}:({1[0]}={2},{1[1]}={3})'.format(self.__class__.__name__,self._latlongnames_,self._lat.d,self._long.d)
    
    def __eq__(self,other):
        if hasattr(other,'_lat') and hasattr(other,'_long'):
            return self._lat==other._lat and self._long==other._long
        else:
            return False
        
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __sub__(self,other):
        if isinstance(other,self.__class__):
            from math import cos
            dcorrlong = self._long.radians * cos(self._lat.radians) \
                 - other._long.radians * cos(other._lat.radians)
            dlat = self._lat.radians-other._lat.radians
            sep = AngularSeperation(degrees((dlat*dlat+dcorrlong*dcorrlong)**0.5))
            sep.start = other
            return sep
        else:
            raise "unsupported operand type(s) for -: '%s' and '%s'"%(self.__class__,other.__class__)
        
    def transform(self,matrix,apply=True,unitarycheck=False):
        """
        Applies the supplied transformation matrix (in cartesian coordinates) to
        these coordinates.  The transform matrix is assumed to be unitary, 
        although this is not checked unless `unitarycheck` is True
        
        if `apply` is True, the transform will be applied to this coordinate as
        well as being returned
        
        *returns*
        lat,long after the supplied transform is applied
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
        
        #spherical w/ r=1 > cartesian
        x = cos(lat)*cos(long)
        y = cos(lat)*sin(long) 
        z = sin(t)
        
        #do transform
        xp,yp,zp = m*np.matrix((x,y,z))
        
        #cartesian > spherical
        latp = atan2(zp,sqrt(xp*xp+yp*yp))
        longp = atan2(yp,xp)
        
        if apply:
            self.lat = latp
            self.long = longp
        
        return latp,longp

class HorizontalPosition(LatLongPosition):
    """
    This object represents an angular location on the unit sphere, with the 
    north pole of the coordinate position fixed to the local zenith
    
    To convert from other fixed position types to horizontal positions, see 
    :class:`astropysics.obstools.Site`, as site information is required for
    these corrections
    """  
    __slots__ = tuple()
    _latlongnames_ = ('alt','az')
    _longrange_ = (0,360)
    

class EquatorialPosition(LatLongPosition):
    """
    This object represents an angular location on the unit sphere, specified in
    right ascension and declination.  Thus, the fundamental plane is given by 
    the projection of the equator in the sky.
    """
    
    __slots__ = ('_epoch','_refsys')
    _latlongnames_ = ('dec','ra')
    _longrange_ = (0,360)
    
    def __getstate__(self):
        d = super(EquatorialPosition,self).__getstate__
        d['_epoch'] = self._epoch
        return d
    
    def __setstate__(self,d):
        super(EquatorialPosition,self).__setstate__
        self._epoch = d['_epoch']
    
    def __init__(self,*args,**kwargs):
        """
        input may be in the following forms:
        
        * EquatorialPosition()
        * EquatorialPosition(EquatorialPosition)
        * EquatorialPosition('rastr decstr')
        * EquatorialPosition(ra,dec)
        * EquatorialPosition(ra,dec,raerr,decerr)
        * EquatorialPosition(ra,dec,raerr,decerr,epoch)
        * EquatorialPosition(ra,dec,raerr,decerr,epoch,refsys) 
        
        """
        posargs = {}
        if len(args) == 0:
            pass
        if len(args) == 1:
            if isinstance(args[0],EquatorialPosition):
                super(EquatorialPosition,self).__init__(args[0])
                self.epoch = args[0].epoch
                self.refsys = args[0].refsys
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
        
        for k,v in posargs.iteritems():
            if k in kwargs:
                raise ValueError('got multiple values for argument '+k)
            kwargs[k] = v
        
        kwargs.setdefault('ra',0)
        kwargs.setdefault('dec',0)
        kwargs.setdefault('raerr',None)
        kwargs.setdefault('decerr',None)
        kwargs.setdefault('epoch','J2000')
        kwargs.setdefault('refsys','ICRS')
        
        super(EquatorialPosition,self).__init__(kwargs['dec'],kwargs['ra'],kwargs['decerr'],kwargs['raerr'])
        self.epoch = kwargs['epoch']
        self.refsys = kwargs['refsys']
            
    def _getEpoch(self):
        return self._epoch
    def _setEpoch(self,val):
        if hasattr(self,'_epoch') and self._epoch != val:
            from warnings import warn
            warn('epoch transforms not ready yet')
        if not isinstance(val,basestring):
            val = 'J'+str(float(val))
        self._epoch = val
    epoch = property(_getEpoch,_setEpoch,doc=None)
    
    def _getRefsys(self):
        return self._refsys
    def _setRefsys(self,val):
        self._refsys = val
        if hasattr(self,'_refsys') and self._refsys != val:
            from warnings import warn
            warn('refsys transforms not ready yet')
    refsys = property(_getRefsys,_setRefsys,doc=None)
    
    def toGal(self):
        """
        converts this position to Galactic coordinates
        """
        newpos = EquatorialPosition(self)
        newpos.epoch = 'J2000'
        
        latang = GalacticPosition._ngp_J2000.lat.d
        longang = GalacticPosition._ngp_J2000.long.d
        long0 = GalacticPosition._long0_J2000.d
        
        mrot = rotZ(180-long0)*rotY(90-latang)*rotZ(longang)
        newpos.transform(mrot)
        return GalacticPosition(newpos)
    
class EclipticPosition(LatLongPosition):
    """
    Ecliptic Coordinates (beta, lambda) such that the fundamental plane passes
    through the ecliptic at the given epoch.
    """
    
    __slots__ = ('_eclipticepoch','_refsys')
    _latlongnames_ = ('beta','lamb')
    _longrange_ = (0,360)
    
    def __init__(self,beta=0,lamb=0,betaerr=None,lamberr=None,eclipticepoch='J2000',refsys='ICRS'):
        raise NotImplementedError
    
    def __getstate__(self):
        d = super(EquatorialPosition,self).__getstate__
        d['_eclipticepoch'] = self._eclipticepoch
        return d
    
    def __setstate__(self,d):
        super(EquatorialPosition,self).__setstate__
        self._eclipticepoch = d['_eclipticepoch']
        
    def _getEclipticepoch(self):
        return self._epoch
    def _setEclipticepoch(self,val):
        if hasattr(self,'_eclipticepoch') and self._eclipticepoch != val:
            from warnings import warn
            warn('warning: epoch transforms not ready yet')
        if not isinstance(val,basestring):
            val = 'J'+str(float(val))
        self._eclipticepoch = val
    eclipticepoch = property(_getEclipticepoch,_setEclipticepoch,doc=None)
    
    def _getRefsys(self):
        return self._refsys
    def _setRefsys(self,val):
        self._refsys = val
        if hasattr(self,'_refsys') and self._refsys != val:
            from warnings import warn
            warn('warning: refsys transforms not ready yet')
    refsys = property(_getRefsys,_setRefsys,doc=None)
    
    def __init__(self,*args,**kwargs):
        raise NotImplementedError
    
class GalacticPosition(LatLongPosition):
    __slots__ = tuple()
    _latlongnames_ = ('b','l')
    _longrange_ = (0,360)
    
    _ngp_J2000 = EquatorialPosition(192.859508, 27.128336,epoch='J2000')
    _long0_J2000 = AngularCoordinate(122.932)
    
    def toSGal(self):
        """
        converts this position to Supergalactic coordinates
        """
        latang = SupergalacticPosition._nsgp_gal.lat.d
        longang = SupergalacticPosition._nsgp_gal.long.d
        mrot = rotation_matrix(90-latang,'x')*rotation_matrix(90+longang,'z')
        sgp = SupergalacticPosition(self)
        sgp.transform(mrot)
        return sgp
    
    
class SupergalacticPosition(LatLongPosition):   
    __slots__ = tuple()
    _latlongnames_ = ('sgb','sgl')
    _longrange_ = (0,360)
    
    _nsgp_gal = GalacticPosition(6.32,47.37)
    _sglong0_gal = 137.37
    _nsgp_J2000 = EquatorialPosition(283.75420420,15.70894043,epoch='J2000')
    _sglong0_J2000 = 42.30997710
    

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
                if isinstance(o,EquatorialPosition):
                    coords.append((o.ra.d,o.dec.d))
                else:
                    coords.append((o.lat.d,o.long.d))
            else:
                coords.append([getattr(o,c).d for c in coordnames])
    else:
        for o in posobjs:
            if coordnames is None:
                if isinstance(o,EquatorialPosition):
                    coords.append((o.ra.r,o.dec.r))
                else:
                    coords.append((o.lat.r,o.long.r))
            else:
                coords.append([getattr(o,c).r for c in coordnames])

    return np.array(coords).T

def rotation_matrix(angle,axis='z',degrees=True):
    """
    generates a 3x3 rotation matrix in cartesian coordinates for rotation about
    the requested axis.
    
    `axis` can either be a string 'x','y', or 'z', or a 3-sequence specifying
    an axis to rotate about.
    
    if `degrees` is True the input angle is in degrees, otherwise it is radians
    """
    from math import sin,cos,radians,sqrt
    if degrees:
        angle = radians(angle)
        
    
    
    if axis == 'z':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c,-s, 0),
                          ( s, c, 0),
                          ( 0, 0, 1)))
    elif axis == 'y':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, 0,-s),
                          ( 0, 1, 0),
                          ( s, 0, c)))
    elif axis == 'x':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( 1, 0, 0),
                          ( 0, c,-s),
                          ( 0, s, c)))
    else:
        x,y,z = axis
        w = cos(angle/2)
        
        #normalize
        if w == 1:
            x=y=z=0
        else:
            l = sqrt((x*x + y*y + z*z)/(1-w*w))
            x /= l
            y /= l
            z /= l
        
        wsq = w*w
        xsq = x*x
        ysq = y*y
        zsq = z*z
        return np.matrix((( wsq+xsq-ysq-zsq, 2*x*y-2*w*z, 2*x*z+2*w*y),
                          ( 2*x*y+2*w*z, wsq-xsq+ysq-zsq,2*y*z-2*w*x),
                          ( 2*x*z-2*w*y, 2*y*z+2*w*x, wsq-xsq-ysq+zsq)))
def angle_axis(matrix,degrees=True):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix
    
    *returns*
    angle,axis where angle is in degrees if `degrees` is True, otherwise in
    radians
    """
    from math import sin,cos,acos,degrees,sqrt
    
    m = np.asmatrix(matrix)
    if m.shape != (3,3):
        raise ValueError('matrix is not 3x3')
    
    
    
    angle = acos((m[0,0] + m[1,1] + m[2,2] - 1)/2)
    denom = sqrt(2*((m[2,1]-m[1,2])+(m[0,2]-m[2,0])+(m[1,0]-m[0,1])))
    axis = np.array((m[2,1]-m[1,2],m[0,2]-m[2,0],m[1,0]-m[0,1]))/denom
    axis /= sqrt(np.sum(axis**2)) 
    
    if degrees:
        return degrees(angle),axis
    else:
        return angle,axis
    
    
    
#<--------------------canonical coordinate transforms-------------------------->
def cartesian_to_polar(x,y,degrees=False):
    """
    Converts two arrays in 2D rectangular Cartesian coordinates to
    polar coordinates.
    
    if degrees is True, the output theta angle will be in degrees, 
    otherwise radians
    
    returns (r,theta) where theta is measured from the +x axis increasing
    towards the +y axis
    """
    r = (x*x+y*y)**0.5
    t = np.arctan2(y,x)
    if degrees:
        t = np.degrees(t)
    
    return r,t

def polar_to_cartesian(r,t,degrees=False):
    """
    Converts two arrays in polar coordinates to
    2D rectangular Cartesian coordinates.
    
    if degrees is True, the input angle t will be interpreted as given
    in degrees, otherwise radians
    
    theta is measured from the +x axis increasing
    towards the +y axis
    
    returns (x,y)
    """
    if degrees:
        t=np.radians(t)
        
    return r*np.cos(t),r*np.sin(t)

def cartesian_to_spherical(x,y,z,degrees=False):
    """
    Converts three arrays in 3D rectangular cartesian coordinates to
    spherical polar coordinates.
    
    returns r,theta,phi in PHYSICIST convention - (1,0,pi/2) is x-axis
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
    Converts three arrays in 3D spherical polar coordinates to
    rectangular cartesian coordinates.
    
    if degrees is true, converts from degrees to radians for theta and phi
    returns x,y,z in PHYSICIST convention - (1,0,pi/2) is x-axis
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
    converts from colatitude/inclination (i.e. "theta" in physicist convention) 
    to latitude (i.e. 0 at the equator).
    """
    if degrees:
        return 90 - theta
    else:
        return pi/2 - theta

def cartesian_to_cylindrical(x,y,z,degrees=False):
    """
    Converts three arrays in 3D rectangular Cartesian coordinates
    to cylindrical polar coordinates. 
    
    if degrees is true, outputtheta is in degrees, otherwise radians
    returns s,theta,z (theta increasing from +x to +y, 0 at x-axis)
    """
    s,t = cartesian_to_polar(x,y)
    return s,t,z
    
def cylindrical_to_cartesian(s,t,z,degrees=False):
    """
    Converts three arrays in cylindrical polar coordinates to
    3D rectangular Cartesian coordinates. 
    
    if degrees is true, converts from degrees to radians for theta
    returns x,y,z (theta increasing from +x to +y, 0 at x-axis)
    """
    x,y = poloar_to_cartesian(s,t,degrees)
    return x,y,z
    
def proj_sep(rx,ty,pz,offset,spherical=False):
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




def spherical_distance(ra1,dec1,ra2,dec2,degrees=True):
    ra1,dec1 = np.array(ra1,copy=False),np.array(dec1,copy=False)
    ra2,dec2 = np.array(ra2,copy=False),np.array(dec2,copy=False)
    
    x1,y1,z1 = spherical_to_cartesian(1.0,90-dec1 if degrees else dec1,ra1,degrees)
    x2,y2,z2 = spherical_to_cartesian(1.0,90-dec2 if degrees else dec2,ra2,degrees)
    
    dp = x1*x2+y1*y2+z1*z2
    
    sep = np.arccos(dp)
    if degrees:
        sep = np.degrees(sep)
    return sep

def seperation3d(d1,d2,pos1,pos2):
    """
    Compute the full 3d seperation between two objects at distances `d1` and 
    `d2` and angular positions `pos1` and `pos2` (LatLongPositions, or an
    argument to initialize a new EquatorialPosition)
    """
    from numpy import sin,cos
    if not isinstance(pos1,LatLongPosition):
        pos1=EquatorialPosition(pos1)
    if not isinstance(pos2,LatLongPosition):
        pos2=EquatorialPosition(pos2)
    
    theta1,phi1=(pi/2-pos1.lat.radians,pos1.long.radians)
    theta2,phi2=(pi/2-pos2.lat.radians,pos2.long.radians)
    
    dx=d2*sin(theta2)*cos(phi2)-d1*sin(theta1)*cos(phi1)
    dy=d2*sin(theta2)*sin(phi2)-d1*sin(theta1)*sin(phi1)
    dz=d2*cos(theta2)-d1*cos(theta1)
    
    return (dx*dx+dy*dy+dz*dz)**0.5

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
    Match one set of coordinates to another within a tolerance eps
    e.g. ra1,dec1,ra2,dec2
    
    returns (mask of matches for array 1, mask of matches for array 2)
    
    if multi is 'warn', a arning will be issued if more than 1 match is found
    if multi is 'print', a statement will be printed  if more than 1 match is found
    if multi is 'full', the 2D array with matches as booleans along the axes will be returned
    if multi is 'count', a count of matches will be returned instead of a mask
    if multi is 'index', a list of match indecies will be returned instead of a mask
    if it otherwise evaluates to True, an exception will be raised
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
    

#<-----------------Cosmological distance and conversions ---------->
def cosmo_z_to_dist(z,zerr=None,disttype=0,inttol=1e-6,normed=False,intkwargs={}):
    """
    calculates the cosmolgical distance to some object given a redshift
    note that this uses H0,omegaM,omegaL, and omegaR from the current cosmology
    
    intkwargs are kwargs for the integration func (usually scipy.integrate.quad)
    
    if z is None, returns the limit of z->inf to the requested tolerance (or at 
    maximum if angular) - this may be infinite
    
    disttype can be 'comoving'(0),'luminosity'(1),'angular'(2),'lookback'(3),
    or 
    'distmod'(4)
    for the angular diameter distance, a flat universe is assumed for simple 
    dang=a*chi
    
    output in Mpc or Gyr or if normed is True, normalized to value for z=None, 
    or if normed is a number, normalized to z=normed
    """
    from operator import isSequenceType
    from scipy.integrate import quad as integrate
    from numpy import array,vectorize,abs
    
    from .constants import H0,omegaM,omegaL,omegaR,c
    
    c=c/1e5 #convert to km/s
    if type(disttype) == str:
        disttypemap={'comoving':0,'luminosity':1,'angular':2,'lookback':3,'distmod':4}
        try:
            disttype=disttypemap[disttype]
        except KeyError,e:
            e.message='invalid disttype string'
            raise
    if omegaL+omegaM+omegaR != 1:
        from warnings import warn
        warn('cosmological parameters do not sum to 1 - algorithm may not be correct for non-flat universe')
    
    if z is None:
        if normed:
            return 1.0
        if disttype == 2:
            #need to find maximum value for angular diam dist
            from scipy.optimize import brent
            #-1 on the H0 will make the function negative for purposes of minimizing
            raise NotImplemetedError('need to fix -H0')
            return -1*brent(cosmo_z_to_dist,(None,disttype,-1*H0,omegaM,omegaL,
                            omegaR,inttol,False,intkwargs),tol=inttol,full_output=1)[1]
        else:
            #iterate towards large numbers until convergence achieved
            iterz=1e6
            currval=cosmo_z_to_dist(iterz,None,disttype,inttol,False,intkwargs)
            lastval=currval+2*inttol
            while(abs(lastval-currval)>inttol):
                lastval=currval
                iterz*=10
                currval=cosmo_z_to_dist(iterz,None,disttype,inttol,False,intkwargs)
            return currval
        
    z=array(z)
    a0=1/(z+1)
    omegaK=1-omegaM-omegaL-omegaR
    
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
        integratevec=vectorize(lambda x:integrate(integrand,x,1,args=(H0,omegaR,
                                             omegaM,omegaL,omegaK),**intkwargs))
        res=integratevec(a0)
        intres,interr=res[0],res[1]        
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
    
    if disttype == 0:
        d=c*intres
    elif disttype == 1:
        d=c*intres/a0
    elif disttype == 2:
        d=a0*c*intres
    elif disttype == 3:
        d=c*intres*3.26163626e-3
        #d=chi*3.08568025e19/24/3600/365.25e9
    elif disttype == 4:
        from .phot import distance_modulus
        d=distance_modulus(c*intres/a0*1e6,autocosmo=False)
        
    else:
        raise KeyError('unknown disttype')
    
    if normed:
            nrm=1/cosmo_z_to_dist(None if normed is True else normed,None,disttype,inttol,intkwargs)
    else:
        nrm=1
        
    if zerr is None:
        return nrm*d
    else: 
        upper=cosmo_z_to_dist(z+zerr,None,disttype,inttol,intkwargs)
        lower=cosmo_z_to_dist(z-zerr,None,disttype,inttol,intkwargs)
        return nrm*d,nrm*(upper-d),nrm*(lower-d)
    
def cosmo_dist_to_z(d,derr=None,disttype=0,inttol=1e-6,normed=False,intkwargs={}):
    from scipy.optimize import brenth
    maxz=10000.0
    
    if derr is not None:
        raise NotImplementedError
    
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
    calculate the hubble constant as a function of redshift for the selected
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
    converts an observed angular size (in arcsec or as an AngularSeperation 
    object) to a physical size (in pc)
    
    if usez is True, zord is interpreted as a redshift, and cosmo_z_to_dist 
    is used to determine the distance, with kwargs passed into cosmo_z_to_dist 
    otherwise, zord is taken directly as a angular diameter distance (in pc) 
    and kwargs should be absent
    """
    if usez:
        d = cosmo_z_to_dist(zord,disttype=2,**kwargs)*1e6 #pc
    else:
        if len(kwargs)>0:
            raise TypeError('if not using redshift, kwargs should not be provided')
        d = zord
    
    if hasattr(angsize,'arcsec'):
        angsize = angsize.arcsec
    sintheta = np.sin(angsize/206265)
    return d*(1/sintheta/sintheta-1)**-0.5
    #return angsize*d/206265

def physical_to_angular_size(physize,zord,usez=True,objout=False,**kwargs):
    """
    converts a physical size (in pc) to an observed angular size (in arcsec or 
    as an AngularSeperation object if objout is True)
    
    if usez is True, zord is interpreted as a redshift, and cosmo_z_to_dist 
    is used to determine the distance, with kwargs passed into cosmo_z_to_dist 
    otherwise, zord is taken directly as a angular diameter distance (in pc) 
    and kwargs should be absent
    """
    if usez:
        d = cosmo_z_to_dist(zord,disttype=2,**kwargs)*1e6 #pc
    else:
        if len(kwargs)>0:
            raise TypeError('if not using redshift, kwargs should not be provided')
        d = zord
        
    r=physize
    res = 206265*np.arcsin(r*(d*d+r*r)**-0.5)
    if objout:
        return AngularSeperation(res/3600)
    else:
        return res
    
    
#<---------------------DEPRECATED transforms----------------------------------->

#galactic coordate reference positions from IAU 1959 and wikipedia
_galngpJ2000=EquatorialPosition('12h51m26.282s','+27d07m42.01s')
_galngpB1950=EquatorialPosition('12h49m0s','27d24m0s')
_gall0J2000=122.932
_gall0B1950=123

def celestial_transforms(ai,bi,transtype=1,epoch='J2000',degin=True,degout=True):
    """
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
