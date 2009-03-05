#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains objects and functions for specifying locations as well 
as calculating distances and similar tools.

Some of the calculations involved make use of the currently selected cosmology 
(see astropysics.constants) and hence may not function properly if a 
particularly strange cosmology is in use.

Note that some of the functions use the PyEphem (http://rhodesmill.org/pyephem/)
project for conversions and calculations - some functions will work poorly or 
not at all without this package installed
"""

from __future__ import division
from math import pi
import numpy as np


try:
    import ephem
except ImportError:
    from warnings import warn
    warn('PyEphem not found - some astropysics.coords functions will not work correctly')
    ephem = None

#<----------------coordinate classes and related functions------------------>


class AngularCoordinate(object):
    import re as _re
    __slots__=('__decval')
    
    def __setdegminsec(self,dms):
        if not hasattr(dms, '__iter__') or len(dms)!=3:
            raise ValueError('Must set as a length-3 iterator')
        self.degrees=abs(dms[0])+abs(dms[1])/60.+abs(dms[2])/3600.
        if dms[0]<0:
            self.__decval*=-1
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
        self.__decval=deg*pi/180.
    def __getdegdec(self):
        return self.__decval*180/pi
    
    def __setrad(self,deg):
        self.__decval=deg
    def __getrad(self):
        return self.__decval
    
    def __sethrdec(self,hr):
        self.__decval=hr*pi/12
    def __gethrdec(self):
        return self.__decval*12/pi
    
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
    __hmsre=_re.compile(r'.*?(\d{1,2})(?:h|hr)\s*(\d{1,2})(?:m|min)\s*(\d+(?:\.?\d*))(?:s|sec).*')
    __puredre=_re.compile(r'.*?([+-]?\s*\d+(?:\.?\d+))(?:d|deg).*')
    #__dmsre=_re.compile(r'.*?([+-]?\s*\d{1,2})(?:d|deg)\s*(\d{1,2})(?:m|min)\s*(\d+(?:\.?\d+))(?:s|sec).*')
    __dmsre=_re.compile(r'.*?([+-])?(\d{1,2})(?:d|deg)\s*(\d{1,2})(?:m|min)\s*(\d+(?:\.?\d*))(?:s|sec).*')
    __sexre=_re.compile(r'.*?(\+|\-)?(\d{1,3}):(\d{1,2}):(\d+(?:.\d+)?).*')
    def __init__(self,inpt=None,sghms=None):
        """
        If an undecorated 3-element iterator, inpt is taken to be deg,min,sec, 
        othewise, input is cast to a float and treated as decimal degrees
        
        if sexigesimal, input will be interpreted as h:m:s if sghms
        is True, or d:m:s if sghms is False.  If None, a +/- will indicate
        d:m:s and nothing indicates r:m:s
        """
        if inpt.__class__.__name__=='AngularCoordinate':
            self.__decval=inpt.__decval
        elif isinstance(inpt,basestring):
            sexig=self.__sexre.match(inpt)
            hm=self.__purehre.match(inpt)
            hmsm=self.__hmsre.match(inpt)
            dm=self.__puredre.match(inpt)
            dmsm=self.__dmsre.match(inpt)
            if sexig:
                t=sexig.group(2,3,4)
                if sghms is None:
                    if sexig.group(1) is None:
                        self.hrsminsec=int(t[0]),int(t[1]),float(t[2])
                    else:
                        sgn = 1 if sexig.group(1) == '+' else -1
                        self.degminsec=sgn*int(t[0]),int(t[1]),float(t[2])
                else:
                    sgn = -1 if sexig.group(1) == '-' else 1
                    if sghms:
                        self.hrsminsec=sgn*int(t[0]),int(t[1]),float(t[2])
                    else:
                        self.degminsec=sgn*int(t[0]),int(t[1]),float(t[2]) 
            elif hmsm:
                t=hmsm.group(1,2,3)
                self.hrsminsec=int(t[0]),int(t[1]),float(t[2])
            elif hm:
                self.hours=float(hm.group(1))
            elif dmsm:
                sgn = -1 if dmsm.group(1) =='-' else 1
                t=dmsm.group(2,3,4)
                self.degminsec=sgn*int(t[0]),int(t[1]),float(t[2])
            elif dm:
                self.degrees=float(hm.group(1))
            else:
                raise Exception('Unrecognized string format')
            
        elif hasattr(inpt,'__iter__') and len(inpt)==3:
            self.degminsec=inpt
        elif inpt is None:
            self.__decval=0
        else:
            self.__decval=float(inpt)*pi/180.
            
            
    def __eq__(self,other):
        if type(other) == AngularCoordinate:
            return self.__decval==other.__decval
        else:
            return self.__decval==other
    def __ne__(self,other):
        return not self.__eq__(other)
        
    def __add__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval+other.__decval)
        else:
            return AngularCoordinate(self.__decval+other)
    def __sub__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval-other.__decval)
        else:
            return AngularCoordinate(self.__decval-other)
        
    def __mul__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval*other.__decval)
        else:
            return AngularCoordinate(self.__decval*other)
        
    def __div__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval/other.__decval)
        else:
            return AngularCoordinate(self.__decval/other)
        
    def __truediv__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval/other.__decval)
        else:
            return AngularCoordinate(self.__decval/other)
        
    def __pow__(self,other):
        if type(other) == AngularCoordinate:
            return AngularCoordinate(self.__decval**other.__decval)
        else:
            return AngularCoordinate(self.__decval**other)
        
    def __float__(self):
        return self.degrees
        
        
    def getDmsStr(self,secform = None,sep = None, sign=True, canonical=False):
        """
        gets the string representation of this AngularCoordinate as degrees,
        minutes, and seconds
        
        secform is the formatter for the seconds component
        
        sep is the seperator between components - defaults to degree symbol,
        ' and " symbols 
        
        sign forces sign to be present before degree component
        
        canonical forces [+/-]dd:mm:ss.ss , overriding other arguments
        """
        d,m,s = self.degminsec
        
        if canonical:
            return '%0+3.i:%02.i:%05.2f'%(d,m,s)
        
        d,m=str(d),str(m)
        
        if secform is None:
            s = str(s)
        else:
            s = secform%s
        
        tojoin = []
        
        if sign:
            d='+'+d if d >= 0 else d
        
        if d is not '0':
            tojoin.append(d)
            if sep is None:
                tojoin.append(unichr(176).encode("latin-1"))
                
        if m is not '0':
            tojoin.append(m)
            if sep is None:
                tojoin.append("'")
                
        tojoin.append(s)
        if sep is None:
            tojoin.append('"')
            
        if sep is None:
            sep = ''
        return sep.join(tojoin)
        
    def getHmsStr(self,secform = None,sep = None, canonical = False):
        """
        gets the string representation of this AngularCoordinate as hours,
        minutes, and seconds
        
        secform is the formatter for the seconds component
        
        sep is the seperator between components - defaults to h, m, and s
        
        canonical forces [+/-]dd:mm:ss.ss , overriding other arguments
        """
        
        h,m,s = self.hrsminsec
        
        if canonical:
            return '%02.i:%02.i:%05.2f'%(h,m,s)
        
        h,m=str(h),str(m)
        if secform is None:
            s = str(s)
        else:
            s = secform%s
        
        tojoin = []
        
        tojoin.append(h)
        if sep is None:
            tojoin.append('h')
            
        tojoin.append(m)
        if sep is None:
            tojoin.append('m')
#        if h is not '0':
#            tojoin.append(h)
#            if sep is None:
#                tojoin.append('h')
                
#        if m is not '0':
#            tojoin.append(m)
#            if sep is None:
#                tojoin.append('m')
                
        tojoin.append(s)
        if sep is None:
            tojoin.append('s')
            
        if sep is None:
            sep = ''
        return sep.join(tojoin)
    
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
    
class AngularPosition(object):
    __slots__=('__ra','__dec','__raerr','__decerr','__epoch')
    
    def __init__(self,*args,**kwargs):
        """
        can either specify kwargs ra,dec,raerr,decerr,epoch or follow form:
        AngularPosition() (default)
        AngularPosition(AngularPosition | 'ra dec')
        AngularPosition(ra,dec)
        AngularPosition(ra,dec,epoch)
        AngularPosition(ra,dec,raerr,decerr)
        AngularPosition(ra,dec,raerr,decerr,epoch)
        """
        if len(args)  == 0:
            pass
        elif len(args) == 1:
            if isinstance(args[0],AngularPosition):
                kwargs['ra'] = AngularCoordinate(args[0].__ra)
                kwargs['dec'] = AngularCoordinate(args[0].__dec)
                kwargs['raerr'] = AngularCoordinate(args[0].__raerr)
                kwargs['decerr'] = AngularCoordinate(args[0].__decerr)
                kwargs['epoch'] = args[0].__epoch
            aspl = args[0].split()
            #TODO: make smarter
            kwargs['ra'] = AngularCoordinate(aspl[0])
            kwargs['dec'] = AngularCoordinate(aspl[1],sghms=False)
        elif len(args) == 2:
            kwargs['ra'] = AngularCoordinate(args[0])
            kwargs['dec'] = AngularCoordinate(args[1],sghms=False)
        elif len(args) == 3:
            kwargs['ra'] = AngularCoordinate(args[0])
            kwargs['dec'] = AngularCoordinate(args[1],sghms=False)
            kwargs['epoch'] = args[2]
        elif len(args) == 4:
            kwargs['ra'] = AngularCoordinate(args[0])
            kwargs['dec'] = AngularCoordinate(args[1],sghms=False)
            kwargs['raerr'] = AngularCoordinate(args[2])
            kwargs['decerr'] = AngularCoordinate(args[3])
        elif len(args) == 5:
            kwargs['ra'] = AngularCoordinate(args[0])
            kwargs['dec'] = AngularCoordinate(args[1],sghms=False)
            kwargs['raerr'] = AngularCoordinate(args[2])
            kwargs['decerr'] = AngularCoordinate(args[3])
            kwargs['epoch'] = args[4]
        else:
            raise ValueError('Unrecognized format for coordiantes')
        
        self.__ra=kwargs.pop('ra',AngularCoordinate(0))
        self.__dec=kwargs.pop('dec',AngularCoordinate(0))
        self.__raerr=kwargs.pop('raerr',AngularCoordinate(0))
        self.__decerr=kwargs.pop('decerr',AngularCoordinate(0))
        self.__epoch=kwargs.pop('epoch','J2000')
        if len(kwargs) > 0:
            raise ValeError('unrecognized keyword'+'s ' if len(kwargs)> 1 else ' '+','.join(kwargs.keys()))
        
    
    def __setra(self,ra):
        if type(ra) == AngularCoordinate:
            self.__ra=ra
        else:
            self.__ra=AngularCoordinate(ra)
    
    def __getra(self):
        return self.__ra
    
    def __setdec(self,dec):
        if type(dec) == AngularCoordinate:
            self.__dec=dec
        else:
            self.__dec=AngularCoordinate(dec)
    
    def __getdec(self):
        return self.__dec
    
    def __setraerr(self,raerr):
        if type(raerr) == AngularCoordinate:
            self.__raerr=raerr
        else:
            self.__raerr=AngularCoordinate(raerr)
    
    def __getraerr(self):
        return self.__raerr
    
    def __setdecerr(self,decerr):
        if type(decerr) == AngularCoordinate:
            self.__decerr=decerr
        else:
            self.__dec=AngularCoordinate(decerr)
    
    def __getdecerr(self):
        return self.__decerr
    
    def __setepoch(self,epoch):
        ra,dec=epoch_transform(self.__ra.degrees,self.__dec.degrees,self.__epoch,epoch)
        self.__ra=AngularCoordinate(ra)
        self.__dec=AngularCoordinate(ra)
        #TODO:should do some sort of error propogation, but that only matters for large precessions
        self.__epoch=epoch
        
    def __getepoch(self):
        return self.__epoch
    
    def __str__(self):
        rastr=self.__ra.getHmsStr()+('\\pm'+self.__raerr.getHmsStr() if self.__raerr != 0 else '')
        decstr=self.__dec.getDmsStr()+('\\pm'+self.__decerr.getDmsStr() if self.__decerr != 0 else '')
        return ' '.join((rastr,decstr))
    
    ra=property(fget=__getra,fset=__setra)
    dec=property(fget=__getdec,fset=__setdec)
    raerr=property(fget=__getraerr,fset=__setraerr)
    decerr=property(fget=__getdecerr,fset=__setdecerr)
    epoch=property(fget=__getepoch,fset=__setepoch)
    e=epoch
    
    def __correctedra(self):
        from math import cos
        return self.ra.radians * cos(self.dec.radians)
    
    def __eq__(self,other):
        from operator import isSequenceType
        if type(other) == AngularPosition:
            return self.__ra==other.__ra and  self.__dec==other.__dec
        elif isSequenceType(other):
            return self.__ra==other[0] and  self.__dec==other[1]
        else:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)
    
    #TODO:add sequence operations
    def __add__(self,other):
        return AngularPosition(self.__ra+other,self.__dec+other)
    def __sub__(self,other):
        from math import degrees
        print type(other)
        if type(other) == AngularPosition:
            dcra=self.__correctedra()-other.__correctedra()
            ddec=self.__dec.radians-other.__dec.radians
            res=AngularCoordinate(degrees((dcra*dcra+ddec*ddec)**0.5))
            if res.degrees > 1:
                print 'WARNING: small-angle approximation incorrect'
            return res
        return AngularPosition(self.__ra-other,self.__dec-other)
    def __mul__(self,other):
        return AngularPosition(self.__ra*other,self.__dec*other)
    def __div__(self,other):
        return AngularPosition(self.__ra/other,self.__dec/other)
    def __truediv__(self,other):
        return AngularPosition(self.__ra/other,self.__dec/other)
    def __pow__(self,other):
        return AngularPosition(self.__ra**other,self.__dec**other)

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

def seperation3d(d1,d2,ap1,ap2):
    from numpy import sin,cos
    if type(ap1) != AngularPosition:
        ap1=AngularPosition(ap1)
    if type(ap2) != AngularPosition:
        ap2=AngularPosition(ap2)
    
    theta1,phi1=(pi/2-ap1.dec.radians,ap1.ra.radians)
    theta2,phi2=(pi/2-ap2.dec.radians,ap2.ra.radians)
    
    from math import degrees
    dx=d2*sin(theta2)*cos(phi2)-d1*sin(theta1)*cos(phi1)
    dy=d2*sin(theta2)*sin(phi2)-d1*sin(theta1)*sin(phi1)
    dz=d2*cos(theta2)-d1*cos(theta1)
    
    return (dx*dx+dy*dy+dz*dz)**0.5
    


def galactic_to_equatorial(l,b,epoch='J2000',strout=None):
    """
    convinience function for celestial_transforms
    if strout is None, will automatically decide based on inputs
    """
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


#<-------------basic transforms-------------------->
def cartesian_to_spherical(x,y,z,degrees=False):
    """
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
    if degrees is true, converts from degrees to radians for theta and phi
    returns x,y,z in PHYSICIST convention - (1,0,pi/2) is x-axis
    """
    if degrees:
        t,p=np.radians(t),np.radians(p)
    x=r*np.sin(t)*np.cos(p)
    y=r*np.sin(t)*np.sin(p)
    z=r*np.cos(t)
    
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



#galactic coordate reference positions from IAU 1959 and wikipedia
_galngpJ2000=AngularPosition('12h51m26.282s','+27d07m42.01s')
_galngpB1950=AngularPosition('12h49m0s','27d24m0s')
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
    if inepoch != 'B1950' and inepoch != 'J2000':
        raise ValueError('unrecognized epoch '+inepoch)
    if outepoch != 'B1950' and outepoch != 'J2000':
        raise ValueError('unrecognized epoch '+outepoch)
    if degrees:
        ra,dec=np.radians(ra),np.radians(dec)
    else:
        ra,dec=np.array(ar),np.array(dec)
    
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
    
    v=matrix((x,y,z))
    xp,yp,zp=trans*v
    
    rap=np.atan2(yp,xp)
    decp=np.arcsin(zp)
    
    return rap,decp

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

def angular_to_physical_size(angsize,z,**kwargs):
    """
    converts an observed angular size (in arcsec) to a physical size (in kpc)
    
    (assumes small angle approximation)
    
    kwargs go into cosmo_z_to_dist
    """
    d=cosmo_z_to_dist(z,disttype=2,**kwargs)*1e3 #kpc
    return angsize*d/206265

def physical_to_angular_size(physize,z,**kwargs):
    """
    converts a physical size (in kpc) to an observed angular size (in arcsec)
    
    (assumes small angle approximation)
    
    kwargs go into cosmo_z_to_dist
    """
    d=cosmo_z_to_dist(z,disttype=2,**kwargs)*1e3 #kpc
    return physize*206265/d