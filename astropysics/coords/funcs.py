#Copyright 2010 Erik Tollerud
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
This module contains functions for coordinate transforms and coordinate system
calculations.  It also includes distance-related calculations, including 
distances in expanding cosmologies.
   
Module API
^^^^^^^^^^

"""

from __future__ import division,with_statement
from ..constants import pi
import numpy as np

    
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
    from ..obstools import jd2000
    
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
    from ..obstools import jd2000
    
    d = jd - jd2000 #days since 2000
    res = (0.7790572732640 + 0.00273781191135448*d + (d%1.0))%1.0
    
    if degrees is None:
        return res
    elif degrees:    
        return res*360
    else:
        return res*2*pi
    
#polynomials for function greenwich_sidereal_time
from ..constants import asecperrad
_gmst_poly_circular179 = np.poly1d(np.array([-0.0000000368,
                                             -0.000029956,
                                             -0.00000044,
                                             1.3915817,
                                             4612.156534,
                                             0.014506])/asecperrad)
_gmst_poly_sofa = np.poly1d(np.array([0.00001882,
                                      -0.00009344,
                                      1.39667721,
                                      4612.15739966,
                                      0.014506])/asecperrad)
del asecperrad

def greenwich_sidereal_time(jd,apparent=True):
    """
    Computes the Greenwich Sidereal Time for a given Julian Date.
    
    :param jd: The Julian Date or a sequence of JDs, UT1
    :type jd: scalar or array-like
    :param apparent: 
        If True, the Greenwich Apparent Sidereal Time (GAST) is returned,
        using the method in the SOFA function iauGst00b, which
        computes nutation from the IAU 2000B nutation model, leaves out
        complementary terms in the equation of the equinox and uses
        UT1 instead of TT in the expression for GMST. In the special case that
        'simple' is given, a faster (but much lower precision) nutation model
        will be used. If False, the Greenwich Mean Sidereal Time (GMST) is
        returned, instead.
    :type apparent: 
    
    :returns: GMST or GAST in hours or an array of times (if `jd` is an array) 
        
    .. seealso:: 
        :func:`equation_of_the_equinoxes`, USNO Circular 179,  
        http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time, 
        IERS Technical Note No. 32 (esp. 5.10 Method (2B)), and SOFA functions iauGmst00 and iauGst00b 
    
    """
    from ..constants import asecperrad
    
    era = earth_rotation_angle(jd,False) #in radians
    
    t = (jd - 2451545.0)/36525
    #Why is Circular 179 different from SOFA? Using the SOFA value
    #gmst = era + _gmst_poly_circular179(t)
    gmst = era + _gmst_poly_sofa(t)
    
            
    if apparent:
        if apparent == 'simple':
            d = jd - 2451545.0
            eps =  np.radians(23.4393 - 0.0000004*d) #obliquity
            L = np.radians(280.47 + 0.98565*d) #mean longitude of the sun
            omega = np.radians(125.04 - 0.052954*d) #longitude of ascending node of moon
            dpsi = -0.000319*np.sin(omega) - 0.000024*np.sin(2*L) #nutation longitude
            coor = 0
        else:
            from .coordsys import _nutation_components2000B
            eps,dpsi,deps = _nutation_components2000B(jd,False)
            dpsi = dpsi
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
    Computes equation of the equinoxes GAST-GMST. That is, the difference
    between GMT computed using the mean equinox instead of the true equinox
    (i.e. including nutation).
    
    :param jd: The Julian Date or a sequence of JDs.
    :type jd: scalar or array-like
    
    :returns: the equation of the equinoxes for the provided date in hours.
    
    """
    return greenwich_sidereal_time(jd,True) - greenwich_sidereal_time(jd,False)

def equation_of_the_origins(jd):
    """
    Computes the equation of the origins ERA - GAST (ERA = Earth Rotation Angle,
    GAST = Greenwich Apparent Sidereal Time). This quantity is also the
    difference in RA between the Celestial Intermediate Origin and the Equinox.
    
    :param jd: The Julian Date or a sequence of JDs.
    :type jd: scalar or array-like
    
    :returns: the equation of the origins for the provided date in hours.
    
    """
    return earth_rotation_angle(jd,None)*24. - greenwich_sidereal_time(jd,True)

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
    computes the projected separation for a list of points in galacto-centric
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
    Compute the full 3D separation between two objects at distances `d1` and
    `d2` and angular positions `pos1` and `pos2`
    (:class:`~astropysics.coords.coordsys.LatLongCoordinates` objects, or an
    argument that will be used to generate a
    :class:`~astropysics.coords.coordsys.EquatorialCoordinatesEquinox` object)
    
    :param pos1: on-sky position of first object
    :type pos1: :class:`LatLongCoordinates` or initializer
    :param pos2: on-sky position of second object
    :type pos2: :class:`LatLongCoordinates` or initializer
    :param d1: distance to first object
    :type d1: scalar
    :param d2: distance to second object
    :type d2: scalar
    
    
    >>> from coordsys import LatLongCoordinates
    >>> p1 = LatLongCoordinates(0,0)
    >>> p2 = LatLongCoordinates(0,10)
    >>> '%.10f'%sky_sep_to_3d_sep(p1,p2,20,25)
    '6.3397355613'
    >>> '%.10f'%sky_sep_to_3d_sep('0h0m0s +0:0:0','10:20:30 +0:0:0',1,2)
    '2.9375007333'
        
    """    
    from .coordsys import LatLongCoordinates,EquatorialCoordinatesEquinox
    
    if not isinstance(pos1,LatLongCoordinates):
        pos1 = EquatorialCoordinatesEquinox(pos1)
    if not isinstance(pos2,LatLongCoordinates):
        pos2 = EquatorialCoordinatesEquinox(pos2)
        
    return (pos1-pos2).separation3d(d1,d2)

def radec_str_to_decimal(*args):
    """
    Convert a sequence of string coordinate specifiers to decimal degree arrays.
    
    Two input forms are accepted:
    
    * `radec_str_to_decimal(rastrs,decstrs)`
        In this form, `rastrs` and `decstrs` are sequences of strings with the
        RA and Dec, respectively.  
    * `radec_str_to_decimal(radecstrs)`
        In this form, `radecstrs` is a sequence of strings in any form accepted
        by the :class:`EquatorialCoordinatesBase` constructor. (typically
        canonical from like 17:43:54.23 +32:23:12.3)
    
    :returns: 
        (ras,decs) where `ras` and `decs` are  :class:`ndarrays <numpy.ndarray>`
        specifying the ra and dec in decimal degrees.
    
    """
    from .coordsys import AngularCoordinate,ICRSCoordinates
    from itertools import izip
    
    if len(args)==1:
        for s in args[0]:
            c = ICRSCoordinates(s)
            ras.append(c.ra.d)
            decs.append(c.dec.d)
    elif len(args)==2:
        ra,dec = args
        if len(ra) != len(dec):
            raise ValueError("length of ra and dec don't match")
        
        ras,decs=[],[]
        for r,d in izip(ra,dec):
            ras.append(AngularCoordinate(r,sghms=True).d)
            decs.append(AngularCoordinate(d,sghms=False).d)
    else:
        raise ValueError('radec_str_to_decimal only accepts (rastr,decstr) or (radecstr)')
    
    return np.array(ras),np.array(decs)

def match_coords(a1,b1,a2,b2,eps=1,mode='mask'):
    """
    Match one pair of coordinate :class:`arrays <numpy.ndarray>` to another
    within a specified tolerance (`eps`).
    
    Distance is determined by the cartesian distance between the two arrays,
    implying the small-angle approximation if the input coordinates are
    spherical. Units are arbitrary, but should match between all coordinates
    (and `eps` should be in the same units)
    
    :param a1: the first coordinate for the first set of coordinates
    :type a1: array-like
    :param b1: the second coordinate for the first set of coordinates
    :type b1: array-like
    :param a2: the first coordinate for the second set of coordinates
    :type a2: array-like
    :param b2: the second coordinate for the second set of coordinates
    :type b2: array-like
    :param eps: 
        The maximum separation allowed for coordinate pairs to be considered
        matched.
    :type eps: float
    :param mode:
        Determines behavior if more than one coordinate pair matches.  Can be:
        
        * 'mask'
            Returns a 2-tuple of boolean arrays (mask1,mask2) where `mask1`
            matches the shape of the first coordinate set (`a1` and `b1`), and
            `mask2` matches second set (`a2` and `b2`). The mask value is True
            if a match is found for that coordinate pair, otherwise, False.
        * 'maskexcept' 
            Retuns the same values as 'mask', and will raise an exception if
            more than one match is found.
        * 'maskwarn'
            Retuns the same values as 'mask', and a warning will be issued if
            more than one match is found. 
        * 'count'
            Returns a 2-tuple (nmatch1,nmatch2) with the number of objects that
            matched for each of the two sets of coordinate systems.
        * 'index'
            Returns a 2-tuple of integer arrays (ind1,ind2). `ind1` is a set of
            indecies into the first coordinate set, and `ind2` indexes the 
            second.  The two arrays match in shape and each element is the 
            index for a matched pair of coordinates - e.g. a1[ind1[i]] and 
            a2[ind2[i]] will give the "a" coordinate for a matched pair
            of coordinates.
        * 'match2D'
            Returns a 2-dimensional bool array. The array element M[i,j] is True
            if the ith coordinate of the first coordinate set
            matches the jth coordinate of the second set.
        * 'nearest'
            Returns (nearestind,distance,match). `nearestind` is an int array
            such that nearestind holds indecies into the *second* set of
            coordinates for the nearest object to the ith object in the first
            coordinate set (hence, it's shape matches the *first* coordinate
            set). `distance` is a float array of the same shape giving the
            corresponding distance, and `match` is a boolean array that is True
            if the distance is within `eps`, and is the same shape as the other
            outputs. Note that if a1 and b1 are the same object (and a2 and b2),
            this finds the second-closest match (because the first will always
            be the object itself if the coordinate pairs are the same) This mode
            is a wrapper around :func:`match_nearest_coords`.
    
    :returns: See `mode` for a description of return types.
    
    **Examples**
    
    >>> from numpy import array
    >>> ra1 = array([1,2,3,4])
    >>> dec1 = array([0,0,0,0])
    >>> ra2 = array([4,3,2,1])
    >>> dec2 = array([3.5,2.5,1.5,.5])
    >>> match_coords(ra1,dec1,ra2,dec2,1)
    (array([ True, False, False, False], dtype=bool), array([False, False, False,  True], dtype=bool))
    """
    
    identical = a1 is a2 and b1 is b2
    
    a1 = np.array(a1,copy=False).ravel()
    b1 = np.array(b1,copy=False).ravel()
    a2 = np.array(a2,copy=False).ravel()
    b2 = np.array(b2,copy=False).ravel()
    
    #bypass the rest for 'nearest', as it calls match_nearest_coords
    if mode == 'nearest':
        #special casing so that match_nearest_coords dpes second nearest
        if identical: 
            t = (a1,b1)
            seps,i2 = match_nearest_coords(t,t)
        else:
            seps,i2 = match_nearest_coords((a1,b1),(a2,b2))
        return i2,seps,(seps<=eps)
        
    def find_sep(A,B):
        At = np.tile(A,(len(B),1))
        Bt = np.tile(B,(len(A),1))
        return At.T-Bt
    sep1 = find_sep(a1,a2)
    sep2 = find_sep(b1,b2)
    sep = np.hypot(sep1,sep2)
    
    matches =  sep <= eps
    
    if mode == 'mask':
        return np.any(matches,axis=1),np.any(matches,axis=0)
    elif mode == 'maskexcept':
        s1,s2 = np.sum(matches,axis=1),np.sum(matches,axis=0) 
        if np.all(s1<2) and np.all(s2<2):
            return s1>0,s2>0
        else:
            raise ValueError('match_coords found multiple matches')
    elif mode == 'maskwarn':
        s1,s2 = np.sum(matches,axis=1),np.sum(matches,axis=0) 
        from warnings import warn
        
        for i in np.where(s1>1)[0]:
            warn('1st index %i has %i matches!'%(i,s1[i]))
        for j in np.where(s2>1)[0]:
            warn('2nd index %i has %i matches!'%(j,s2[j]))
        return s1>0,s2>0
    elif mode == 'count':
        return np.sum(np.any(matches,axis=1)),np.sum(np.any(matches,axis=0)) 
    elif mode == 'index':
        return np.where(matches)
    elif mode == 'match2D':
        return matches.T
    elif mode == 'nearest':
        assert False,"'nearest' should always return above this - code should be unreachable!"
    else:
        raise ValueError('unrecognized mode')
    
def match_nearest_coords(c1,c2=None,n=None):
    """
    Match a set of coordinates to their nearest neighbor(s) in another set of
    coordinates.
    
    :param c1: 
        A D x N array with coordinate values (either as floats or
        :class:`AngularPosition` objects) or a sequence of
        :class:`LatLongCoordinates` objects for the first set of coordinates.
    :param c2: 
        A D x N array with coordinate values (either as floats or
        :class:`AngularPosition` objects) or a sequence of
        :class:`LatLongCoordinates` objects for the second set of coordinates.
        Alternatively, if this is None, `c2` will be set to `c1`, finding the 
        nearest neighbor of a point in `c1` to another point in `c1`.
    :param int n: 
        Specifies the nth nearest neighbor to be returned (1 means the closest
        match). If None, it will default to 2 if `c1` and `c2` are the same
        object (just equality is not enough - they must actually be the same
        in-memory array), or 1 otherwise. This is because if `c1` and `c2` are
        the same, a coordinate matches to *itself* instead of the nearest other
        coordinate.
    
    :returns: 
        (seps,ind2) where both are arrays matching the shape of `c1`. `ind2` is
        indecies into `c2` to find the nearest to the corresponding `c1`
        coordinate, and `seps` are the distances.
    """
    try:
        from scipy.spatial import cKDTree as KDTree
    except ImportError:
        from warnings import warn
        warn('C-based scipy kd-tree not available - match_nearest_coords will be much slower!')
        from scipy.spatial import KDTree
        
    if c2 is None:
        c2 = c1
    if n is None:    
        n = 2 if c1 is c2 else 1
        
    c1 = np.array(c1,ndmin=1,copy=False)
    c2 = np.array(c2,ndmin=1,copy=False)
    
    if len(c1.shape)==1:
        a1 = np.empty(c1.size)
        b1 = np.empty(c1.size)
        a2 = np.empty(c1.size)
        b2 = np.empty(c1.size)
        for i in range(len(c1)):
            a1[i] = c1[i].long.d
            b1[i] = c1[i].lat.d
            a2[i] = c2[i].long.d
            b2[i] = c2[i].lat.d
        c1 = np.array((a1,b1))
        c2 = np.array((a2,b2))
    elif len(c1.shape)!=2:
        raise ValueError('match_nearest_coords inputs have incorrect number of dimensions')
    
    if c1.shape[0] != c2.shape[0]:
        raise ValueError("match_nearest_coords inputs don't match in first dimension")
    
    kdt = KDTree(c2.T)
    if n==1:
        return kdt.query(c1.T)
    else:
        dist,inds = kdt.query(c1.T,n)
        return dist[:,n-1],inds[:,n-1]
        
        
    
def separation_matrix(v,w=None,tri=False):
    """
    Computes a matrix of the separation between each of the components of the
    first dimension of an array. That is, A[i,j] = v[i]-w[j]. 
    
    :param v: The first array with first dimension n
    :param w: 
        The second array with first dimension m, and all following dimensions
        matched to `v`. If None, `v` will be treated as `w` (e.g. the separation
        matrix of `v` with itself will be returned).
    :param bool tri: 
        If True, the lower triangular part of the matrix is set to 0 (this is
        really only useful if w is None) 
        
    :returns: 
        Separation matrix with dimension nXmX(whatever the remaining dimensions
        are)
        
    .. seealso::
    
        :mod:`scipy.spatial.distance`, in particular the
        :func:`scipy.spatial.distance.pdist` function. It is much more efficient
        and flexible at computing distances if individual components and sign
        information is unnecessary.
        
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
     
    >>> from astropysics.constants import choose_cosmology
    >>> cosmo = choose_cosmology('wmap7baoh0')
    
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
    
    from ..constants import H0,omegaM,omegaL,omegaR,c
    
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
            if np.any(interr/intres > inttol):
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
            from ..phot import distance_modulus
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
    from ..constants import get_cosmology
    c = get_cosmology()
    if zerr is None:
        return c.H(z)
    else:
        H=c.H(z)
        upper=c.H(z+zerr)
        lower=c.H(z-zerr)
        return H,upper-H,lower-H

def angular_to_physical_size(angsize,zord,usez=False,**kwargs):
    """
    Converts an observed angular size (in arcsec or as an AngularSeparation 
    object) to a physical size.
    
    :param angsize: Angular size in arcsecond.
    :type angsize: float or an :class:`AngularSeparation` object
    :param zord: Redshift or distance
    :type zord: scalar number
    :param usez:
        If True, the input will be interpreted as a redshift, and kwargs
        will be passed into the distance calculation. The result will be in
        pc. Otherwise, `zord` will be interpreted as a distance.
    :type usez: boolean
    
    kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
    
    :returns: 
        A scalar value for the physical size (in pc if redshift is used,
        otherwise in `zord` units)
    """
    from ..constants import asecperrad
    
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
    as an AngularSeparation object if objout is True)
    
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
        If True, return value is an :class:`AngularSeparation` object,
        otherwise, angular size in arcsec.
    :type: bool
    
    kwargs are passed into :func:`cosmo_z_to_dist` if `usez` is True.
    
    :returns: 
        The angular size in acsec, or an :class:`AngularSeparation` object if
        `objout` is True.
        
    """
    from ..constants import asecperrad
    
    if usez:
        d = cosmo_z_to_dist(zord,disttype=2,**kwargs)*1e6 #pc
    else:
        if len(kwargs)>0:
            raise TypeError('if not using redshift, kwargs should not be provided')
        d = zord
        
    r = physize
    res = asecperrad*np.arcsin(r*(d*d+r*r)**-0.5)
    
    if objout:
        return AngularSeparation(res/3600)
    else:
        return res
    

def geographic_to_geocentric_latitude(geoglat):
    """
    Converts a geographic/geodetic latitude to a geocentric latitude.
    
    :param geoglat:
        An :class:`astropysics.coords.AngularCoordinate` object (or arguments to
        create one) or an angle in degrees for the geographic latitude.
        
    :returns: 
        An :class:`astropysics.coords.AngularCoordinate` object with the
        geocentric latitude.
    """
    from astropysics.constants import Rea,Reb
    from astropysics.coords import AngularCoordinate
    from operator import isSequenceType
    
    if not isinstance(geoglat,AngularCoordinate):
        if isSequenceType(geoglat):
            rads = AngularCoordinate(*geoglat).radians
        else:
            rads = AngularCoordinate(geoglat).radians
    else:
        rads = geoglat.radians
    
    boasq = (Reb/Rea)**2
    return AngularCoordinate(np.arctan(boasq*np.tan(rads)),radians=True)

def geocentric_to_geographic_latitude(geoclat):
    """
    Converts a geocentric latitude to a geographic/geodetic latitude.
    
    :param geoclat:
        An :class:`astropysics.coords.AngularCoordinate` object (or arguments to
        create one) or an angle in degrees for the geocentric latitude.
        
    :returns: 
        An :class:`astropysics.coords.AngularCoordinate` object with the
        geographic latitude.
    """
    from astropysics.constants import Rea,Reb
    from astropysics.coords import AngularCoordinate
    from operator import isSequenceType
    
    if not isinstance(geoclat,AngularCoordinate):
        if isSequenceType(geoclat):
            rads = AngularCoordinate(*geoclat).radians
        else:
            rads = AngularCoordinate(geoclat).radians
    else:
        rads = geoclat.radians
        
    boasq = (Reb/Rea)**2
    return AngularCoordinate(np.arctan((1/boasq)*np.tan(rads)),radians=True)

    
#<---------------------DEPRECATED transforms----------------------------------->

#These will be removed eventually

#galactic coordate reference positions from IAU 1959
#from astropysics.coords.coordsys import EquatorialCoordinatesEquinox as _EquatorialCoordinatesEquinox
#_galngpJ2000 = _EquatorialCoordinatesEquinox('12h51m26.282s','+27d07m42.01s')
#_galngpB1950 = _EquatorialCoordinatesEquinox('12h49m0s','27d24m0s')
#_gall0J2000=122.932
#_gall0B1950=123

#def celestial_transforms(ai,bi,transtype=1,epoch='J2000',degin=True,degout=True):
#    """
#    :deprecated:
    
#    Use this to transform between Galactic,Equatorial, and Ecliptic coordinates
    
#    transtype can be a number from the table below, or 'ge','eg','gq','qg','gc',
#    'cg','cq','qc'
    
#    transtype   From           To       |  transtype    From         To
#        1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec    
#        2     Galactic       RA-DEC     |     5       Ecliptic      Galactic  
#        3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
        
#    adapted from IDL procedure EULER 
#    (http://astro.uni-tuebingen.de/software/idl/astrolib/astro/euler.html)
#    """
#    #   J2000 coordinate conversions are based on the following constants
#    #   (see the Hipparcos explanatory supplement).
#    #  eps = 23.4392911111d              Obliquity of the ecliptic
#    #  alphaG = 192.85948d               Right Ascension of Galactic North Pole
#    #  deltaG = 27.12825d                Declination of Galactic North Pole
#    #  lomega = 32.93192d                Galactic longitude of celestial equator  
#    #  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
#    #  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
#    #  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator   
    
#    from warnings import warn
#    warn('celestial_transforms function is deprecated - use general coordinate transform framework',DeprecationWarning)           

#    if epoch == 'B1950':
#            psi   = ( 0.57595865315, 4.9261918136,0, 0,0.11129056012, 4.7005372834)     
#            stheta =( 0.88781538514,-0.88781538514, 0.39788119938,-0.39788119938, 0.86766174755,-0.86766174755)    
#            ctheta =( 0.46019978478, 0.46019978478,0.91743694670, 0.91743694670, 0.49715499774, 0.49715499774)    
#            phi =   ( 4.9261918136,  0.57595865315,  0, 0,  4.7005372834, 0.11129056012)

#    elif epoch == 'J2000':
#            psi   = ( 0.57477043300,4.9368292465,0,0,0.11142137093, 4.71279419371)     
#            stheta =( 0.88998808748,-0.88998808748,  0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025)   
#            ctheta =( 0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207,  0.49714719172, 0.49714719172)    
#            phi  =  ( 4.9368292465,  0.57477043300,  0, 0,        4.71279419371, 0.11142137093)
#    else:
#            raise ValueError('unknown epoch')
            
#    from math import pi
#    from numpy import array,sin,cos,arcsin,arctan2
#    twopi   =   2.0*pi
#    fourpi  =   4.0*pi
#    deg_to_rad = 180.0/pi
    
    
#    if degin:
#        ai,bi=array(ai),array(bi)
#    else:
#        ai,bi=np.degrees(ai),np.degrees(bi)
    
#    if type(transtype) == int:
#        i = transtype - 1
#    else:
#        transd={'ge':1,'eg':0,'gq':1,'qg':0,'gc':5,'cg':4,'cq':3,'qc':2}
#        i  = transd[transtype]
#    a  = ai/deg_to_rad - phi[i]
#    b = bi/deg_to_rad
#    sb = sin(b) 
#    cb = cos(b)
#    cbsa = cb * sin(a)
#    b  = -stheta[i] * cbsa + ctheta[i] * sb
#    try:
#            b[b>1.0]=1.0
#    except TypeError: #scalar
#            if b > 1:
#                    b=array(1.0)
#    bo = arcsin(b)*deg_to_rad
#    a =  arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a) )
#    ao = ( (a+psi[i]+fourpi) % twopi) * deg_to_rad
    
#    if not degout:
#        ao,bo = np.radians(ao),np.radians(bo)
    
#    return ao,bo

#_B1950toJ2000xyz=np.matrix([[0.999926,  -0.011179,  -0.004859],
#                            [0.011179,   0.999938,  -0.000027],
#                            [0.004859,   0.000027,   0.999988]])

#def epoch_transform(ra,dec,inepoch='B1950',outepoch='J2000',degrees=True):
#    """
#    :deprecated:
#    """
#    from warnings import warn
#    warn('epoch_transform function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
#    if inepoch != 'B1950' and inepoch != 'J2000':
#        raise ValueError('unrecognized epoch '+inepoch)
#    if outepoch != 'B1950' and outepoch != 'J2000':
#        raise ValueError('unrecognized epoch '+outepoch)
#    if degrees:
#        ra,dec=np.radians(ra),np.radians(dec)
#    else:
#        ra,dec=np.array(ra),np.array(dec)
    
#    if inepoch == outepoch:
#        trans=np.matrix(np.eye(3))
#    elif inepoch == 'B1950' and outepoch == 'J2000':
#        trans=_B1950toJ2000xyz
#    elif inepoch == 'J2000' and outepoch == 'B1950':
#        trans=_B1950toJ2000xyz.I
#    else:
#        raise ('unrecognized epochs')
    
#    x=np.cos(ra)*np.cos(dec)
#    y=np.sin(ra)*np.cos(dec)
#    z=np.sin(dec)
    
#    v=np.matrix((x,y,z))
#    xp,yp,zp=trans*v
    
#    rap=np.arctan2(yp,xp)
#    decp=np.arcsin(zp)
    
#    return rap,decp

#def galactic_to_equatorial(l,b,epoch='J2000',strout=None):
#    """
#    :deprecated:
    
#    convinience function for celestial_transforms
#    if strout is None, will automatically decide based on inputs
#    """
#    from warnings import warn
#    warn('galactic_to_equatorial function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
#    from operator import isSequenceType
    
#    if type(l) == str:
#        l=AngularCoordinate(l).degrees
#        if strout is None:
#            strout=True
#    if type(b) == str:
#        b=AngularCoordinate(b).degrees
#        if strout is None:
#            strout=True
#    ra,dec = celestial_transforms(l,b,transtype='ge',epoch=epoch)
#    if strout:
#        if not isSequenceType(ra):
#            ra=[ra]
#        if not isSequenceType(dec):
#            dec=[dec]
#        rao,deco=[],[]
#        for rai in ra:
#            rao.append(AngularCoordinate(rai).getHmsStr())
#        for deci in dec:
#            deco.append(AngularCoordinate(deci).getDmsStr())
#        return rao,deco
#    else:
#        return ra,dec
    
    
#def equatorial_to_galactic(ra,dec,epoch='J2000',strout=None):
#    """
#    :deprecated:
    
#    convinience function for celestial_transforms
#    if strout is None, will automatically decide based on inputs
#    """
#    from warnings import warn
#    warn('equatorial_to_galactic function is deprecated - use general coordinate transform framework',DeprecationWarning)
    
#    from operator import isSequenceType
    
#    if type(ra) == str:
#        ra=AngularCoordinate(ra).degrees
#        if strout is None:
#            strout=True
#    if type(dec) == str:
#        dec=AngularCoordinate(dec).degrees
#        if strout is None:
#            strout=True
#    l,b = celestial_transforms(ra,dec,transtype='eg',epoch=epoch)
#    if strout:
#        if not isSequenceType(l):
#            l=[l]
#        if not isSequenceType(b):
#            b=[b]
#        lo,bo=[],[]
#        for li in l:
#            lo.append(AngularCoordinate(li).getDmsStr())
#        for bi in b:
#            bo.append(AngularCoordinate(bi).getDmsStr())
#        return lo,bo
#    else:
#        return l,b
    