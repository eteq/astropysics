#Erik Tollerud (etolleru@uci.edu) 2008

"""
This module stores tools for oberving (pre- and post-) as well as functioning as
a "miscellaneous" bin for various corrections and calculations that don't have
a better place to live.

The focus is currently on optical astronomy, as that is what the author does.
"""

from __future__ import division
from math import pi
import numpy as np

def extinction_correction(lineflux,linewl,EBmV,Rv=3.1,exttype='MW'):
    """
    Extinction correct a la Cardelli et al 89 from the supplied line data and
    a given E(B-V) along the sightline 
    
    inputs may be numpy arrays
    
    linewl is in angstroms, lineflux in erg s^-1 cm^-2
    
    if the input lineflux is None (or NaN but NOT False or 0) , Alambda is returned instead
    """
    from numpy import array,logical_and,logical_or,polyval,log10,isscalar,where
    from .phot import bandwl
    
    if exttype=='LMC' or exttype=='SMC':
        raise NotImplementedError 
    
    if type(linewl) == str:
        linewl=bandwl[linewl]
        
    if lineflux is None or lineflux is False:
        lineflux=np.nan
        
    
    lineflux=np.atleast_1d(lineflux).astype(float)
    linewl=np.atleast_1d(linewl)
    EBmV=np.atleast_1d(EBmV)
    n=np.max((lineflux.size,linewl.size,EBmV.size))
    
    if n!=1:
        if lineflux.size == 1:
            lineflux = np.ones(n)*lineflux[0]
        elif lineflux.size != n:
            raise ValueError("lineflux is incorrect length")
        if linewl.size == 1:
            linewl = np.ones(n)*linewl[0]
        elif linewl.size != n:
            raise ValueError("linewl is incorrect length")
        if EBmV.size == 1:
            EBmV = np.ones(n)*EBmV[0]
        elif EBmV.size != n:
            raise ValueError("EBmV is incorrect length")
        
    x=1e4/linewl #CCM x is 1/microns
    
    a=array(x)
    b=array(x)
    
    if any(logical_or(x<0.3,10<x)):
        raise ValueError('some wavelengths outside CCM 89 extinction curve range')
    
    irs=where(logical_and(0.3 <= x,x <= 1.1))
    opts=where(logical_and(1.1 <= x,x <= 3.3))
    nuv1s=where(logical_and(3.3 <= x,x <= 5.9))
    nuv2s=where(logical_and(5.9 <= x,x <= 8))
    fuvs=where(logical_and(8 <= x,x <= 10))
    
    #CCM Infrared
    a[irs]=.574*x[irs]**1.61
    b[irs]=-0.527*x[irs]**1.61
    
    #CCM NIR/optical
    a[opts]=polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
    b[opts]=polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)
    
    #CCM NUV
    y=x[nuv1s]-5.9
    Fa=-.04473*y**2-.009779*y**3
    Fb=-.2130*y**2-.1207*y**3
    a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)+Fa
    b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)+Fb
    
    a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)
    b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)
    
    #CCM FUV
    a[fuvs]=polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
    b[fuvs]=polyval((.374,-.42,4.257,13.67),x[fuvs]-8)
    
    AloAv=a+b/Rv #Al/Av
    Al=AloAv*Rv*EBmV #from Rv=Av/E(B-V)
    
    
    finalval=Al
    magi= ~np.isnan(lineflux)
    realmag=-2.5*log10(lineflux[magi])-Al[magi]
    finalval[magi]=10**(realmag/-2.5)
    
    if n==1:
        return finalval[0]
    else:
        return finalval
        

def extinction_from_flux_ratio(frobs,frexpect,outlambda=None,Rv=3.1,tol=1e-4):
    """
    frobs is the observed flux ratio f1/f2
    
    frexpect is either a string code specifying a hydrogen transition
    (e.g. 'Hab' is Halpha/Hbeta, 'Hde' is Hdelta/Hepsilon, from Ostriker 2E)    
    or a tuple of the form (expected f1/f2,lambda1,lambda2) wl in angstroms
    
    outputs E(B-V) to a tolerance specified by tol if outlambda is 0/False/None,
    otherwise outputs Alambda (tol still determines E(B-V)) (outlambda can be
    UBVRI or ugriz as from B&M)
    
    frobs can be an array, but other values cannot
    """
    from scipy.optimize import fmin
    
    scalarout=np.isscalar(frobs)
    
    frobs=np.atleast_1d(frobs)
    
    hd={
    'Hab':(2.86,6562.82,4861.33),
    'Hag':(6.16,6562.82,4340.46),
    'Had':(11.21,6562.82,4101.74),
    'Hae':(18.16,6562.82,3970.07),
    'Hbg':(2.15,4861.33,4340.46),
    'Hbd':(3.91,4861.33,4101.74),
    'Hbe':(6.33,4861.33,3970.07),
    'Hgd':(1.82,4340.46,4101.74),
    'Hge':(2.95,4340.46,3970.07),
    'Hde':(1.62,4101.74,3970.07)
    }
    if type(frexpect) == str:
        frexpect=hd[frexpect]
    
    fr0,l1,l2=frexpect
    
    A2mA1=-2.5*np.log10(fr0/frobs)
    
    fres = lambda x,Y: np.abs(Y - (extinction_correction(None,l2,x,Rv)-extinction_correction(None,l1,x,Rv)))
    
    EBmV=np.ndarray(A2mA1.shape)
    EBmVr=EBmV.ravel()
    for i,A in enumerate(A2mA1.ravel()):
        EBmVr[i]=fmin(fres,0,args=(A,),xtol=tol,disp=0)[0]
        
    if scalarout:
        EBmV=EBmV[0]
        
    if outlambda:
        if type(outlambda) == str:
            from .phot import bandwl
            outlambda=bandwl[outlambda]
        return extinction_correction(None,outlambda,EBmV,Rv)
    else:
        return EBmV
    
def get_dust_gal(dustmapfn,long,lat,interpolate=True):
    """
    gets map values from Schlegel, Finkbeiner, and Davis 1998 extinction maps
    
    input coordinates are in degrees
    
    if interpolate is an integer, it can be used to specify the order of the
    interpolating polynomial
    """
    from numpy import sin,cos,round,isscalar,array
    from pyfits import open
    if not dustmapfn.startswith('SFD') or not dustmapfn.endswith('.fits'):
        raise ValueError('file name not SFD98 data file')
    
    f=open(dustmapfn)
    try:
        mapd=f[0].data
    finally:
        f.close()
    assert mapd.shape[0] == mapd.shape[1],'maps not symmetrical - incorrect map file?'
    npix=mapd.shape[0]
    
    if isscalar(long):
        l=array([long])*pi/180
    else:
        l=array(long)*pi/180
    if isscalar(lat):
        b=array([lat])*pi/180
    else:
        b=array(lat)*pi/180
        
    if not len(l)==len(b):
        raise ValueError('input coordinate arrays are of different length')
        
    polename=dustmapfn.split('.')[0].split('_')[-1].lower()
    if polename=='ngp':
        n=1
    elif polename=='sgp':
        n=-1
    else:
        raise ValueError("couldn't determine South/North from filename")
    #project from galactic longitude/latitude to lambert pixels (see SFD98)
    x=npix/2*cos(l)*(1-n*sin(b))**0.5+npix/2-0.5
    y=-npix/2*n*sin(l)*(1-n*sin(b))**0.5+npix/2-0.5
    #now remap indecies - numpy arrays have y and x convention switched
    x,y=y,x
    
    if interpolate:
        from scipy.ndimage import map_coordinates
        if type(interpolate) is int:
            retval=map_coordinates(mapd,[x,y],order=interpolate)
        else:
            retval=map_coordinates(mapd,[x,y])
    else:
        x=round(x).astype(int)
        y=round(y).astype(int)
        retval=mapd[x,y]
        
    if isscalar(long) or isscalar(lat):
        return retval[0]
    else:
        return retval
        
    
def get_dust_radec(dustmapfn,ra,dec,interpolate=True):
    from .coords import equatorial_to_galactic_flt
    l,b=equatorial_to_galactic_flt(ra,dec,interpolate)
    return get_dust_gal(dustmapfn,l,b)