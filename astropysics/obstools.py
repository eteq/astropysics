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

#<----------------------Extinction and dust-related---------------------------->

class Extinct(object):
    def  __init__(self):
        raise NotImplementedError

def get_SFD_dust(long,lat,dustmap='ebv',interpolate=True):
    """
    Gets map values from Schlegel, Finkbeiner, and Davis 1998 extinction maps
    
    dustmap can either be a filename (if '%s' appears in the string, it will be
    replaced with 'ngp' or 'sgp'), or one of:
        i100: 100-micron map in MJy/Sr
        x   : X-map, temperature-correction factor
        t   : Temperature map in degrees Kelvin for n=2 emissivity
        ebv : E(B-V) in magnitudes
        mask: Mask values
    
    input coordinates are in degrees of galactic latiude and logitude - they can
    be scalars or arrays
    
    if interpolate is an integer, it can be used to specify the order of the
    interpolating polynomial
    """
    from numpy import sin,cos,round,isscalar,array,ndarray
    from pyfits import open
    
    if type(dustmap) is not str:
        raise ValueError('dustmap is not a string')
    dml=dustmap.lower()
    if dml == 'ebv' or dml == 'eb-v' or dml == 'e(b-v)' :
        dustmapfn='SFD_dust_4096_%s.fits'
    elif dml == 'i100':
        dustmapfn='SFD_i100_4096_%s.fits'
    elif dml == 'x':
        dustmapfn='SFD_xmap_%s.fits'
    elif dml == 't':
        dustmapfn='SFD_temp_%s.fits'
    elif dml == 'mask':
        dustmapfn='SFD_mask_4096_%s.fits'
    else:
        dustmapfn=dustmap
    
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
    
    
    
    if '%s' not in dustmapfn:
        f=open(dustmapfn)
        try:
            mapds=[f[0].data]
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        
        polename=dustmapfn.split('.')[0].split('_')[-1].lower()
        if polename=='ngp':
            n=[1]
            if sum(b > 0) > 0:
                print 'using ngp file when lat < 0 present... put %s wherever "ngp" or "sgp" should go in filename'
        elif polename=='sgp':
            n=[-1]
            if sum(b < 0) > 0:
                print 'using sgp file when lat > 0 present... put %s wherever "ngp" or "sgp" should go in filename'
        else:
            raise ValueError("couldn't determine South/North from filename - should have 'sgp' or 'ngp in it somewhere")
        masks = [ones_like(b).astype(bool)]
    else: #need to do things seperately for north and south files
        nmask = b >= 0
        smask = ~nmask
        
        masks = [nmask,smask]
        ns = [1,-1]
        
        mapds=[]
        f=open(dustmapfn%'ngp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        f=open(dustmapfn%'sgp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
    
    retvals=[]
    for n,mapd,m in zip(ns,mapds,masks):
        #project from galactic longitude/latitude to lambert pixels (see SFD98)
        npix=mapd.shape[0]
        
        x=npix/2*cos(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        y=-npix/2*n*sin(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        #now remap indecies - numpy arrays have y and x convention switched from SFD98 appendix
        x,y=y,x
        
        if interpolate:
            from scipy.ndimage import map_coordinates
            if type(interpolate) is int:
                retvals.append(map_coordinates(mapd,[x,y],order=interpolate))
            else:
                retvals.append(map_coordinates(mapd,[x,y]))
        else:
            x=round(x).astype(int)
            y=round(y).astype(int)
            retvals.append(mapd[x,y])
            
            
    
        
    if isscalar(long) or isscalar(lat):
        for r in retvals:
            if len(r)>0:
                return r[0]
        assert False,'None of the return value arrays were populated - incorrect inputs?'
    else:
        #now recombine the possibly two arrays from above into one that looks like  the original
        retval=ndarray(l.shape)
        for m,val in zip(masks,retvals):
            retval[m] = val
        return retval
        
    
def get_dust_radec(ra,dec,dustmap,interpolate=True):
    from .coords import equatorial_to_galactic_flt
    l,b = equatorial_to_galactic_flt(ra,dec)
    return get_dust_gal(l,b,dustmap,interpolate)



#<---------------Deprecated---------------------------------------------------->
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
    
    from warnings import warn
    warn('extinction_correction function is deprecated - use Extinct class instead',DeprecationWarnings)
    
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
    
    from warnings import warn
    warn('extinction_from_flux_ratio function is deprecated - use Extinct class instead',DeprecationWarnings)
    
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
    
