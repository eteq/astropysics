#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module stores tools for oberving (pre- and post-) as well as functioning as
a "miscellaneous" bin for various corrections and calculations that don't have
a better place to live.

The focus is currently on optical astronomy, as that is what the primary author 
does.
"""

from __future__ import division,with_statement
from math import pi
import numpy as np

from .utils import PipelineElement

#<----------------------Extinction and dust-related---------------------------->

class Extinction(PipelineElement):
    """
    This is the base class for extinction-law objects. Extinction laws can be
    passed in as a function to the initializer, or subclasses should override 
    the function f with the preferred law as f(self,lambda), and their 
    initializers should set the zero point
    
    Note that functions are interpreted as magnitude extinction laws ... if 
    optical depth is desired, f should return (-2.5/log(10))*tau(lambda)
    
    A0 is the normalization factor that gets multiplied into the reddening law.
    """
    def  __init__(self,f=None,A0=1):
        if f is not None:
            if not callable(f):
                self.f = f
            else:
                raise ValueError('function must be a callable')
        self.A0 = A0
        
        self._plbuffer = None
    
    def f(self,lamb):
        raise NotImplementedError('must specify an extinction function ')
    
    def __call__(self,*args,**kwargs):
        return self.A0*self.f(*args,**kwargs)
    
    def correctPhotometry(self,mags,band):
        """
        Uses the extinction law to correct a magnitude (or array of magnitudes)
        
        bands is either a string specifying the band, or a wavelength to use
        """
        from .phot import bandwl
        #TODO:better band support
        if isinstance(band,basestring):
            wl = bandwl[band]
        else:
            wl = band
        
        return mags-self(wl)
        
    def Alambda(self,band):
        """
        determines the extinction for this extinction law in a given band or bands
        
        band can be a wavelength or a string specifying a band 
        """
        from operator import isSequenceType
        if isSequenceType(band) and not isinstance(band,basestring):
            return -1*np.array([self.correctPhotometry(0,b) for b in band])
        else:
            return -1*self.correctPhotometry(0,band)
        
    def correctColor(self,colors,bands):
        """
        Uses the supplied extinction law to correct a color (or array of colors)
        where the color is in the specified bands
        
        bands is a length-2 sequence with either a band name or a wavelength for
        the band, or of the form 'bandname1-bandname2' or 'E(band1-band2)'
        """
        
        if isinstance(bands,basestring):
            b1,b2=bands.replace('E(','').replace(')','').split(',')
        else:
            b1,b2=bands
            
        from .phot import bandwl
        #TODO:better band support
        wl1 = bandwl[b1] if isinstance(b1,basestring) else b1
        wl2 = bandwl[b2] if isinstance(b1,basestring) else b2
    
        return colors-self(wl1)+self(wl2)
    
    def correctSpectrum(self,spec,newspec=True):
        """
        Uses the supplied extinction law to correct a spectrum for extinction.
        
        if newspec is True, a copy of the supplied spectrum will have the 
        extinction correction applied
        
        returns the corrected spectrum
        """
    
        if newspec:
            spec = spec.copy()
            
        oldunit = spec.unit
        spec.unit = 'wavelength-angstrom'
        corr = 10**(self(spec.x)/2.5)
        spec.flux *= corr
        spec.err *= corr
        
        spec.unit = oldunit
        return spec
    
    __builtinlines ={
            'Ha':6562.82,
            'Hb':4861.33,
            'Hg':4340.46,
            'Hd':4101.74,
            'He':3970.07
            }
    def correctFlux(self,flux,lamb):
        if isinstance(lamb,basestring) and lamb in Extinction.__builtinlines:
            lamb = Extinction.__builtinlines[lamb]
        return flux*10**(self(lamb)/2.5)
    
    
    __balmerratios={
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
    def computeA0FromFluxRatio(self,measured,expected,lambda1=None,lambda2=None,filterfunc=None):
        """
        This derives the normalization of the Extinction function from provided
        ratios for theoretically expected fluxes.  If multiple measurements are
        provided, the mean will be used
        
        measured is the measured line ratio (possibly an array), while expected
        is either the expected line ratios, or a string specifying the
        appropriate balmer flux ratio as "Hab","Hde",etc. (for Halpha/Hbeta or
        Hdelta/Hepsilon) to assume case B recombination fluxes (see e.g. 
        Osterbrock 2007).
        
        lambda1 and lambda2 are the wavelengths for the ratios F1/F2, or None if
        a string is provided 
    
        filterfunc is a function to be applied as the last step - it can 
        either be used to contract an array to (e.g. np.mean), or filter out
        invalid values (e.g. lambda x:x[np.isfinite(x)]).  Default
        does nothing
        
        returns A0,standard deviation of measurements
        """
        if isinstance(expected,basestring):
            R = measured
            balmertuple = Extinction.__balmerratios[expected]
            R0=balmertuple[0]
            lambda1,lambda2=balmertuple[1],balmertuple[2]
        elif isinstance(expected[0],basestring):
            if not np.all([isinstance(e,basestring) for e in expected]):
                raise ValueError('expected must be only transitions or only numerical')
            R0,lambda1,lambda2=[],[],[]
            for e in expected:
                balmertuple = Extinction.__balmerratios[e]
                R0.append(balmertuple[0])
                lambda1.append(balmertuple[1])
                lambda2.append(balmertuple[2])
            R0,lambda1,lambda2=np.array(R0),np.array(lambda1),np.array(lambda2)    
        else:
            if lambda1 and lambda2:
                R,R0 = np.array(measured,copy=False),np.array(expected,copy=False)
                lambda1,lambda2 = np.array(lambda1,copy=False),np.array(lambda2,copy=False)
            else:
                raise ValueError('need to provide wavelengths if not specifying transitions')
        
        A0 = -2.5*np.log10(R/R0)/(self.f(lambda1)-self.f(lambda2))
        
        if filterfunc is None:
            filterfunc = lambda x:x
        self.A0 = A0 = filterfunc(A0)
        return A0,np.std(A0)
    
    #PipelineElement methods
    def _plFeed(self,data,src):
        from .utils import  PipelineError
        from .spec import Spectrum
        if self._plbuffer is None:
            self._plbuffer = {'in':[],'out':[]} 
    
        if isinstance(data,Spectrum):
            self._plbuffer['in'].append(('spec',data))
        else:
            raise PipelineError('unrecognized Extinction correction input data')
        
    def _plProcess(self):
        from .utils import  PipelineError
        if self._plbuffer is None:
            self._plbuffer = {'in':[],'out':[]} 
        
        type,data = self._plbuffer['in'].pop(0)
        try:
            if type=='spec':
                newspec = self.correctSpectrum(data)
                self._plbuffer['out'].append(newspec)
            else:
                assert False,'Impossible point - code error in Extinction pipeline'
        except:
            self.insert(0,spec(type,data))
            
    def _plExtract(self):
        if self._plbuffer is None:
            self._plbuffer = {'in':[],'out':[]} 
            
        if len(self._plbuffer['out']<1):
            return None
        else:
            return self._plbuffer['out'].pop(0)
        
    def _plClear(self):
        self._plbuffer = None
    
class CalzettiExtinction(Extinction):
    """
    This is the average extinction law derived in Calzetti et al. 1994:
    x=1/lambda in mu^-1
    Q(x)=-2.156+1.509*x-0.198*x**2+0.011*x**3
    """
    _poly=np.poly1d((0.011,-0.198,1.509,-2.156)) 
    def __init__(self,A0=1):
        super(CalzettiExtinction,self).__init__(A0=A0)
        
    def f(self,lamb):
        return (-2.5/np.log(10))*self._poly(1e4/lamb)
    
class _EBmVExtinction(Extinction):
    """
    Base class for Extinction classes that get normalization from E(B-V)
    """
    from .phot import bandwl
    __lambdaV=bandwl['V']
    del bandwl
    
    def __init__(self,EBmV=1,Rv=3.1):
        super(_EBmVExtinction,self).__init__(f=None,A0=1)
        self.Rv=Rv
        self.EBmV=EBmV
    
    def _getEBmV(self):
        return self.A0*self.f(self.__lambdaV)/self.Rv
    def _setEBmV(self,val):
        Av=self.Rv*val
        self.A0=Av/self.f(self.__lambdaV)
    EBmV = property(_getEBmV,_setEBmV)
    
    def f(self,lamb):
        raise NotImplementedError
    
class FMExtinction(_EBmVExtinction):
    """
    Base class for Extinction classes that use the form from 
    Fitzpatrick & Massa 90
    """
    
    def __init__(self,C1,C2,C3,C4,x0,gamma,EBmV=1,Rv=3.1):
        self.x0 = x0
        self.gamma = gamma
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        
        super(FMExtinction,self).__init__(EBmV=EBmV,Rv=Rv)
        
    def f(self,lamb):
        x=1e4/np.array(lamb,copy=False)
        C1,C2,C3,C4 = self.C1,self.C2,self.C3,self.C4
        gamma,x0 = self.gamma,self.x0
        
        xsq=x*x
        D=xsq*((xsq-x0*x0)**2+xsq*gamma*gamma)**-2
        FMf=C1+C2*x+C3*D
        
        if np.isscalar(FMf):
            if x>=5.9:
                FMf+=C4*(0.5392*(x-5.9)**2+0.05644*(x-5.9)**3)
        else:
            C4m=x>=5.9
            FMf[C4m]+=C4*(0.5392*(x[C4m]-5.9)**2+0.05644*(x[C4m]-5.9)**3)
        
        return FMf+self.Rv #EBmV is the normalization and is multiplied in at the end
    
class CardelliExtinction(_EBmVExtinction):
    """
    Milky Way Extinction law from Cardelli et al. 1989
    """
    def f(self,lamb):
        scalar=np.isscalar(lamb)
        x=1e4/np.atleast_1d(lamb) #CCM x is 1/microns
        a,b=np.ndarray(x.shape,x.dtype),np.ndarray(x.shape,x.dtype)
        
        if any((x<0.3)|(10<x)):
            raise ValueError('some wavelengths outside CCM 89 extinction curve range')
        
        irs=(0.3 <= x) & (x <= 1.1)
        opts = (1.1 <= x) & (x <= 3.3)
        nuv1s = (3.3 <= x) & (x <= 5.9)
        nuv2s = (5.9 <= x) & (x <= 8)
        fuvs = (8 <= x) & (x <= 10)
        
        #TODO:pre-compute polys
        
        #CCM Infrared
        a[irs]=.574*x[irs]**1.61
        b[irs]=-0.527*x[irs]**1.61
        
        #CCM NIR/optical
        a[opts]=np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
        b[opts]=np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)
        
        #CCM NUV
        y=x[nuv1s]-5.9
        Fa=-.04473*y**2-.009779*y**3
        Fb=-.2130*y**2-.1207*y**3
        a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)+Fa
        b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)+Fb
        
        a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)
        b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)
        
        #CCM FUV
        a[fuvs]=np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
        b[fuvs]=np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)
        
        AloAv = a+b/self.Rv
        
        if scalar:
            return AloAv[0]
        else:
            return AloAv
    
class LMCExtinction(FMExtinction):
    """
    LMC Extinction law from Gordon et al. 2003 LMC Average Sample
    """
    def __init__(self,EBmV=.3,Rv=3.41):
        super(LMCExtinction,self).__init__(-.890,0.998,2.719,0.400,4.579,0.934,EBmV,Rv)
    
class SMCExtinction(FMExtinction):
    """
    SMC Extinction law from Gordon et al. 2003 SMC Bar Sample
    """
    def __init__(self,EBmV=.2,Rv=2.74):
        super(SMCExtinction,self).__init__(-4.959,2.264,0.389,0.461,4.6,1,EBmV,Rv)

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
    (in which case the files are assumed to lie in the current directory)
    
    input coordinates are in degrees of galactic latiude and logitude - they can
    be scalars or arrays
    
    if interpolate is an integer, it can be used to specify the order of the
    interpolating polynomial
    """
    from numpy import sin,cos,round,isscalar,array,ndarray,ones_like
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
    from .coords import equatorial_to_galactic
    l,b = equatorial_to_galactic(ra,dec)
    return get_SFD_dust(l,b,dustmap,interpolate)


  
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
    warn('extinction_correction function is deprecated - use Extinction class instead',DeprecationWarning)
    
    if exttype=='LMC':
        eo = LMCExtinction(EBmV=EBmV,Rv=Rv)
        return eo.correctFlux(lineflux,linewl)
    elif exttype == 'SMC':
        eo = SMCExtinction(EBmV=EBmV,Rv=Rv)
        return eo.correctFlux(lineflux,linewl)
    
    if isinstance(linewl,basestring):
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
    warn('extinction_from_flux_ratio function is deprecated - use Extinct class instead',DeprecationWarning)
    
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
    if isinstance(frexpect,basestring):
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
        if isinstance(outlambda,basestring):
            from .phot import bandwl
            outlambda=bandwl[outlambda]
        return extinction_correction(None,outlambda,EBmV,Rv)
    else:
        return EBmV
    
