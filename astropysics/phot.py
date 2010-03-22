#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 
"""

====
phot
====

The :mod:`phot` module contains classes and funtions for photometry or other
image-based flux measurements.


.. todo:: examples/tutorials


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.phot
   :parts: 1

Module API
----------

"""

#TODO: c-based kcorrect
#TODO: reduction framework that outputs objcat catalogs and properly uses CCD noise models
#TODO: atmospheric extinction corrections
#TODO: flux calibration

from __future__ import division,with_statement
from math import pi
import numpy as np

from .spec import HasSpecUnits as _HasSpecUnits
from .utils import DataObjectRegistry
try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    
    
#<-----------------------Magnitude systems------------------------------------->

class Magnitude(object):
    __metaclass__ = ABCMeta
    @abstractmethod
    def magToFlux(self,mag):
        raise NotImplementedError
    @abstractmethod
    def magerrToFluxerr(self,err,mag):
        raise NotImplementedError
    @abstractmethod
    def fluxToMag(self,flux):
        raise NotImplementedError
    @abstractmethod
    def fluxerrToMagerr(self,err,flux):
        raise NotImplementedError
    
class PogsonMagnitude(Magnitude):
    a = -2.5/np.log(10.0)
    
    def magToFlux(self,mag):
        return np.exp(mag/self.a)
    def magerrToFluxerr(self,err,mag):
        return -err*self.magToFlux(mag)/self.a
    def fluxToMag(self,flux):
        return self.a*np.log(flux)
    def fluxerrToMagerr(self,err,flux):
        return -self.a*err/flux

class AsinhMagnitude(Magnitude):
    a = -2.5/np.log(10.0)
    
    def __init__(self,b=1e-10):
        self.b = b
        
    def _getb(self):
        return self._b
    def _setb(self,val):
        self._b = val
        self._logb = np.log(b)
    b= property(_getb,_setb)
    
    def magToFlux(self,mag):
        return 2*b*np.sinh(mag/self.a-self._logb)
    def magerrToFluxerr(self,err,mag):
        b = self.b
        #TODO:fix/test
        return 2*b*err/-self.a*(1+self.magToFlux(mag)**2)**0.5
    def fluxToMag(self,flux):
        b = self.b
        return self.a*(np.arcsinh(flux/2/b)+self._logb)
    def fluxerr_to_magerr(self,err,flux):
        #TODO:fix/test
        return -self.a*err/2/b/(1 + (flux/2/b)**2)**0.5


def choose_magnitude_system(system):
    """
    This function is used to change the magnitude system used where magnitudes
    are used in astropysics. It can be:
    
    * 'pogson'
        :math:`M = -2.5 \\log_{10}(f/f_0)`
    * 'asinh'
        :math:`M = -2.5 \\log_10(e) [{\\rm asinh}(x/2b)+\\ln(b)]` 
    """
    #TODO: choose b from the band or something instead of fixed @ 25 zptmag
    global _magsys,_mag_to_flux,_magerr_to_fluxerr,_flux_to_mag,_fluxerr_to_magerr
    if isinstance(system,Magnitude):
        _magsys = system
    elif system == 'pogson':
        _magsys = PogsonMagnitude()
    elif system == 'asinh':
        _magsys = AsinhMagnitude(b=1e-10)
    else:
        raise ValueError('unrecognized ')
    
    _mag_to_flux = _magsys.magToFlux
    _magerr_to_fluxerr = _magsys.magerrToFluxerr
    _flux_to_mag = _magsys.fluxToMag
    _fluxerr_to_magerr = _magsys.fluxerrToMagerr   
    
choose_magnitude_system('pogson')

#<---------------------Photometric Classes------------------------------------->

class Band(_HasSpecUnits):
    """
    This class is the base of all photometric band objects.  Bands are expected
    to be immutable once created except for changes of units.
    
    The response functions should be photon-counting response functions (see
    e.g. Bessell et al 2000).
    
    Subclasses must implement the following methods:
    
    * :meth:`__init__`
        initialize the object
    * :meth:`_getCen`
        return the center of band 
    * :meth:`_getFWHM`
        return the FWHM of the band
    * :meth:`_getx`
        return an array with the x-axis
    * :meth:`_getS`
        return the photon sensitivity (shape should match x)
    * :meth:`_applyUnits`
        units support (see :class:`astropysics.spec.HasSpecUnits`)
    
    the following can be optionally overriden:
    
    * :meth:`alignBand`
        return the sensitivity of the Band interpolated to the provided x-array
    * :meth:`alignToBand`: 
        return the provided sensitivity interpolated onto the x-array of the
        Band

    """
    __metaclass__ = ABCMeta
    
    #set default units to angstroms
    _phystype = 'wavelength'
    _unit = 'angstrom'
    _xscaling = 1e-8
    
    @abstractmethod
    def __init__(self):
        raise NotImplementedError
    
    @abstractmethod
    def _getCen(self):
        raise NotImplementedError
    
    @abstractmethod
    def _getFWHM(self):
        raise NotImplementedError
    
    @abstractmethod
    def _getx(self):
        raise NotImplementedError
    
    @abstractmethod
    def _getS(self):
        raise NotImplementedError
    
    #comparison operators TODO:unit matching
    def __eq__(self,other):
        try:
            return np.all(self.x == other.x) and np.all(self.S == other.S)
        except:
            return False
    def __ne__(self,other):
        return not self.__eq__(other)
    def __lt__(self,other):
        oldu = None
        try:
            oldu = other.unit
            other.unit = self.unit
            return self.cen < other.cen
        except:
            return NotImplemented
        finally:
            if oldu is not None:
                other.unit = oldu 
    def __le__(self,other):
        oldu = None
        try:
            oldu = other.unit
            other.unit = self.unit
            return self.cen <= other.cen
        except:
            return NotImplemented
        finally:
            if oldu is not None:
                other.unit = oldu 
    def __gt__(self,other):
        oldu = None
        try:
            oldu = other.unit
            other.unit = self.unit
            return self.cen > other.cen
        except:
            return NotImplemented
        finally:
            if oldu is not None:
                other.unit = oldu 
    def __ge__(self,other):
        oldu = None
        try:
            oldu = other.unit
            other.unit = self.unit
            return self.cen >= other.cen
        except:
            return NotImplemented
        finally:
            if oldu is not None:
                other.unit = oldu 
    
    def __interp(self,x,x0,y0,interpolation):
        """
        interpolate the function defined by x0 and y0 onto the coordinates in x
        
        see Band.alignBand for interpolation types/options
        """
        if interpolation == 'linear':
            return np.interp(x,x0,y0)
        elif 'spline' in interpolation:
            if not hasattr(self,'_splineobj'):
                self._splineobj = None
                
            deg = interpolation.replace('spline','').strip()
            deg = 3 if deg == '' else int(deg)
            if self._splineobj is None or self._splineobj._data[5] != deg:
                from scipy.interp import InterpolatedUnivariateSpline
                self._splineobj = InterpolatedUnivariateSpline(k=deg)
            return self._splineobj(x)
                
        else:
            raise ValueError('unrecognized interpolation type')
    
    
    cen = property(lambda self:self._getCen())
    FWHM = property(lambda self:self._getFWHM())
    x = property(lambda self:self._getx())
    S = property(lambda self:self._getS())
    def _getSrc(self):
        if not hasattr(self._src):
            from .objcat import Source
            self._src = Source(None)
        return self._src
    def _setSrc(self,val=None):
        pass
        #from .objcat import Source
        #self._src = Source(val)
    source = property(_getSrc,_setSrc)
    def _getName(self):
        if self._name is None:
            return 'Band@%4.4g'%self.cen
        else:
            return self._name
    def _setName(self,val):
        self._name = val
    name = property(_getName,_setName)
    _name = None #default is nameless
    
    zptflux = 1 #flux at mag=0
    
    def alignBand(self,x,interpolation='linear'):
        """
        interpolates the Band sensitivity onto the provided x-axis
        
        x can either be an array with the wanted x-axis or a ``Spectrum``
        if it is a spectrum, the units for the band will be converted to 
        the Spectrum units before alignment
        
        interpolation can be 'linear' or 'spline'
        """
        if hasattr(x,'x') and hasattr(x,'unit'):
            units = x.unit
            x = x.x
        oldunit = self.unit
        try:
            self.unit = units
            sorti = np.argsort(self.x)
            return self.__interp(x,self.x[sorti],self.S[sorti],interpolation)
        finally:
            self.unit = oldunit
        
    def alignToBand(self,*args,**kwargs):
        """
        Interpolates the provided function onto the :class:`Band` x-axis. Input
        arguments can be of the following forms:
        
        * alignToBand(arr)
            assume that the provided sequence fills the band and interpolate it
            onto the :class:`Band` coordinates
        * alignToBand(specobj)
            if specobj is of type :class:`astropyics.spec.Spectrum`, the
            spectrum's flux will be returned aligned to the band, in the units
            of the :class:`Band`
        * alignToBand(x,y)
            the function with coordinates (x,y) will be aligned to the
            :class:`Band` x-axis
        
        *keywords*
            * `interpolation`
                see :meth:`alignBand` for interpolation types/options
        """
        if len(args) == 1:
            from operator import isSequenceType
            from .spec import Spectrum,FunctionSpectrum
            specin = args[0]
            
            if isSequenceType(specin):
                y = np.array(specin,copy=False)
                x = np.linspace(self.x[0],self.x[-1])
            elif isinstance(specin,FunctionSpectrum):
                oldx = specin.x
                try:
                    specin.x = self.x
                    return specin.getUnitFlux(self.unit)
                finally:
                    specin.x = oldx
            elif hasattr(specin,'x') and hasattr(specin,'flux'):
                x,y = specin.getUnitFlux(self.unit)
                #x = specin.x
                #y = specin.flux
        elif len(args) == 2:
            x,y = args
        else:
            raise TypeError('alignToBand() takes 1 or 2 arguments')
        
        if 'interpolation' not in kwargs:
            kwargs['interpolation'] = 'linear'
        
        sorti = np.argsort(x)    
        return self.__interp(self.x,x[sorti],y[sorti],kwargs['interpolation'])
    
    def isOverlapped(self,xorspec):
        """
        determes if this Band overlaps on a significant part of the provided 
        spectrum (that is, if center+/-fwhm overlaps) 
        """
        
        if hasattr(xorspec,'x') and hasattr(xorspec,'flux'):
            x = xorspec.x
        else:
            x = xorspec
            
        mx,mn = np.max(x),np.min(x)
        c,w = self.cen,self.FWHM/2
        return (mx >= c+w) and (mn <= c-w)
    
    def computeMag(self,*args,**kwargs):
        """
        Compute the magnitude of the supplied Spectrum in this band using the 
        band's ``zptflux`` attribute
        
        args and kwargs go into Band.computeFlux
        """
        return _flux_to_mag(self.computeFlux(*args,**kwargs)/self.zptflux)
    
    def computeFlux(self,spec,interpolation='linear',aligntoband=None,overlapcheck=True):
        """
        compute the flux in this band for the supplied Spectrum
        
        the Spectrum is aligned to the band if aligntoband is True, or vice 
        versa, otherwise, using the specified interpolation (see Band.alignBand
        for interpolation options).  If aligntoband is None, the higher of the
        two resolutions will be left unchanged 
        
        if overlapcheck is True, a ValueError will be raised if most of the 
        Spectrum does not lie within the band
        
        the spectrum units will be converted to match those of the Band
        
        the spectrum can be an array, but then aligntoband and 
        interpolation are ignored and must match the band's x-axis
        """
        from operator import isSequenceType
        from scipy.integrate import simps as integralfunc
        if isSequenceType(spec):
            x = self.x
            y = np.array(spec)
            if y.shape != x.shape:
                raise ValueError('input array shape does not match Band x-axis')
            units = self.unit
        elif hasattr(spec,'x') and hasattr(spec,'flux'):
            oldunits = spec.unit
            try:
                spec.unit = self.unit
                if aligntoband is None:
                    lx,px = self.x,spec.x
                    aligntoband = lx.size/(lx.max()-lx.min()) > px.size/(px.max()-px.min())
                
                if aligntoband:
                    x = self.x
                    y = self.S*self.alignToBand(spec.x,spec.flux,interpolation=interpolation)
                else:
                    x = spec.x.copy()
                    y = self.alignBand(spec)*spec.flux
            finally:
                spec.unit=oldunits
        else:
            raise ValueError('unrecognized input spectrum')
        
        if overlapcheck and not self.isOverlapped(x):
            raise ValueError('provided input does not overlap on this band')
            
        if 'wavelength' in self.unit:
            y*=x
        else:
            y/=x
            #TODO:check energy factor unit-wise
            
        sorti=np.argsort(x)
        return integralfunc(y[sorti],x[sorti])
    
    def computeZptFromSpectrum(self,*args,**kwargs):
        """
        use the supplied Spectrum as a zero-point Spectrum to 
        define the zero point of this band as well as the flux 
        for that zero-point
        
        args and kwargs are the same as those for computeMag
        """
        #fluxfactor = self.computeFlux(np.ones_like(self.x)) #scale by this factor to get actual flux
        #mag = self.computeMag(*args,**kwargs)
        self.zptflux = self.computeFlux(*args,**kwargs)
        
    
    def plot(self,spec=None,bandscaling=1,bandoffset=0,labelband=True,clf=True,**kwargs):
        """
        this plots the requested band and passes kwargs into matplotlib.pyplot.plot
        OR
        if spec is provided, the spectrum will be plotted along with the band and 
        the convolution
        
        *keywords*
            
            * `clf`
                clear the figure before plotting
            * `leg` 
                show legend where appropriate
        
        """
        from matplotlib import pyplot as plt
        band = self
        
        if clf:
            plt.clf()
                
        if spec:
            from .spec import Spectrum
            if not isinstance(spec,Spectrum):
                raise ValueError('input not a Spectrum')
            bandscaling = np.max(spec.flux)*bandscaling
            convflux = band.alignToBand(spec)*band.S
            
            plt.plot(band.x,bandscaling*band.S+bandoffset,c='k',label='$S_{\\rm Band}$')
            bandmax = bandscaling*band.S.max()+bandoffset
            plt.plot(spec.x,spec.flux,c='g',ls=':',label='$f_{\\rm spec}$')
            plt.plot(band.x,convflux,c='b',ls='--',label='$f_{\\rm spec}S_{\\rm Band}$')
            
            plt.ylabel('Spec flux or S*max(spec)')
            if kwargs.pop('leg',True):
                plt.legend(loc=0)
            mi,ma=band.x.min(),band.x.max()
            rng = ma-mi
            plt.xlim(mi-rng/3,ma+rng/3)
        else:
            bandmax = bandscaling*band.S.max()+bandoffset
            plt.plot(band.x,bandscaling*band.S+bandoffset,**kwargs)
            plt.ylabel('$S$')
            plt.ylim(0,bandscaling*band.S.max())
            
        if labelband and self.name is not None:
            plt.text(band.cen,bandmax,self.name)
        plt.xlabel('$\\lambda$')
            

class GaussianBand(Band):
    def __init__(self,center,width,A=1,sigs=6,n=100,unit='angstroms',name=None):
        """
        center is the central wavelength of the band, while width is the sigma 
        (if positive) or FWHM (if negative)
        """
        self._cen = center
        self._sigma = width if width > 0 else -width*(8*np.log(2))**-0.5
        self._fwhm = self._sigma*(8*np.log(2))**0.5
        self._A = A
        
        self._n = n
        self._sigs = sigs
        self._updateXY(self._cen,self._sigma,self._A,self._n,self._sigs)
        
        _HasSpecUnits.__init__(self,unit)
        self._unittrans = (lambda t:t,lambda t:t)
        
        self.name = name
        
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        self._unittrans = (xtrans,xftrans,xitrans)
        self._updateXY(self._cen,self._sigma,self._A,self._n,self._sigs)
        
    def _updateXY(self,mu,sigma,A,n,sigs):
        xtrans,xftrans,xitrans = self._unittrans
        
        x = np.linspace(-sigs*sigma,sigs*sigma,n)+mu
        xp = (self._x-mu)/sigma
        y = A*np.exp(-xp*xp/2)
        self._x,self._y = xftrans(x,y)
        
    def _getSigs(self):
        return self._sigs
    def _setSigs(self,val):
        self._updateXY(self._cen,self._sigma,self._A,self._n,val)
        self._sigs = val
    sigs = property(_getSigs,_setSigs)
    
    def _getn(self):
        return self._n
    def _setn(self,val):
        self._updateXY(self._cen,self._sigma,self._A,val,self._sigs)
        self._n = val
    n = property(_getn,_setn)

    @property
    def sigma(self):
        return self.unittrans[0](self._sigma)
    
    def _getCen(self):
        return self.unittrans[0](self._cen)
    
    def _getFWHM(self):
        return self.unittrans[0](self._fwhm) #=sigma*(8*log(2))**0.5
    
    def _getx(self):
        return self._x
    
    def _getS(self):
        return self._y
    
    def alignBand(self,x):
        xtrans,xftrans,xitrans = self._unittrans
        oldx = xitrans(x)
        xp = (x-self._cen)/self._sigma 
        return xftrans(oldx,self._A*np.exp(-xp*xp/2))[1]
        
class ArrayBand(Band):
    def __init__(self,x,S,copyonaccess=True,normalized=True,unit='angstroms',name=None):
        """
        generate a Band from a supplied quantum (e.g. photonic) response 
        function 
        
        if the response function is ever negative, the band will be offset so
        that the lowest point is 0 .  If the band x-axis is not sorted, it 
        will be sorted from lowest to highest.
        """
        self._x = np.array(x,copy=True)
        self._S = np.array(S,copy=True)
        sorti = np.argsort(self._x)
        self._x,self._S = self._x[sorti],self._S[sorti]

        if self._x.shape != self._S.shape or len(self._x.shape) != 1:
            raise ValueError('supplied x and S are not matching 1D arrays')
        if self._S.min() < 0:
            self._S -= self._S.min()
        self._N = self._S.max()
        if normalized:
            self._S /= self._N
        self.copyonaccess = copyonaccess
        self._norm = normalized
        
        self._computeMoments()
        
        _HasSpecUnits.__init__(self,unit)
        
        self.name = name
        
    #units support
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        xfinplace(self._x,self._S) 
        mx = self._S.max()
        self._S/=mx
        if not self._norm:
            self._S/=self._N
        self._computeMoments() #need to recompute moments
        
    def _computeMoments(self):
        from scipy.integrate import  simps as integralfunc
        #trapz vs. simps: trapz faster by factor ~5,simps somewhat more accurate
        x=self._x
        y=self._S/self._S.max()
        N=1/integralfunc(y,x)
        yn = y*N #normalize so that overall integral is 1
        self._cen = integralfunc(x*yn,x)
        xp = x - self._cen
        
        
        if y[0] > 0.5 or y[-1] > 0.5:
            self._fwhm = (integralfunc(xp*xp*yn,xp)*8*np.log(2))**0.5 
        else:
            #from scipy.interpolation import interp1d
            #from scipy.optimize import newton
            #ier = interp1d(x,y-0.5)
            #lval-uval=newton(ier,x[0])-newton(ier,x[-1])
            #need peak to be 1 for this algorithm
            yo = y-0.5
            edgei = np.where(yo*np.roll(yo,1)<0)[0]-1
            li,ui = edgei[0]+np.arange(2),edgei[1]+np.arange(2)
            self._fwhm = np.interp(0,yo[ui],x[ui])-np.interp(0,yo[li],x[li])
            
    
    def _getNorm(self):
        return self._norm
    def _setNorm(self,val):
        val = bool(val)
        if val != self._norm:
            if self._norm:
                self.S *= self._N
            else:
                self.S /= self._N
        self._norm = val
    normalized = property(_getNorm,_setNorm)
        
    def _getCen(self):
        return self._cen
    
    def _getFWHM(self):
        return self._fwhm
    
    def _getx(self):
        if self.copyonaccess:
            return self._x.copy()
        else:
            return self._x
    
    def _getS(self):
        if self.copyonaccess:
            return self._S.copy()
        else:
            return self._S
        
        
class FileBand(ArrayBand):
    def __init__(self,fn,type=None):
        """
        type can be 'txt', or 'fits', or inferred from extension
        
        if txt, first column should be lambda, second should be response
        """
        from os import path
        
        self.fn = fn
        if type is None:
            ext = path.splitext(fn)[-1].lower()
            if ext == 'fits' or ext == 'fit':
                type = 'fits'
            else:
                type = 'txt'
        
        if type == 'txt':
            x,S = np.loadtxt(fn).T
        elif type == 'fits':
            import pyfits
            f = pyfits.open(fn)
            try:
                #TODO: much smarter/WCS
                x,S = f[0].data
            finally:
                f.close()
        else:
            raise ValueError('unrecognized type')
        
        ArrayBand.__init__(self,x,S)
        
        
def plot_band_group(bandgrp,**kwargs):
    """
    Plot a group of bands on the same plot using :mod:`matplotlib`
    
    `group` can be a sequnce of bands, strings or a dictionary mapping names to
    bands
    
    kwargs go into :meth:`matplotlib.pyplot.plot` except for:
    
    * `clf`
        clear figure before plotting
    * `leg`
        show legend with band names
    * `unit`
        unit to use on plots (defaults to angstroms)
    * `c`
        (do not use)
    * `spec`
        (do not use)
    """
    from matplotlib import pyplot as plt
#    from operator import isMappingType,isSequenceType
    
    
#    if isMappingType(bandgrp):
#        names = bandgrp.keys()
#        bs = bandgrp.values()
#    elif isinstance(bandgrp,basestring):
#        names = bands.getGroupDataNames(bandgrp)
#        bs = str_to_bands(bandgrp)
#    elif isSequenceType(bandgrp):
#        names = bandgrp
#        bs = str_to_bands(bandgrp)
#    else:
#        raise TypeError('invalid bandgrp input')
    bs = str_to_bands(bandgrp)
    names = [b.name for b in bs]

    d = dict([(n,b) for n,b in  zip(names,bs)])
    sortednames = sorted(d,key=d.get)
    
    clf = kwargs.pop('clf',True)
    leg = kwargs.pop('leg',True)
    kwargs.pop('spec',None)
    cs = kwargs.pop('c',None)
    unit = kwargs.pop('unit','wl')
       
    oldunits = [b.unit for b in bs]
    try:
        for b in bs:
            b.unit = unit
    except NotImplementedError:
        pass
    
    
    if clf:
        plt.clf()
    
    bluetored = 'wavelength' in bs[0].unit
    
    if cs is None:
        
        x = np.arange(len(bs))
        b = 1-x*2/(len(x)-1)
        g = 1-2*np.abs(0.5-x/(len(x)-1))
        g[g>0.5] = 0.5
        r = (x*2/(len(x)-1)-1)
        r[r<0] = b[b<0] = 0
        
        if not bluetored:
            b,r = r,b
        cs = np.c_[r,g,b]
        
        #custom coloring for UBVRI and ugriz
        if len(sortednames) == 5:
            if sortednames[:3] == ['U','B','V']:
                cs = ['c','b','g','r','m']
            elif np.all([(n in sortednames[i]) for i,n in enumerate(('u','g','r','i','z'))]):
                cs = ['c','g','r','m','k']
        
    
    kwargs['clf'] = False
    for n,c in zip(sortednames,cs):
        kwargs['label'] = n
        kwargs['c'] = c
        d[n].plot(**kwargs)
    
    if leg:
        plt.legend(loc=0)
    
    try:
        for u,b in zip(oldunits,bs):
            b.unit = u
    except (NameError,NotImplementedError),e:
        pass
    
    
def compute_band_zpts(specorscale,bands):
    """
    this will use the supplied spectrum to determine the magnitude zero points
    for the specified bands 
    OR
    scale the zeropoint 
    """
    bnds = str_to_bands(bands)
    
    if np.isscalar(specorscale):
        for b in bnds:
            b.zptflux *= specorscale
    else:
        for b in bnds:
            b.computeZptFromSpectrum(specorscale)
            
def set_zeropoint_system(system,bands='all'):
    """
    this uses a standard system to apply the zeropoints for the specified band
    names.
    
    system can be:
    * 'AB' - use AB magnitudes
    * 'vega' - use vega=0 magnitudes based on Kurucz '93 Alpha Lyrae models
    * 'vega##' - use vega=## magnitudes based on Kurucz '93 Alpha Lyrae models
    """
    
    from .spec import Spectrum,FunctionSpectrum,ABSpec
    bands = str_to_bands(bands)
    
    if system == 'AB':
        s = ABSpec.copy()
        for b in bands:
            s.unit = b.unit
            s.x = b.x
            b.computeZptFromSpectrum(s,aligntoband=True,interpolation='linear')
        
    elif 'vega' in system.lower():
        from cPickle import loads
        from .io import _get_package_data
        
        try:
            offset = float(system.lower().replace('vega',''))
        except ValueError:
            offset = 0

        vegad = loads(_get_package_data('vega_k93.pydict'))['data']
        s = Spectrum(vegad['WAVELENGTH'],vegad['FLUX']*10**(offset/2.5))
        
        for b in bands:
            b.computeZptFromSpectrum(s)
    else:
        raise ValueError('unrecognized magnitude system')
    
    
class PhotObservation(object):
    """
    A photometric measurement (or array of measurements) in a fixed 
    group of bands.  These bands must be present in the band registry.
    
    the first index of the values must have length equal to nbands
    
    if asmags is True, values and errs will be interpreted as magnitudes,
    otherwise, flux
    """
    
    def __init__(self,bands,values,errs=None,asmags=True):
        self._bandnames = tuple([b.name for b in str_to_bands(bands)])
            
        if asmags:
            self.mag = values
            if errs is None:
                errs = np.zeros_like(values) 
            elif np.isscalar(errs):
                errs = errs*np.ones_like(values) 
            self.magerr = errs
        else:
            self.flux = values
            if errs is None:
                errs = np.zeros_like(values) 
            elif np.isscalar(errs):
                errs = errs*np.ones_like(values) 
            self.fluxerr = errs
            
    def __str__(self):
        ss = ['Band %s:%s'%(bn,v) for bn,v in zip(self.bandnames,self._values)]
        ss.append('(mag)' if self._mags else '(flux)')
        return ' '.join(ss)
            
    def __len__(self):
        return len(self._bandnames)
    
    @property
    def bandnames(self):
        return self._bandnames
    
    @property
    def bands(self):
        return str_to_bands(self._bandnames)
    
    def _zeroPoints(self):
        return np.array([b.zptflux for b in self.bands])
    
    def _getMags(self):
        if self._mags:
            return self._values.copy()
        else:  
            zptfluxes = self._zeroPoints()
            zptfluxes = zptfluxes.reshape((zptfluxes.size,1))
            if len(self._values.shape)==1:
                return _flux_to_mags(self._values.reshape(zpts.shape)/zptfluxes).ravel()
                #return (-2.5*np.log10(self._values.reshape(zpts.shape))-zpts).ravel()
            else:
                return _flux_to_mags(self._values)
                #return -2.5*np.log10(self._values)-zpts
    def _setMags(self,val):
        if len(val) != len(self._bandnames):
                raise ValueError('input length does not match number of bands')
        self._mags = True
        self._values = np.array(val)
    mag = property(_getMags,_setMags,doc='photmetric measurements in magnitudes')
    def _getMagsErr(self):
        if self._mags:
            return self._err
        else:
            return _fluxerr_to_magerr(self._err,self._values)
            #return -2.5/np.log(10)*self._err/self._values
    def _setMagsErr(self,val):
        val = np.array(val)
        if val.shape != self._values.shape:
            raise ValueError("Errors don't match values' shape")
        if self._mags:
            self._err = val
        else:
            self._err = _magerr_to_fluxerr(val,self.mag)
            #self._err = self._values*val*np.log(10)/-2.5
    magerr = property(_getMagsErr,_setMagsErr,doc='photmetric errors in magnitudes')
    
    def _getFlux(self):
        if self._mags:
            zptfluxes = self._zeroPoints()
            zptfluxes = zptfluxes.reshape((zptfluxes.size,1))
            if len(self._values.shape)==1:
                return zptfluxes[:,0]*_mag_to_flux(self._values)
                #return (10**((self._values.reshape(zpts.shape)+zpts)/-2.5)).ravel()
            else:
                return zptfluxes*_mag_to_flux(self._values)
                #return 10**((self._values+zpts)/-2.5)
        else:  
            return self._values.copy()
    def _setFlux(self,val):
        if len(val) != len(self._bandnames):
                raise ValueError('input length does not match number of bands')
        self._mags = False
        self._values = np.array(val)
    flux = property(_getFlux,_setFlux,doc='photmetric measurements in flux units')
    def _getFluxErr(self):
        if self._mags:
            zptfluxes = self._zeroPoints()
            zptfluxes = zptfluxes.reshape((zptfluxes.size,1))
            if len(self._values.shape)==1:
                return zptfluxes[:,0]*_magerr_to_fluxerr(self._err,self._values)
            else:
                return zptfluxes*_magerr_to_fluxerr(self._err,self._values)
            #return _magerr_to_fluxerr(self._err,self._values)
            #return self.flux*self._err*np.log(10)/-2.5
        else:
            return self._err.copy()
    def _setFluxErr(self,val):
        val = np.array(val)
        if val.shape != self._values.shape:
            raise ValueError("Errors don't match values' shape")
        if self._mags:
            self._err = _fluxerr_to_magerr(val,self.mag)
            #self._err = -2.5/np.log(10)*val/self.flux
        else:
            self._err = val
    fluxerr = property(_getFluxErr,_setFluxErr,doc='photmetric errors in flux units')
    
    def getFluxDensity(self,unit='angstroms'):
        """
        compute an approximate flux density (e.g. erg cm^-2 s^-1 angstrom^-1)
        assuming flat over the FWHM        
        """
        f,e = self.flux,self.fluxerr
        x,w = self.getBandInfo(unit)
        if len(f.shape) > 1:
            w = np.tile(w,f.size/len(x)).reshape((len(x),f.size/len(x)))
        
        return f/w,e/w
    
    def getBandInfo(self,unit='angstroms'):
        """
        Extracts the center and FWHM of the bands in these observations
        
        returns x,w
        """
        x=[]
        w=[]
        for b in self.bands:
            oldunit = b.unit
            b.unit = unit
            x.append(b.cen)
            w.append(b.FWHM)
            b.unit = oldunit
        return np.array(x),np.array(w)
    
    def toSpectrum(self,unit='angstroms'):
        """
        converts this data set to a Spectrum
        
        TODO: fill in gaps assuming a flat spectrum
        """
        from .spec import Spectrum
        
        y = self.flux
        err = self.fluxerr
        
        x=[]
        w=[]
        for b in self.bands:
            oldunit = b.unit
            b.unit = unit
            x.append(b.cen)
            w.append(b.FWHM)
            utype = b._phystype
            b.unit = oldunit
            
        x = np.tile(x,np.prod(y.shape)/y.shape[0]).reshape((y.shape))
        w = np.tile(w,np.prod(y.shape)/y.shape[0]).reshape((y.shape))
        return Spectrum(x,y/w,err/w,unit=unit)
    
    def plot(self,includebands=True,fluxtype=None,unit='angstroms',clf=True,**kwargs):
        """
        Plots these photometric data points.
        
        `fluxtype` can be:
        
        * None
            either magnitude or flux depending on the input type
        * 'mag'
            magnitudes
        * 'flux'
            flux in :math:`erg cm^{-2} s^{-1}`
        * 'fluxden'
            flux spectral density e.g. :math:`erg cm^{-2} s^{-1} angstroms^{-1}`
        
        kwargs go into matplotlib.pyplot.errorbar
        """
        from matplotlib import pyplot as plt
        
        if fluxtype is None:
            asmags = self._mags
            den = False
        elif fluxtype == 'mag':
            asmags = True
            den = False
        elif fluxtype == 'flux':
            asmags = False
            den = False
        elif fluxtype == 'fluxden':
            asmags = False
            den = True
        else:
            raise ValueError('unrecognized fluxtype')
            
        bands = self.bands
        nb = len(bands) 
        
        if asmags:
            y = self.mag
            yerr = self.magerr
        elif den:
            y,yerr = self.getFluxDensity(unit=unit)
        else:
            y = self.flux
            yerr = self.fluxerr
            
        npts = np.prod(y.shape)/nb
        y,yerr=y.reshape((nb,npts)),yerr.reshape((nb,npts))
        
        x=[]
        #w=[]
        for b in bands:
            oldunit = b.unit
            b.unit = unit
            x.append(b.cen)
            #w.append(b.FWHM)
            utype = b._phystype
            b.unit = oldunit
        x = np.tile(x,npts).reshape((nb,npts))
        #w = np.tile(w,npts).reshape((nb,npts))
            
        if clf:
            plt.clf()
            
        if 'fmt' not in kwargs:
            kwargs['fmt'] = 'o'
        if 'ecolor' not in kwargs:
            kwargs['ecolor'] = 'k'
        plt.errorbar(x.ravel(),y.ravel(),yerr.ravel() if np.any(yerr) else None,**kwargs)
       
        if includebands:
            if includebands is True:
                includebands = {}
            if 'bandscaling' not in includebands:
                includebands['bandscaling'] = 0.5*(y.max()-y.min())
            if 'bandoffset' not in includebands:
                includebands['bandoffset'] = y.min()
            if 'ls' not in includebands:
                includebands['ls'] = '--'
            if 'leg' not in includebands:
                includebands['leg'] = False
            includebands['clf'] = False
            
            xl = plt.xlim()
            plot_band_group(bands,**includebands)
            plt.xlim(*xl)
            
        if utype == 'wavelength':
            plt.xlabel('$\\lambda$')
        elif utype=='frequency':
            plt.xlabel('$\\nu$')
        else:
            plt.xlabel(unit)
        plt.ylabel('Mag' if asmags else 'Flux[${\\rm erg} {\\rm s}^{-1} {\\rm cm}^{-2} {%s}^{-1}$]'%unit)
        rng = y.max() - y.min()
        plt.ylim(y.min()*1.05-y.max()*.05,y.min()*-0.05+y.max()*1.05)

#<---------------------Analysis Classes/Tools---------------------------------->

class CMDAnalyzer(object):
    """
    This class is intended to take multi-band photometry and compare it
    to fiducial models to find priorities/matches for a set of data from the
    same bands
    """
    
    @staticmethod
    def _dataAndBandsToDicts(bandsin,arr):
        from warnings import warn
        from operator import isSequenceType
        
        if len(bandsin) != arr.shape[0]:
            raise ValueError('bandsin does not match fiducial bands')
        
        n = arr.shape[1]
        
        bandnames,bandobjs,colors=[],[],[]
        for b in bandsin:
            if type(b) is str:
                if '-' in b:
                    bandnames.append(None)
                    bandobjs.append(None)
                    colors.append(b)
                else:
                    bandnames.append(b)
                    bandobjs.append(bands[b])
                    colors.append(None)
            elif hasattr(b,'name'):
                bandnames.append(b.name)
                bandobjs.append(b)
                colors.append(None)
            else:
                raise ValueError("Can't understand input band "+str(b))
            
        bd={}
        datad={}
        maskd = {}
        for i,(n,b) in enumerate(zip(bandnames,bandobjs)):
            if b is not None:
                if n is None:
                    bd[b.name] = b
                    datad[b.name] = arr[i]
                    maskd[b.name] = True
                else:
                    bd[n] = b
                    datad[n] = arr[i]
                    maskd[n] = True
                    
        for i,c in enumerate(colors):
            if c is not None:
                b1,b2 = (b.strip() for b in c.split('-'))
                if b1 in bd and not b2 in bd:
                    bd[b2] = bands[b2]
                    datad[b2] = datad[b1] - arr[i]
                    maskd[b2] = True
                elif b2 in bd and not b1 in bd:
                    bd[b1] = bands[b1]
                    datad[b1] = arr[i] + datad[b2]
                    maskd[b1] = True
                elif b1 in bd and b2 in bd:
                    warn('Bands %s and %s overedetermined due to color %s - using direct band measurements'%(b1,b2,c))
                    #do nothing because bands are already defined
                else:
                    warn('Bands %s and %s underdetermined with only color %s - assuming zeros for band %s'%(b1,b2,c,b1))
                    bd[b1] = bands[b1]
                    datad[b1] = np.zeros(len(arr[i]))
                    maskd[b1] = False
                    bd[b2] = bands[b2]
                    datad[b2] = -1*arr[i]
                    maskd[b2] = True
                
        return bd,datad,maskd
    
    def __init__(self,fiducial,fbands,fidname='fiducial'):
        """
        fiducial is a B x N array where B is the number of bands/colors
        and N is the number of fiducial points.  
        
        fbands specifies the band (either as a string or a Band object) or 
        color name (seperated by "-") matching the fiducial
        
        fidname is a string or a dictionary that maps fiducial names to 
        indecies in the fidcuial array.  If a dictionary, each index must 
        be in one of the values (but not in more than one)
        """
        from warnings import warn
        from operator import isSequenceType
        
        arr = np.array(fiducial,ndmin=2)
        bandd,datad,maskd = self._dataAndBandsToDicts(fbands,arr)
                
        self._banddict = bandd
        self._fdatadict = datad 
        self._dm0fdd = datad #distance modulus 0 data
        self._dmod = 0
        self._bandnames = tuple(sorted(bandd.keys(),key=bandd.get))
        self._fidmask = np.array([maskd[b] for b in self._bandnames])
        self._fmaskd = maskd
        self._nb = len(self._bandnames)
        self._nfid = arr.shape[1]
        
        if isinstance(fidname,basestring):
            self._fidnamedict={fidname:np.arange(self._nfid)}
        else:
            d,vals = {},[]
            vals = []
            for k,v in fidname.iteritems():
                d[k] =  np.array(v,dtype=int)
                vals.extend(d[k])
            if not np.all(np.sort(vals) == np.arange(self._nfid)):
                raise ValueError('fidname dictionary does not account for all fiducial array elements or duplicate names')
            self._fidnamedict = d
        
        self._data = None
        self._datamask = np.zeros(self._nb,dtype=bool)
        self._nd = 0
        self._cen = (0,0,0)
        self._locs = None
        self._dnames = None
        self._dgrp = None
        
        self._offsets  = None
        
        self._offbands = None
        self._offws = None
        self._locw = 1
        
        
        
    def _calculateOffsets(self):
        #use self._offbands and self._dmod to set self._offsets
        if self._offbands is None:
            bands = [self._bandnames[i] for i,b in enumerate(self._fidmask & self._datamask) if b]
            fids = np.array([self._fdatadict[b]+self._dmod for b in bands],copy=False).T
            #fids = np.array([self.getFiducialData(b) for b in bands],copy=False).T
            lbds = list(self._bandnames)
            data = np.array([self._data[lbds.index(b)] for b in bands],copy=False).T
        else:
            lbns=list(self._bandnames)
            
            #first check to make sure offset bands are valid
            for b in self._offbands:
                if isinstance(b,tuple):
                    i1 = lbns.index(b[0])
                    i2 = lbns.index(b[1])
                    assert (isinstance(i1,int) and i1 >=0 and i1 < self._nb),(i1,b1)
                    assert (isinstance(i2,int) and i2 >=0 and i2 < self._nb),(i2,b2)
                    if not ((self._fidmask[i1] and self._datamask[i1]) or  (self._fidmask[i2] and self._datamask[i2])):
                        raise ValueError('masks incompatible with bands at positions %i and %i'%(i1,i2))
                else:
                    i = lbns.index(b)
                    assert (isinstance(i,int) and i >=0 and i < self._nb),(i,b)
                    if not self._datamask[i]:
                        raise ValueError('data mask incompatible with band at position %i'%i)
                    if not self._fidmask[i]:
                        raise ValueError('fiducial mask incompatible with band at position %i'%i)
            
            
            #now do the calculations
            fids,dats = [],[]
            
            for b in self._offbands:
                if isinstance(b,tuple):
                    fids.append(self._fdatadict[b[0]]-self._fdatadict[b[1]])
                    dats.append(self._data[lbns.index(b[0])]-self._data[lbns.index(b[1])])
                else:
                    fids.append(self._fdatadict[b]+self._dmod)
                    #fids.append(self.getFiducialData(b))
                    dats.append(self._data[lbns.index(b)])
            fids,data = np.array(fids,copy=False).T,np.array(dats,copy=False).T
            
        #TODO:interpolate along fiducials to a reasonable density
        #dims here: fids = nfXnb and data = ndXnb
        datat = np.tile(data,(fids.shape[0],1,1)).transpose((1,0,2))
        diff = (fids - datat)
        
        if self._offws is not None:
            ws = self._offws
            m = ws<0
            if np.any(m):
                ws = ws.copy().astype(float)
                ws[m] = -ws[m]/(diff[:,:,m].max(axis=1).max(axis=0)-diff[:,:,m].min(axis=1).min(axis=0))
            diff *= ws
        
        sepsq = np.sum(diff*diff,axis=2)
        sepsq = np.min(sepsq,axis=1) #CMD offset
        if self._locw and self._locs is not None:
            locsep = self.locs.T-self.center[:self.locs.shape[0]]
            sepsq = sepsq + self._locw*np.sum(locsep*locsep,axis=1)
        self._offsets = sepsq**0.5
        
    def getBand(self,i):
        if isinstance(i,int):
            return self._banddict[self._bandnames[i]]
        elif isinstance(i,basestring):
            return self._banddict[i]
        raise TypeError('Band indecies must be band names or integers')
    
    def getFiducialData(self,i,fidname=None):
        """
        Retreives fiducial data values in a particular band.  
        
        i is either an index into the bands, or a band name.
    
        fidname specifies which fiducial sequence to use (or if None,
        the full array for the requested band will be returned)
        """
        if fidname is None:
            mask = np.arange(self._nfid)
        else:
            mask = self._fidnamedict[fidname]
        if isinstance(i,int):
            return self._fdatadict[self._bandnames[i]][mask]+self._dmod
        elif isinstance(i,basestring):
            if '-' in i:
                b1,b2 = i.split('-')
                b1,b2 = b1.strip(),b2.strip()
                if (b1 not in self._fmaskd and b2 not in self._fmaskd) or (not self._fmaskd[b1] and not self._fmaskd[b2]):
                    raise ValueError('%s and %s are not valid fiducial bands'%(b1,b2))
                
                c = self._fdatadict[b1][mask] - self._fdatadict[b2][mask]
                return c
            else:
                b = i.strip()
                if b not in self._fmaskd or not self._fmaskd[b]:
                    raise ValueError('%s is not a valid fiducial band'%b)
                return self._fdatadict[b][mask]+self._dmod
        else:
            return self._fdatadict[i][mask]+self._dmod
        raise TypeError('Band indecies must be band names or integers')
    
    @property
    def fidnames(self):
        return self._fidnamedict.keys()
    
    @property
    def bandnames(self):
        return self._bandnames
    
    @property
    def validbandnames(self):
        return [b for i,b in enumerate(self._bandnames) if self._datamask[i] and self._fidmask[i]]
    
    @property
    def validbandinds(self):
        return [i for i,b in enumerate(self._bandnames) if self._datamask[i] and self._fidmask[i]]
    
    def getData(self,i=None):
        if i is None:
            return self._data
        elif isinstance(i,int):
            return self._data[i]
        elif isinstance(i,basestring):
            lbn=list(self._bandnames)
            if '-' in i:
                b1,b2 = i.split('-')
                b1,b2=b1.strip(),b2.strip()
                if b1 not in self.validbandnames and b2 not in self.validbandnames:
                    raise ValueError('%s and %s are not valid bands'%(b1,b2))
                #if b2 not in self.validbandnames:
                #    raise ValueError('%s is not a valid band'%b2)
                return self._data[lbn.index(b1)] - self._data[lbn.index(b2)]
            else:
                if i not in self.validbandnames:
                    raise ValueError('%s is not a valid band'%i)
                return self._data[lbn.index(i)]
        else:
            raise ValueError('unrecognized data type')
    def setData(self,val):
        """
        This loads the data to be compared to the fiducial sequence - the input
        can be an array (B X N) or an 
        """
        from operator import isMappingType
        if isMappingType(val):
            bandd,datad,maskd = self._dataAndBandsToDicts(val.keys(),np.array(val.values(),ndmin=2))
            
            arrinddict = {}
            currnd = None
            for b,d in datad.iteritems():
                if b in self._bandnames:
                    bi = list(self._bandnames).index(b)
                    if maskd[b]:
                        if currnd and len(d) != currnd:
                            raise ValueError('Band %s does not match previous band sizes'%b)
                        currnd = len(d)
                        arrinddict[bi] = d
                else:
                    from warnings import warn
                    warn('input data includes band %s not present in fiducial'%b)
                        
            mask,dat = [],[]
            for delem in [arrinddict.pop(i,None) for i in range(self._nb)]:
                if delem is None:
                    mask.append(False)
                    dat.append(np.NaN*np.empty(currnd))
                else:
                    mask.append(True)
                    dat.append(delem)
                    
            self._data = np.array(dat)
            self._datamask = np.array(mask,dtype=bool)
            self._nd = currnd
            
        else:
            if len(val) != self._nb:
                raise ValueError('input sequence does not match number of bands')
            self._data = arr = np.array(val,ndmin=2)
            self._nd = arr.shape[1]
            self._datamask = np.ones(arr.shape[0],dtype=bool)
            
        #check to make sure the offset bands are ok for the new data
        try:
            self.offsetbands = self.offsetbands
        except ValueError:
            print 'new data incompatible with current offset bands - setting to use all'
            self._offbands = None
        
        #need to recalculate offsets with new data
        self._offsets = None 
        #similarly for locs and names - assume new data means new objects
        self._locs = None
        self._dnames = None
    data = property(getData,setData,doc='\nsee ``getData`` and ``setData``')
        
        
    def _getLocs(self):
        return self._locs
    def _setLocs(self,val):
        
        if val is None:
            self._locs = None
        else:
            arr = np.array(val,ndmin=2)
            if arr.shape[0] < 2 or arr.shape[0] > 3:
                raise ValueError('locations must have either two componets or three')
            if arr.shape[1] != self._nd:
                raise ValueError('location data must match photometric data')
            self._locs = arr
    locs = property(_getLocs,_setLocs,doc="""
    The physical location of the data for prioritizing based on distance
    """)
    def _getDName(self):
        if self._dnames is None:
            return [str(i+1) for i in range(self._nd) ]
        else:
            return list(self._dnames)
    def _setDName(self,val):
        if len(val) != self._nd:
            raise ValueError('number of names do not match number of data points')
        self._dnames = tuple(val)
    datanames = property(_getDName,_setDName,doc="""
    Names for the objects - default is in numerical order starting at 1
    """)
    
    def _getDGrp(self):
        return self._dgrp
    def _setDGrp(self,val):
        if val is None:
            self._dgrp = None
            return
        if len(val) != self._nd:
            raise ValueError('number of group numbers do not match number of data points')
        self._dgrp = np.array(val,dtype=int)
    datagroups = property(_getDGrp,_setDGrp,doc="""
    Group numbers for the objects - will be cast to ints.
    Typical meanings:
    0 -> don't use in spectrum
    -1 -> alignment or guide star
    """)
    
    def _getCen(self):
        return self._cen
    def _setCen(self,val):
        if val is None:
            self._cen = (0,0,0)
        else:
            if len(val) == 2:
                self._cen = (val[0],val[1],0)
            elif len(val) == 3:
                self._cen = tuple(val)
            else:
                raise ValueError('locations must have either two componets or three')
    center = property(_getCen,_setCen,doc="""
    The center of the data for prioritizing based on distance
    """)
            
    def _getOffBands(self):
        if self._offbands is None:
            return None
        return ['%s-%s'%e if isinstance(e,tuple) else e for e in self._offbands]
    
    def _setOffBands(self,val):
        from operator import isSequenceType
        
        if val is None:
            self._offbands = None
        else:
            op = []
            lbn = list(self._bandnames)
            for i,v in enumerate(val):
                if isinstance(v,basestring):
                    try:
                        vs = v.split('-')
                        if len(vs) == 2:
                            b1,b2 = vs
                            op.append((lbn.index(b1.strip()),lbn.index(b2.strip())))
                        else:
                            op.append(lbn.index(v.strip()))
                    except:
                        raise ValueError('uninterpretable band type in position %i'%i)
                elif isinstance(v,int):
                    op.append(v)
                else:
                    try:
                        b1,b2 = v
                        if isinstance(b1,basestring):
                            op.append(lbn.index(b1))
                        else:
                            op.append(b1)
                        if isinstance(b2,basestring):
                            op.append(lbn.index(b2))
                        else:
                            op.append(b2)
                    except:
                        raise ValueError('uninterpretable band type in position %i'%i)
                    
            bandkeys = []
            for v in op:
                if isinstance(v,tuple):
                    bandkeys.append((self._bandnames[v[0]],self._bandnames[v[1]]))
                else:
                    bandkeys.append(self._bandnames[v])
    
            self._offbands = bandkeys
        #need to recalculate offsets with new data
        self._offsets = None 
        #offset weights are now invalid
        self._offws = None
        
        
    offsetbands = property(_getOffBands,_setOffBands,doc="""
    this selects the bands or colors to use in determining offsets
    as a list of band names/indecies (or for colors, 'b1-b2' or a 2-tuple)
    if None, all the bands will be used directly.
    """)
        
    def getOffsets(self):
        """
        computes and returns the CMD offsets
        """
        if self._data == None:
            raise ValueError("data not set - can't compute offsets")
        if self._offsets is None:
            self._calculateOffsets()
        return self._offsets
    def clearOffsets(self):
        """
        clears the offsets for setting of properties without recalculation
        """
        self._offsets = None
    
    def _getDist(self):
        return distance_from_modulus(self.distmod)
    def _setDist(self,val):
        self.distmod = distance_modulus(val)
    distance = property(_getDist,_setDist,doc="""
    The distance to use in calculating the offset between the fiducial values 
    and the data
    """)
    def _getDMod(self):
        return self._dmod
    def _setDMod(self,val):
#        if val == 0:
#            self._fdatadict = self._dm0fdd
#        else:
#            newd = {}
#            for b,arr in self._dm0fdd.iteritems():
#                newd[b] =  arr+val
#            self._fdatadict = newd
        self._dmod = val
        #need to recalculate offsets with new data
        self._offsets = None
    distmod = property(_getDMod,_setDMod,doc="""
    The distance modulusto use in calculating the offset 
    between the fiducial values and the data 
    """)
    def _getOffWs(self):
        if self._offws is None:
            return np.ones(len(self._offbands))
        else:
            return self._offws
    def _setOffWs(self,val):
        if val is None:
            self._offws = None
            return
        if (self._offbands is None and len(val) != self._nb) or (len(val)!=len(self._offbands)):
            raise ValueError('offset length does not match weights')
        self._offws = np.array(val)
        self._offsets = None
    offsetweights = property(_getOffWs,_setOffWs,doc="""
    Weights to apply to bands while calculating the offset. 
    Note that this will be invalidated whenever the offsetbands
    change
    
    if None, no weighting will be done.  otherwise, offsetweights
    must be a sequence of length equal to offsetbands.  positibe
    values will be interpreted as raw scalings (i.e. offset=w*(dat-fid) )
    and negative values will be compared to the total range (i.e. 
    offset=-w*(dat-fid)/(max-min) )
    """)
    def _getLocW(self):
        if self._locw is None:
            return 0
        return self._locw
    
    def _setLocW(self,val):
        self._locw = val
        self._offsets = None
    locweight = property(_getLocW,_setLocW,doc="""
    Weights to apply to the location while calculating the offset. 
    """)
    
    
    def plot(self,bx,by,clf=True,skwargs={},lkwargs={}):
        """
        this plots the two bands specified by bx and by
        
        clf determines if the figure should be cleared before plotting
        
        skwargs are passed into pylab.scatter and lkwargs are passed
        into pylab.plot
        """
        from matplotlib import pyplot as plt
            
        fx,fy = self.getFiducialData(bx),self.getFiducialData(by)
        dx,dy = self.getData(bx),self.getData(by)
        
        if clf:
            plt.clf()
        
        if 'c' not in skwargs:
            skwargs['c'] = 'b'
        if 'c' not in lkwargs:
            skwargs['c'] = 'r'
        plt.scatter(dx,dy,**skwargs)
        plt.plot(fx,fy,**lkwargs)
        
        if isinstance(bx,int):
            bx = self._bandnames[bx]
        if isinstance(by,int):
            by = self._bandnames[by]
        
        plt.xlabel(bx)
        plt.ylabel(by)
        
        if '-' not in bx:
            plt.xlim(*reversed(plt.xlim()))
        if '-' not in by:
            plt.ylim(*reversed(plt.ylim()))
        
    
def _CMDtest(nf=100,nd=100,xA=0.3,yA=0.2,plot=False):
    from matplotlib import pyplot as plt
    x=np.linspace(0,1,nf)
    y=5-(x*2)**2+24
    x=x*0.6+.15
    
    fdi = (np.random.rand(nd)*nf).astype(int)
    dx = x[fdi]+xA*(2*np.random.rand(nd)-1)
    dy = y[fdi]+yA*np.random.randn(nd)
    if plot:
        plt.clf()
        plt.plot(x,y,c='r')
        plt.scatter(dx,dy)
        plt.ylim(*reversed(plt.ylim()))
        
    cmda = CMDAnalyzer((x,y),('g-r','r'),{'a':np.arange(50),'b':np.arange(50)+50})
    cmda.setData({'g-r':dx,'r':dy})
    cmda.center = (10,10)
    cmda.locs = np.random.randn(2,nd)/10+cmda.center[0]
    cmda.offsetbands=['g-r','r']
        
    return x,y,dx,dy,cmda


class PhotometryBase(object):
    """
    A base class for objects that perform photometric measurements from
    an image.  
    
    Subclasses should:
    
    * override `_compute` to take an image and perform the computation
      steps necessary to generate the parameters
    * pass any outputs through `_fluxToMag` to do flux>mag conversion
    
    all outputs will be in magnitudes with the zeropoint set by the
    `zeropoint` attribute if `usemags` is True, otherwise they will
    be 
    """
    __metaclass__ = ABCMeta
    
    image = None
    psf = None
    loc = None
    _band = None
    _zpt = None
    
    @abstractmethod
    def _compute(self,image,loc,psf):
        """
        This method should be overridden in the base class to perform
        the photometric measurement given an image, a psf (may be None) 
        and a location (may be None).  Should set:
        
        * `self._totalflux` - total flux
        
        whatever it returns will be returned from `computePhotometry`
        """
        raise NotImplementedError
    
    def computePhotometry(self,**kwargs):
        """
        Performs the primary computationally-intensive phase of the 
        photometric measurement.  
        
        kwargs will be set as attributes on this object
        """
        for k,v in kwargs.iteritems():
            setattr(self,k,v)
        
        if self.image is None:
            raise ValueError('no image provided for photometry')
        
        return self._compute(image=self.image,loc=self.loc,psf=self.psf)
    
    def _fluxToMag(self,flux):
        """
        internally used to convert fluxes into the appropriate output form
        """
        return _flux_to_mag(flux) + self.zeropoint
    
    def _magToFlux(self,mag):
        return _mag_to_flux(mag - self.zeropoint)
        
    
    def _getZeropoint(self):
        if self._zpt is None:
            if hasattr(self.image,'zeropoint'):
                return self.image.zeropoint
            else:
                return 0
        else:
            return self._zpt
    def _setZeropoint(self,val):
        self._zpt = val
    zeropoint = property(_getZeropoint,_setZeropoint,doc='If set to None, zeropoint comes from image')
    
    def _getBand(self):
        if self._band is None and hasattr(self.image,'band'):
            return self.image.band
        return self._band
    def _setBand(self,val):
        self._band = val
    band = property(_getBand,_setBand,doc=None)
    
    @property
    def totalflux(self):
        return self._totalflux
    
    @property
    def totalmag(self):
        return self._magToFlux(self._totalflux)
    
    
class PointSpreadFunction(object):
    """
    Represents a Point Spread Function (PSF) for an image. 
    
    Coordinate system is always in pixel coordinates 
    """
    __metaclass__ = ABCMeta
        
    @abstractmethod
    def convolve(self,arr2d,background=0):
        """
        Convolve this psf with the supplied 2D Array.
        
        background sets the value of the background to assume around the
        edges. 
        """
        raise NotImplementedError
        
    @abstractmethod
    def fit(self,arr2d):
        """
        set the PSF parameters from the provided 2D image array 
        """
        raise NotImplementedError
    
class KernelPointSpreadFunction(PointSpreadFunction):
    def __init__(self,kernelarr2d):
        self.kernel = kernelarr2d
        self.fftconvolve = False
        self.convmode = 'constant'
        
    def _getKernel(self):
        return self._kernel
    def _setKernel(self,val):
        kernel = np.array(val)
        if len(self.kernel.shape) != 2:
            raise ValueError('Supplied kernel is not 2D')
        self._kernel = kernel
        self._Kft = None
    kernel = property(_getKernel,_setKernel,doc=None)        
    
    def convolve(self,arr2d,background=0):
        if self.fftconvolve:
            fft = np.fft.fftn
            ifft = np.fft.ifftn
            
            if self._Kft is None:
                self._Kft = K = fft(self._kernel)
            else:
                K = self._Kft
            return np.fft.ifftshift(ifft(fft(arr2d)*K))
        else:
            from scipy.ndimage import convolve
            return convolve(arr2d,self._kernel,mode=self.convmode,cval=background)
    
    def fit(self,arr2d):
        self.kernel = arr2d
    
class ModelPointSpreadFunction(PointSpreadFunction):
    def __init__(self,model):
        self.model = model
        self.convsize = None
        self.convsampling = None
        self.convmode = 'constant'
        
    def _getModel(self):
        return self._mod
    def _setModel(self,val):
        from .model import get_model_instance,FunctionModel2DScalar
        self._mod = get_model_instance(val,FunctionModel2DScalar)
        self._mod.incoordsys = 'polar'
    model = property(_getModel,_setModel,doc=None)
    
    
    def convolve(self,arr2d,background=False):
        from scipy.ndimage import convolve
        if self.convsize is None:
            nx,ny = arr2d.shape
        else:
            if np.isscalar(self.convsize):
                nx = ny = self.convsize
            else:
                nx,ny = self.convsize
            
        k = self.model.pixelize(-nx/2,nx/2,-ny/2,ny/2,nx,ny,sampling=self.convsampling)
        return convolve(arr2d,k,mode=self.convmode,cval=background)
    
    def fit(self,arr2d,**kwargs):
        """
        kwargs are passed into the model `fitData` method
        """
        nx,ny = arr2d.shape
        x = np.arange(nx)-nx/2
        y = np.arange(nx)-ny/2
        
        self.model.fitData([x,y],arr2d,**kwargs)
    
class GaussianPointSpreadFunction(PointSpreadFunction):
    def __init__(self,sigma=1):
        self.sigma = sigma
        self.convmode = 'constant'
    
    __fwhmtosig = np.sqrt(8*np.log(2))    
    def _getFwhm(self):
        return self.__fwhmtosig*self.sigma
    def _setFwhm(self,val):
        self.sigma = val/self.__fwhmtosig
    fwhm = property(_getFwhm,_setFwhm,doc=None)
    
    
    def convolve(self,arr2d,background=0):
        from scipy.ndimage import gaussian_filter
        return gaussian_filter(arr2d,self.sigma,mode=self.convmode,cval=background)
    
    def fit(self,arr2d):
        """
        fits the image using a moment analysis - returns
        (cenx,ceny),[sigx,sigy]
        """
        from .utils import moments
        x0,y0 = moments(arr2d,1)
        sx,sy = moments(arr2d,2)
        self.sigma = (sx+sy)/2
        
        return (x0,y0),[sx,sy]
    
#class AperturePhotometry(object):
#    """
#    This class generates radial surface brightness profiles from
#    various sources and provides plotting and fitting of those 
#    profiles   
    
#    specifying a Band results in use of that band's zero point
    
#    apertures can be a positive number to specify a fixed number of
#    apertures, a negative number specifying the spacing between
#    apertures in the specified units, or a sequence for fixed 
#    aperture sizes in the specified units
#    #TODO:TEST!
    
#    type can be 'circular' or 'elliptical'
#    """
#    def __init__(self,band = None,units='arcsec',apertures=10,usemags=True,type='circular'):
#        self.type = type
#        self.usemags = usemags #do this right
#        self.apertures = apertures
        
#        self.band = band
#        self.units = units
#        raise NotImplementedError('not usable yet')
        
        
#    def _getUnit(self):
#        return self._sunit
#    def _setUnit(self,unit):
#        if unit == 'arcsec':
#            self._sunit = 'arcsec'
#            self._unitconv = 1
#        elif unit == 'radians' or unit == 'sterradians' or unit == 'sr':
#            self._sunit = 'radians'
#            self._unitconv = 60*60*360/2/pi
#        elif 'deg' in unit:
#            self._sunit = 'degrees'
#            self._unitconv = 60*60
#        else:
#            raise ValueError('unrecognized unit')
#    units = property(_getUnit,_setUnit)
    
#    def _getBand(self):
#        return self._band
#    def _setBand(self,band):
#        if isinstance(band,basestring):
#            global bands
#            self._band = bands[band]
#        elif isinstance(band,Band):
#            self._band = band
#        elif band is None:
#            self._band = None
#        else:
#            raise ValueError("input band was not a Band or registered band name")
#    band = property(_getBand,_setBand)
        
#    def _magorfluxToFlux(self,magorflux):
#        if self.usemags:
#            zptf = 0 if (self._band is None) else self.band.zptflux
#            flux = zptf*_mag_to_flux(magorflux)
#        else:
#            flux = magorflux
#        return flux 
    
#    def _outFluxConv(self,flux):
#        """
#        convert flux back into mags if necessary and set zero-points appropriately
#        """
#        if self.usemags:
#            zptf = 0 if (self._band is None) else self.band.zptflux
#            out = _flux_to_mag(flux/zptf)
#        else:
#            out = flux
#        return out
    
#    def _fitPhot(self,x,y,flux,cen,aps=''):
#        """
#        input x,y,cen should be in un-converted units
#        """        
#        if self.type == 'circular':
                
#            r = np.sum((np.array((x,y),copy=False).T-cen)**2,axis=1)**0.5
            
#            if np.isscalar(self.apertures):
#                if self.apertures<0:
#                    rap = np.arange(0,np.max(r),self.apertures)+self.apertures
#                else:
#                    rap = np.linspace(0,np.max(r),self.apertures+1)[1:]
#            else:    
#                rap = np.array(self.apertures)
                
#            r,rap = self._unitconv*r,self._unitconv*rap
            
#            apmask = r.reshape((1,r.size)) <= rap.reshape((rap.size,1)) #new shape (rap.size,r.size)
#            fluxap = np.sum(apmask*flux,axis=1)
#            fluxann = fluxap-np.roll(fluxap,1)
#            fluxann[0] = fluxap[0] 
#            A = pi*rap*rap
#            A = A - np.roll(A,1)
#            A[0] = pi*rap[0]**2
#            sbap = fluxann/A
#            #TODO:better/correct SB?
            
#    #        valid = sbap > 0
#    #        self.enclosed = self._outFluxConv(fluxap[valid])
#    #        self.mu = self._outFluxConv(sbap[valid])
#    #        self.rap = rap   [valid]

#            highest = np.min(np.where(np.concatenate((sbap,(0,)))<=0))
#            self.enclosed = self._outFluxConv(fluxap[:highest])
#            self.mu = self._outFluxConv(sbap[:highest])
#            self.rap = rap[:highest]/self._unitconv
#        else:
#            raise NotImplementedError('elliptical fitting not supported yet')
        
#    def pointSourcePhotometry(self,x,y,magorflux,cen='weighted'):
#        """
#        x,y,and magorflux are arrays of matching shape 
        
#        cen is 'centroid' to use the raw centroid, 'weighted' to weight the
#        centroid by the flux, or a 2-sequence to manually assign a center
#        """
#        x = np.array(x,copy=False)
#        y = np.array(x,copy=False)
#        flux = self._magorfluxToFlux(np.array(magorflux,copy=False))
        
#        if x.shape != y.shape != flux.shape:
#            raise ValueError("shapes don't match!")
        
#        if isinstance(cen,basestring):
#            if cen == 'centroid':
#                cen = np.array((np.average(x),np.average(y)))
#            elif cen == 'weighted':
#                cen = np.array((np.average(x,weights=flux),np.average(y,weights=flux)))
#            else:
#                raise ValueError('unrecognized cen string')
#        else:
#            cen = np.array(cen,copy=False)
#            if cen.shape != (2,):
#                raise ValueError('cen not a 2-sequence')
        
#        return self._fitPhot(x,y,flux,cen)
    
#    def imagePhotometry(self,input,platescale,cen='weighted'):
#        """
#        input must be a 2D array or a CCD
        
#        platescale is linear units as specified by the units property per pixel
#        """
#        from .ccd import CCDImage
        
#        if isinstance(input,CCDImage):
#            raise NotImplementedError
        
#        magorflux = np.array(input,copy=True)
#        if len(magorflux.shape) != 2:
#            raise ValueError('input not 2D')
#        flux = self._magorfluxToFlux(magorflux)
        
#        x,y = np.meshgrid(np.arange(flux.shape[0]),np.arange(flux.shape[1]))
        
#        return self._fitPhot(x.T+1,y.T+1,flux,cen)
    
#    def plot(self,fmt='o-',logx=False,**kwargs):
#        """
#        plot the surface brightness profile (must be generated)
#        using matplotlib
        
#        kwargs go into matplotlib.plot
#        """
        
#        from matplotlib import pyplot as plt
#        if logx:
#            plt.semilogx()
#        if 'fmt':

#            plt.plot(self.rap,self.mu,fmt,**kwargs)
        
#        plt.ylim(*reversed(plt.ylim()))
        
#        plt.xlabel('$r \\, [{\\rm arcsec}]$')
#        if self.band is None:
#            plt.ylabel('$\\mu$')
#        else:
#            plt.ylabel('$\\mu_{%s}$'%self.band.name)
            
class IsophotalEllipse(object):
    """
    Generates best fit ellipses for an image
    """
    
    def __init__(self,imdata,isolevel=None):
        self.imdata = imdata
        
        if isolevel is None:
            self.isolevel = np.mean(self._imdata)
            
        self._fitted = False
        
        self._a = 1 #major axis
        self._b = 1 #minor axis
        self._phi = 0 #rotation
        self._fitpoints = 100
        self._x0 = self._imdata.shape[0]/2.0
        self._y0 = self._imdata.shape[1]/2.0
        self._fixcen = False
        self._keeppix = .01 #this could be auto-tuned but generally seems to give a good answer
            
    def _getImdata(self):
        return self._imdata
    def _setImdata(self,val):
        val = np.array(val)
        if len(val.shape) != 2:
                raise ValueError('2D image data not provided')
        self._imdata = val
        self._fitted = False
    imdata = property(_getImdata,_setImdata)
    
    def _getFitpoints(self):
        return self._fitpoints
    def _setFitpoints(self):
        self._fitted = False
    nfitpoints = property(_getFitpoints,_setFitpoints)
    
    def _getLevel(self):
        return self._level
    def _setLevel(self,val):
        self._level = val
        self._fitted = False
    isolevel = property(_getLevel,_setLevel)
    
    def _getFixcen(self):
        return self._fixcen
    def _setFixcen(self,val):
        self._fixcen = val
        self._fitted = False
    fixcenter = property(_getLevel,_setLevel)
    
    def _getKeeppix(self):
        return self._keeppix
    def _setKeeppix(self,val):
        self._keeppix = val
        self._fitted = False
    keeppix=property(_getKeeppix,_setKeeppix,doc="""
    The pixels to keep for fitting, either as a fraction of the total
    (if less than 1) or as a fixed number of pixels (if greater than 1)
    """)
        
    def polar(self,npoints):
        """
        returns the full ellipse as an (rho,theta) tuple
        """
        th = np.linspace(0,2*pi,npoints)
        return self._th2r(th),th
    
    def _th2r(self,th):
        a,b = self._a,self._b
        aterm = a*np.sin(th-self._phi)
        bterm = b*np.cos(th-self._phi)
        return a*b*(aterm*aterm+bterm*bterm)**-0.5
        
    def cartesian(self,npoints):
        """
        returns the full ellipse as an (x,y) tuple
        """
        r,th = self.polar(npoints)
        x = r*np.cos(th)
        y = r*np.sin(th)
        return x+self._x0,y+self._y0
        
    def _fitEllipse(self):
        from scipy.optimize import leastsq
        from scipy.ndimage import map_coordinates
        
        diff = self._imdata - self._level
        maxdiff = max(np.max(diff)**2,np.min(diff)**2)
        maxsz = np.min(diff.shape)/2
        
        xi,yi = np.meshgrid(np.arange(diff.shape[0]),np.arange(diff.shape[1]))
        xi = xi.T.ravel()
        yi = yi.T.ravel()
        sorti = np.argsort(diff.ravel()**2)
        nkeep = self._keeppix*xi.size if self._keeppix <= 1 else self._keeppix
        xi = xi[sorti][:nkeep]
        yi = yi[sorti][:nkeep]
        keepdiffs = diff.ravel()[sorti][:nkeep]
        
        def f(vals,self):
            self._a = np.abs(vals[0])
            self._b = np.abs(vals[1])
            self._phi = vals[2]
            if not self._fixcen:
                self._x0 = vals[3]
                self._y0 = vals[4]
            
            xo,yo = xi-self._x0,yi-self._y0
            
            
            sep = (xo*xo+yo*yo)**0.5 - self._th2r(np.arctan2(yo,xo))
            
            return sep*keepdiffs
                
        v0 = [self._a,self._b,self._phi] if self._fixcen else [self._a,self._b,self._phi,self._x0,self._y0]
        diag = [1/diff.shape[0],1/diff.shape[1],1] if self._fixcen else [1/diff.shape[0],1/diff.shape[1],1,1,1]
        soln,cov,infodict,mesg,ier = leastsq(f,v0,args=(self,),full_output=1,diag = diag)
        
        if not ier ==1:
            print 'Possible fit problem:',mesg
        if self._fixcen:
            self._a,self._b,self._phi = soln
        else:
            self._a,self._b,self._phi,self._x0,self._y0 = soln
        self.lastier = ier
        self.lastmesg = mesg
        
        if self._a < self._b:
            self._a,self._b = self._b,self._a
            self._phi += pi/2
            
        if self._phi >= 2*pi:
            self.phi -= 2*pi*np.floor(self.phi/2/pi)
        elif self._phi < 0:
            self._phi += 2*pi*np.floor(-self.phi/2/pi)
        
        self._fitted = True
        
    def getDiff(self,fractional=False):
        """
        returns the difference between the fitted ellipse and the flux level
        """
        from scipy.ndimage import map_coordinates
        if not self._fitted:
            self._fitEllipse()
        
        if fractional:
            res = map_coordinates(self._imdata-self._level,self.cartesian(self._fitpoints),mode='nearest')
        else:
            res = map_coordinates(1-self._level/self._imdata,self.cartesian(self._fitpoints),mode='nearest')
        return res
    
    @property    
    def e(self):
        """
        returns eccentricity
        """
        if not self._fitted:
            self._fitEllipse()
        ratio = self._a/self._b
        return (1-ratio*ratio)**0.5
    
    @property
    def axisratio(self):
        """
        returns major/minor axis ratio
        """
        if not self._fitted:
            self._fitEllipse()
        return self._a/self._b
    
    @property    
    def major(self):
        """
        returns major axis length in image coordinates
        """
        if not self._fitted:
            self._fitEllipse()
        return self._a
    a = major
    
    @property    
    def minor(self):
        """
        returns major axis length in image coordinates
        """
        if not self._fitted:
            self._fitEllipse()
        return self._b
    b = minor
    
    @property
    def phi(self):
        """
        returns rotation angle clockwise from x-axis
        """
        if not self._fitted:
            self._fitEllipse()
        return self._phi
    
    def plot(self,clf=True):
        from matplotlib import pyplot as plt
        
        if not self._fitted:
            self._fitEllipse()
        
        if clf:
            plt.clf()
            
        plt.imshow(self.imdata.T)
        plt.plot(*self.cartesian(self._fitpoints),**{'c':'k'})
        
        
        
class ModelPhotometry(object):
    """
    Represents a astropysics.models 2D model of the photometric flux
    from an image.
    
    The normparam attribute defines which model parameter sets
    the overall scaling.  If None, the scaling methods will not function.
    """
    def __init__(self,model,magzpt=0):
        from . import models
        self.model = model
        self.magzpt = magzpt #TODO: get from a band/image?
        
    normparam = ('A','Ae')
        
    def _getModel(self):
        return self._model
    def _setModel(self,val):
        from .models import get_model_instance,FunctionModel2DScalar
        self._model = get_model_instance(val,FunctionModel2DScalar)
    model = property(_getModel,_setModel)
    
    def fitPhot(self,im):
        raise NotImplementedError
    
    def _getTotalflux(self):
        return self._model.integrateCircular(np.inf)
    def _setTotalflux(self,val):
        if self.normparam is None:
            raise NotImplementedError('No auto-mult normalization built in yet')
        elif isinstance(self.normparam,basestring):
            normparam = (self.normparam,)
        else:
            normparam = self.normparam
        normparam = [p in self._model.params for p in self.normparam]
        if np.sum(normparam)!=1:
            raise ValueError('could not find unique normalization parameter for this model')
        normparam = self.normparam[normparam.index(True)]
        
        currflux = self._getTotalflux()
        newnorm = getattr(self._model,normparam)*val/currflux
        setattr(self._model,normparam,newnorm)
    totalflux = property(_getTotalflux,_setTotalflux,doc=None)
    
    #TODO:fix for flux being per unit arcsec
    def _getTotalmag(self):
        return self.magzpt-2.5*np.log10(self.totalflux)
    def _setTotalmag(self,val):
        self.totalflux = 10**((val-self.magzpt)/-2.5)
    totalmag = property(_getTotalmag,_setTotalmag,doc=None)
    
    def getHalfLightRadius(self,**kwargs):
        """
        computes the radius in which half the flux is encircled.  kwargs
        are passed into `FunctionModel2DScalar.getFluxRadius`
        """
        return self.model.getFluxRadius(0.5,True,**kwargs)
    
    def simulate(self,pixels,scale=1,background=0,noise=None,psf=None,sampling=None):
        """
        Simulate how this model would appear on an image. 
        
        `pixels` is a 2D array that described the number of pixels to use, a
        scalar describing the number of pixels in both directions, or a sequence
        (npixx,npixy)
        
        `scale` is either a scalar giving the pixel scale or a 2-sequence
        (xscale,yscale) as units/pixel
        
        `background` is a scalar background level to assume or a 2D array of
        matching shape as the simulation
        
        `noise` can be:
        
        * None 
            no noise is applied
        * 'poisson'
            assumes the flux units are counts and returns the poisson noise
            model
        * a scalar
            uniform gaussian noise with the scalar as sigma
        * a 2d array
            gaussian noise with the value at each array point as sigma - array
            dimension must match simulation
        * a callable 
            will be called as noise(imagearr2d) and should return the value of
            the model+noise
          
        `psf` is the point-spread function to be applied before the noise,
        either a scalar/2-tuple for a gaussian PSF (specifies the FWHM in actual
        units, not pixels), or a 2D measured PSF to convolve with the model. If
        None, no psf will be included.
        
        `sampling` is passed into the model's `pixelize` method (see
        :meth:`FunctionModel2DScalar.pixelize`)
        """
        pixels = np.array(pixels,copy=False)
        if len(pixels.shape)==0:
            nx = ny = int(pixels)
        elif len(pixels.shape)==1:
            nx,ny = pixels.astype(int)
        elif len(pixels.shape)==2:
            nx,ny = pixels.shape
        else:
            raise ValueError('Invalid input for pixels')
        
        if np.isscalar(scale):
            xscale = yscale = scale
        else:
            xscale,yscale = scale
            
        xsize,ysize = nx*xscale,ny*yscale
        
        m = self._model
        modim = m.pixelize(-xsize/2,xsize/2,-ysize/2,ysize/2,nx,ny,sampling)
        
        if psf: #TODO:generalize PSF in superclass
            psf = np.array(psf,copy=False)
            if len(psf.shape)==0:
                from scipy.ndimage import gaussian_filter
                sigtofwhm = 2*np.sqrt(2*np.log(2))
                sigmas = (psf/xscale/sigtofwhm,psf/yscale/sigtofwhm)
                modim = gaussian_filter(modim,sigmas,mode='nearest')
                
            elif len(psf.shape)==1:
                from scipy.ndimage import gaussian_filter
                sigtofwhm = 2*np.sqrt(2*np.log(2))
                xpsf,ypsf = psf
                sigmas = (xpsf/xscale/sigtofwhm,ypsf/yscale/sigtofwhm)
                modim = gaussian_filter(modim,sigmas,mode='nearest')
                
            elif len(psf.shape)==2:
                from scipy.ndimage import convolve
                normpdf = psf/np.sum(psf)
                modim = convolve(modim,normpsf,mode='nearest')
                
            else:
                raise ValueError('input psf not valid')
            
        if background:
            modim += background
            
        if noise:
            if noise == 'poisson':
                modim = np.random.poisson(modim)
            elif callable(noise):
                modim = noise(modim)
            else:
                modim = np.random.normal(modim,noise)
                
        return modim
    


class SExtractor(object):
    """
    This class is an adaptor to the Sextractor program 
    (Bertin & Arnouts 1996, http://astromatic.iap.fr/software/sextractor/).
    Sextractor must be installed and on the path for this class to function.
    
    options are set by changing values in the options dictionary
    
    output parameters are chosen by setting True/False values in the 
    params dictionary

    """
    
    @staticmethod
    def _getSexDefaults():
        from subprocess import Popen,PIPE
        
        optinfo = {}
        opts = {}
        optorder = []
        parinfo = {}
        parorder = []
        
        try:
            pconf = Popen('sex -dd'.split(),executable='sex',stdout=PIPE,stderr=PIPE)
            pparm = Popen('sex -dp'.split(),executable='sex',stdout=PIPE,stderr=PIPE)
            pconf.wait()
            pparm.wait()
            confstr = pconf.communicate()[0]
            parmstr = pparm.communicate()[0]
        except OSError:
            raise OSError('Sextractor not found on system path')
        
        comm = ''
        newk = k = None
        for l in confstr.split('\n'):
            commenti = l.find('#')
            if commenti>0:
                ls = l[:commenti].split()
                if len(ls)>1:
                    newk = ls[0]
                    newoptval = ' '.join(ls[1:])
                elif len(ls)>0:
                    newk = ls[0]
                    newoptval = ''
                newcomm = l[commenti+1:].strip()
            else:
                newcomm = ''
                
            if newk:
                if k:
                    opts[k] = optval
                    optinfo[k] = comm
                    optorder.append(k)
                k = newk
                optval = newoptval
                
                newk = None
                comm = ''
            comm+=newcomm
              
        for l in parmstr.split('\n'):
            ls = l.split()
            if len(ls) > 1:
                k = ls[0].strip().replace('#','')
                unit = ls[-1].replace('[','').replace(']','') if '[' in ls[-1] else None
                info = ' '.join(ls[1:(-1 if unit else None)])
                parinfo[k] = (info,unit if unit else '')
                parorder.append(k)
        
        SExtractor._optinfo = optinfo
        SExtractor._defaultopts = opts
        SExtractor._optorder = optorder #TODO:OrderedDict for 2.7
        SExtractor._parinfo = parinfo
        SExtractor._parorder = parorder #TODO:OrderedDict for 2.7
    
    @staticmethod   
    def getOptInfo(aslist=False):
        """
        returns the  dictionary of input options and the associated information
        
        if aslist is True, returns an list of 
        """
        if aslist:
            return [(k,SExtractor._optinfo[k]) for k in SExtractor._optorder]
        else:
            return dict(SExtractor._optinfo)
    
    @staticmethod    
    def getParamInfo(aslist=False):
        """
        returns the dictionary of parameters and the associated information
        """
        if aslist:
            return [(k,SExtractor._parinfo[k]) for k in SExtractor._parorder]
        else:
            return dict(SExtractor._parinfo)
        
        
    def __init__(self,sexfile=None,parfile=None):
        
        opts = dict(SExtractor._defaultopts)
        pars = dict([(k,False) for k in  SExtractor._parinfo])
                    
        if sexfile:
            with open(sexfile) as f:
                for l in f:
                    commenti = l.find('#')
                    if l > -1:
                        l = l[:commenti]
                    ls = l.split()
                    if len(ls) > 1:
                        k = ls[0].strip()
                        if k not in opts:
                            raise ValueError('sexfile has invalid option %s'%k)
                        opts[k] = ls[1].strip()
            self.name = sexfile.replace('.sex','')
        else:
            self.name = 'astropysics_auto'
                        
                        
        pfile = opts['PARAMETERS_NAME'] if parfile is None else parfile
        if pfile != 'default.param':
            with open(pfile) as f:
                for l in f:
                    if not (l.strip().startswith('#') or l.strip()==''):
                        k = l.split()[0].strip()
                        if k not in pars:
                            raise ValueError('param file has invalid parameter %s'%k)
                        pars[k] = True
        if not np.any(pars.values()):
            #if no outputs provided, use defaults (from v 2.8.6)
            defs286 =['NUMBER','FLUX_ISO', 'FLUXERR_ISO', 'FLUX_AUTO', 'FLUXERR_AUTO', 'X_IMAGE', 'Y_IMAGE', 'FLAGS']
            for p in defs286:
                if p in pars:
                    pars[p] = True
        
        self.options = opts
        self.params = pars
        
        self.overwrite = False
        
            
    def getParamList(self):
        """
        returns a list of all selected parameters
        """    
        return [k for k in self._parorder if self.params[k]]
    
    def _saveFiles(self,fnbase):
        import os
        
        fnbase = fnbase.replace('.sex','')
        
        dir = os.path.split(fnbase)[0]
        dir = '' if dir=='' else dir+os.sep
        
        self.options['PARAMETERS_NAME'] = dir+fnbase+'.param'
        
        ostr = self._makeOptionStr()
        pstr = self._makeParamStr()
        
        with open(fnbase+'.sex','w') as f:
            f.write(ostr)
            
        with open(self.options['PARAMETERS_NAME'],'w') as f:
            f.write(pstr)
                
    
    def _makeOptionStr(self):
        ostr = [o+'\t'+str(self.options[o])+'\t # '+SExtractor._optinfo[o] for o in self._optorder]
        ostr = '\n'.join(ostr)
        return ostr
        
    def _makeParamStr(self):
        pstr = [p+'\t # '+str(SExtractor._parinfo[p]) for p in self._parorder if self.params[p]]
        pstr = '\n'.join(pstr)
        return pstr
    
    
    def getOptionList(self,incval=False):
        """
        returns a list of all options.  If incval is True, returns a list
        of tuples of the form (optionname,optionvalue)
        """
        if incval:
            return [(k,self.options[k]) for k in self._optorder]
        else:
            return [k for k in self._optorder]
        
    def sextractImage(self,detimfn,analysisimfn=None,mode='waiterror'):
        """
        writes configuration files and runs sextractor on the input image
        
        mode can be:
        
        * 'waiterror': waits for sextractor to finish, and raises an 
          SExtractorError if it does not complete sucessfully. stdout 
          and sterr are saved to self.lastout and self.lasterr (returns 0)
        * 'wait': waits for sextractor to finish and returns the return code
          stdout and sterr are saved to self.lastout and self.lasterr
        * 'proc': stars the processes but does not wait - returns the Popen 
          instance of the processes
        """
        from subprocess import Popen,PIPE
        from os.path import exists
        
        fnbase = self.name
        if not self.overwrite:
            fnbase = fnbase.replace('.sex','')
            if exists(fnbase+'.sex'):
                fns = fnbase.split('-')
                try:
                    i = int(fns[-1])
                    i+=1
                except ValueError:
                    i = 2
                if len(fns)<2:
                    fns.append(str(i))
                else:
                    fns[-1] = str(i)
                fnbase = '-'.join(fns)
            self.name = fnbase
                
        self._saveFiles(fnbase)
        if analysisimfn:
            clstr = 'sex {0} {1} -c {2}'.format(detimfn,analysisimfn,self.name+'.sex')
        else:
            clstr = 'sex {0} -c {1}'.format(detimfn,self.name+'.sex')
        proc = Popen(clstr.split(),executable='sex',stdout=PIPE,stderr=PIPE)
        
        if mode == 'waiterror' or mode =='wait':
            res = proc.wait()
            sout,serr = proc.communicate()
            
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                raise SExtractorError(serr,sout)
            return res
        elif mode == 'proc':
            return proc
        else:
            raise ValueError('unrecognized mode argument '+str(mode))
try:
    SExtractor._getSexDefaults()
except OSError:
    from warnings import warn
    warn('SExtractor not found on system - phot.SExtractor class will not function')
    
class SExtractorError(Exception):
    def __init__(self,*args):
        super(SExtractorError,self).__init__(*args)
            
        
#<------------------------conversion functions--------------------------------->
    
def UBVRI_to_ugriz(U,B,V,R,I,ugrizprimed=False):
    """
    transform UBVRcIc magnitudes to ugriz magnitudes as per Jester et al. (2005)
    """
    if not ugrizprimed:
        umg    =    1.28*(U-B)   + 1.13  
        #gmr    =    1.02*(B-V)   - 0.22  
        rmi    =    0.91*(R-I) - 0.20 
        rmz    =    1.72*(R-I) - 0.41
        g      =    V + 0.60*(B-V) - 0.12 
        r      =    V - 0.42*(B-V) + 0.11
        

    else:
        raise NotImplementedError
    
    return umg+g,g,r,r-rmi,r-rmz
    
def ugriz_to_UBVRI(u,g,r,i,z,ugrizprimed=False):
    """
    transform ugriz magnitudes to UBVRcIc magnitudes as per Jester et al. (2005)
    (note that z is unused)
    """
    if not ugrizprimed:
        UmB    =    0.78*(u-g) - 0.88 
        #BmV    =    0.98*(g-r) + 0.22 
        VmR    =    1.09*(r-i) + 0.22
        RmI  =    1.00*(r-i) + 0.21
        B      =    g + 0.39*(g-r) + 0.21
        V      =    g - 0.59*(g-r) - 0.01 

    else:
        raise NotImplementedError
    
    return UmB+B,B,V,V-VmR,V-VmR-RmI

def transform_dict_ugriz_UBVRI(d,ugrizprimed=False):
    ugriz = 'u' in d or 'g' in d or 'r' in d or 'i' in d or 'z' in d
    UBVRI = 'U' in d or 'B' in d or 'V' in d or 'R' in d or 'I' in d
    if ugriz and UBVRI:
        raise ValueError('both systems already present')
    if ugriz:
        u=d['u'] if 'u' in d else 0
        g=d['g'] if 'g' in d else 0
        r=d['r'] if 'r' in d else 0
        i=d['i'] if 'i' in d else 0
        z=d['z'] if 'z' in d else 0
        U,B,V,R,I=ugriz_to_UBVRI(u,g,r,i,z,ugrizprimed)
        if 'g' in d and 'r' in d:
            d['B']=B
            d['V']=V
            if 'u' in d:
                d['U']=U
            if 'i' in d and 'z' in d:
                d['I']=I
                d['R']=R
        else:
            raise ValueError('need g and r to do anything')
        
    if UBVRI:
        U=d['U'] if 'U' in d else 0
        B=d['B'] if 'B' in d else 0
        V=d['V'] if 'V' in d else 0
        R=d['R'] if 'R' in d else 0
        I=d['I'] if 'I' in d else 0
        u,g,r,i,z=UBVRI_to_ugriz(U,B,V,R,I,ugrizprimed)
        if 'B' in d and 'V' in d:
            d['g']=g
            d['r']=r
            if 'U' in d:
                d['u']=u
            if 'R' in d and 'I' in d:
                d['i']=i
                d['z']=z
        else:
            raise ValueError('need B and V to do anything')
        
        
def ML_ratio_from_color(c,color='B-V'):
    """
    uses Bell&DeJong 01 relations - note that this is intended for 
    normal spiral stellar pops
    
    
    color can either be a 'B-V','B-R','V-I','V-J','V-H',or 'V-K'
    
    returns tuple of mass-to-light ratios for each c as
    (mlrb,mlrv,mlrr,mlri,mlrJ,mlrH,mlrK)
    """
    bdj01table1 = np.array([[-.994,1.804,-.734,1.404,-.660,1.222,-.627,1.075,-.621,.794,-.663,.704,-.692,.652],
         [-1.224,1.251,-.916,.976,-.820,.851,-.768,.748,-.724,.552,-.754,.489,-.776,.452],
         [-1.919,2.214,-1.476,1.747,-1.314,1.528,-1.204,1.347,-1.040,0.987,-1.030,.870,-1.027,.800],
         [-1.903,1.138,-1.477,.905,-1.319,.794,-1.209,.700,-1.029,.505,-1.014,.442,-1.005,0.402],
         [-2.181,.978,-1.700,.779,-1.515,.684,-1.383,.603,-1.151,.434,-1.120,.379,-1.1,.345],
         [-2.156,.895,-1.683,.714,-1.501,.627,-1.370,.553,-1.139,.396,-1.108,.346,-1.087,.314]],
         [-0.942,1.737,-0.628,1.305,-0.520,1.094,-0.399,0.824,-0.261,0.433,-0.209,0.210,-0.206,0.135],
         [-0.976,1.111,-0.633,0.816,-0.523,0.683,-0.405,0.518,-0.289,0.297,-0.262,0.180,-0.264,0.138])
    colorrowmap={'B-V':0,'B-R':1,'V-I':2,'V-J':3,'V-H':4,'V-K':5}
    if color not in colorrowmap:
        raise ValueError('Unknown color')
    a = bdj01table1[colorrowmap[color],0::2].reshape((7,1))
    b = bdj01table1[colorrowmap[color],1::2].reshape((7,1))
    
    return tuple(10**(a+b*c))
    

def M_star_from_mags(B,V,R,I,color='B-V'):
    """
    uses Bell&DeJong 01 relations - note that this is intended for 
    normal spiral stellar pops
    
    
    color can either be a 'B-V','B-R','V-I', or 'mean'
    
    returns stellar mass as mean,(B-derived,V-derived,R-derived,I-derived)
    """    
    if color=='B-V':
        mlrb,mlrv,mlrr,mlri = ML_ratio_from_color(B-V,'B-V')[:4]
    elif color=='B-R':
        mlrb,mlrv,mlrr,mlri = ML_ratio_from_color(B-R,'B-R')[:4]
    elif color=='V-I':
        mlrb,mlrv,mlrr,mlri = ML_ratio_from_color(V-I,'V-I')[:4]
    elif color=='mean':
        bv=M_star_from_mags(B,V,R,I,'B-V')
        br=M_star_from_mags(B,V,R,I,'B-R')
        vi=M_star_from_mags(B,V,R,I,'V-I')
        return np.mean((bv[0],br[0],vi[0])),np.mean((bv[1],br[1],vi[1]),axis=0)
    else:
        raise ValueError('Unknown color')
    
    mstar=[]
    mstar.append(mag_to_lum(B,'B')*mlrb)
    mstar.append(mag_to_lum(V,'V')*mlrv)
    mstar.append(mag_to_lum(R,'R')*mlrr)
    mstar.append(mag_to_lum(I,'I')*mlri)
    
    return np.mean(mstar),tuple(mstar)
def ML_ratio_from_color_SDSS(c,color='g-r'):
    """
    uses Bell 03 relations derived from SDSS
    
    color can either be a 'u-g','u-r','u-i','u-z','g-r','g-i','g-z','r-i','r-z'
    
    returns tuple of mass-to-light ratios for each c as
    (mlrg,mlrr,mlri,mlrz,mlrJ,mlrH,mlrK)
    """
    b03table=np.array([[-0.221,  0.485, -0.099,  0.345, -0.053,  0.268, -0.105,  0.226,
        -0.128,  0.169, -0.209,  0.133, -0.26 ,  0.123],
       [-0.39 ,  0.417, -0.223,  0.299, -0.151,  0.233, -0.178,  0.192,
        -0.172,  0.138, -0.237,  0.104, -0.273,  0.091],
       [-0.375,  0.359, -0.212,  0.257, -0.144,  0.201, -0.171,  0.165,
        -0.169,  0.119, -0.233,  0.09 , -0.267,  0.077],
       [-0.4  ,  0.332, -0.232,  0.239, -0.161,  0.187, -0.179,  0.151,
        -0.163,  0.105, -0.205,  0.071, -0.232,  0.056],
       [-0.499,  1.519, -0.306,  1.097, -0.222,  0.864, -0.223,  0.689,
        -0.172,  0.444, -0.189,  0.266, -0.209,  0.197],
       [-0.379,  0.914, -0.22 ,  0.661, -0.152,  0.518, -0.175,  0.421,
        -0.153,  0.283, -0.186,  0.179, -0.211,  0.137],
       [-0.367,  0.698, -0.215,  0.508, -0.153,  0.402, -0.171,  0.322,
        -0.097,  0.175, -0.117,  0.083, -0.138,  0.047],
       [-0.106,  1.982, -0.022,  1.431,  0.006,  1.114, -0.052,  0.923,
        -0.079,  0.650, -0.148,  0.437, -0.186,  0.349],
       [-0.124,  1.067, -0.041,  0.78 , -0.018,  0.623, -0.041,  0.463,
        -0.011,  0.224, -0.059,  0.076, -0.092,  0.019]])
    colorrowmap=dict([(ci,i) for i,ci in enumerate(['u-g','u-r','u-i','u-z','g-r','g-i','g-z','r-i','r-z'])])
    
    a = b03table[colorrowmap[color],0::2].reshape((7,1))
    b = b03table[colorrowmap[color],1::2].reshape((7,1))
    
    return tuple(10**(a+b*c))
    
def M_star_from_mags_SDSS(u,g,r,i,z,J=None,H=None,K=None,color='mean'):
    mags = {'u':u,'g':g,'r':r,'i':i,'z':z,'J':J,'H':H,'K':K}
    
    if color=='mean':
        mstars = []
        for c in ['u-g','u-r','u-i','u-z','g-r','g-i','g-z','r-i','r-z']:
            mstars.append(M_star_from_mags_SDSS(u,g,r,i,z,J,H,K,c)[1])
        mstars = np.array(mstars)
        return mstars.mean(axis=0).mean(axis=0),mstars
    elif '-' in color:
        c1,c2 = color.split('-')
        mlrs = ML_ratio_from_color_SDSS(mags[c1]-mags[c2],color)
    else:
        raise ValueError('invalid Color')
    
    mstar=[mag_to_lum(mag,bstr)*mlr for mag,bstr,mlr in zip(mags.values(),mags.keys(),mlrs) if mag is not None]
    
    return np.mean(mstar),tuple(mstar)

def distance_modulus(x,intype='distance',dx=None,autocosmo=True):
    """
    compute the distance modulus given  a distance or redshift
    
    intype can be:
    * 'distance': will treat x as a distance in pc
    * 'redshift': x is treated as a redshift
     
    if autocosmo is true, the cosmological calculation will be 
    automatically  performed for z > 0.1.  If False, a 
    non-cosmological calculation will be done.  If it is 'warn', 
    a warning will be issued for any inputs where z>0.1
    
    """
    from .coords import cosmo_z_to_dist
    from .constants import H0,c
    
    c=c/1e5 #km/s
    cosmo=False
    
    if intype == 'distance':
        z=x/1e6*H0/c
        if dx is not None:
            dz = dx/1e6*H0/c
        else:
            dz = None
    elif intype == 'redshift':
        z=x
        x=z*c/H0
        if dx is not None:
            dz = dx
            dx = dz*c/H0
        else:
            dz = None
    else:
        raise ValueError('unrecognized intype')
    
    if autocosmo and np.any(z) > 0.1:
        if autocosmo == 'warn':
            from warnings import warn
            warn('redshift < 0.1 - cosmological calculation should be used')
        else:
            cosmo=True
    if cosmo:
        return cosmo_z_to_dist(z,dz,4)
    elif dx is None:
        return 5*np.log10(x)-5
    else:
        dm = 5*np.log10(x)-5
        ddm = 5*dx/x/np.log(10)
        return dm,ddm,ddm
    
def distance_from_modulus(dm):
    """
    compute the distance given the specified distance modulus.  Currently
    non-cosmological
    """
    return 10**(1+dm/5.0)
    

def abs_mag(appmag,x,intype='distance',autocosmo=True):
    """
    computes absolute magnitude from apparent magnitude and distance.
    See astro.phot.distance_modulus for details on arguments
    """
    from operator import isSequenceType
    if isSequenceType(appmag):
        appmag=np.array(appmag)
    
    distmod = distance_modulus(x,intype,None,autocosmo)
    return appmag-distmod

def app_mag(absmag,x,intype='distance',autocosmo=True):
    """
    computes apparent magnitude from absolute magnitude and distance.
    See astro.phot.distance_modulus for details on arguments
    """
    from operator import isSequenceType
    if isSequenceType(absmag):
        absmag=np.array(absmag)
        
    distmod = distance_modulus(x,intype,None,autocosmo)
    return absmag+distmod

def rh_to_surface_brightness(totalm,rh):
    """
    Compute the surface brightness given a half-light radius in arcsec.  Note
    that the totalm is the integrated magnitude, not just the half-light
    """
    return area_to_surface_brightness(totalm+2.5*np.log10(2),pi*rh*rh)

def area_to_surface_brightness(m,area):
    """
    Compute the surface brightness given a particular magnitude and area in 
    mag/sq arc seconds
    """
    return m+2.5*np.log10(area)

def intensity_to_flux(radius,distance):
    """
    this converts a specific intensity measurement to a flux assuming 
    the source is a sphere of a specified radius at a specified
    distance (in the same units)
    """
    ratio=radius/distance
    return pi*ratio*ratio

def flux_to_intensity(radius,distance):
    """
    this converts a flux measurement to a specific intensity assuming 
    the source is a sphere of a specified radius at a specified
    distance (in the same units)
    """
    ratio=radius/distance
    return 1.0/(pi*ratio*ratio)

##Bell&DeJong & Blanton 03
#_band_to_msun={'B':5.47,
#               'V':4.82,
#               'R':4.46,
#               'I':4.14,
#               'K':3.33,
#               'u':6.80,
#               'g':5.45,
#               'r':4.76,
#               'i':4.58,
#               'z':4.51}

#B&M and http://www.ucolick.org/~cnaw/sun.html
_band_to_msun={'U':5.61,
               'B':5.48,
               'V':4.83,
               'R':4.42,
               'I':4.08,
               'J':3.64,
               'H':3.32,
               'K':3.28,
               'u':6.75,
               'g':5.33,
               'r':4.67,
               'i':4.48,
               'z':4.42}
               
def mag_to_lum(M,Mzpt=4.83,Lzpt=1,Merr=None):
    """
    calculate a luminosity from a magnitude
    
    input M can be either a magnitude value, a sequence/array of magnitudes, or
    a dictionary where the keys will be interpreted as specifying the bands  
    (e.g. they can be 'U','B','V', etc. to be looked up as described 
    below or a value for the zero-point)
    
    if Merr is given, will return (L,dL)
    
    Mzpt specifies the magnitude that matches Lzpt, so Lzpt is a unit conversion
    factor for luminosity - e.g. 4.64e32 for ergs in V with solar zero points,
    or 1 for solar
    
    Mzpt can also be 'U','B','V','R','I','J','H', or 'K' and will use B&M values 
    for solar magnitudes or 'u','g','r','i', or 'z' from http://www.ucolick.org/~cnaw/sun.html
    """
    from operator import isMappingType,isSequenceType
    
    dictin = isMappingType(M) and not np.isscalar(M) #numpy scalars are mapping types
    if dictin:
        dkeys = Mzpt = M.keys()
        M = M.values()
    elif type(M) is not np.ndarray:
        M=np.array(M)
        
    if isinstance(Mzpt,basestring):
        Mzpt=_band_to_msun[Mzpt]
    elif isSequenceType(Mzpt):    
        Mzpt=np.array(map(lambda x:_band_to_msun.get(x,x),Mzpt))
        
    #L=(10**((Mzpt-M)/2.5))*Lzpt
    L = _mag_to_flux(M-Mzpt)*Lzpt
    
    if dictin:
        L = dict([t for t in zip(dkeys,L)])
    
    if np.any(Merr):
        #dL=Merr*L/1.0857362047581294 #-1.0857362047581294 = -2.5/ln(10)
        dL = _magerr_to_fluxerr(Merr,M)
        return L,dL
    else:
        return L

def lum_to_mag(L,Mzpt=4.83,Lzpt=1,Lerr=None):
    """
    calculate a magnitude from a luminosity
    
    see mag_to_lum() for syntax details
    """
    from operator import isMappingType,isSequenceType
        
    dictin = isMappingType(L) and not np.isscalar(M) #numpy scalars are mapping types    
    if dictin:
        dkeys = Mzpt = L.keys()
        L = L.values()
    elif type(L) is not np.ndarray:
        L=np.array(L)     
        
    if isinstance(Mzpt,basestring):
        Mzpt=_band_to_msun[Mzpt]
    elif isSequenceType(Mzpt):    
        Mzpt=np.array(map(lambda x:_band_to_msun.get(x,x),Mzpt))
        
    #M=Mzpt-2.5*np.log10(L/Lzpt)
    M = Mzpt+_flux_to_mag(L/Lzpt)

    if dictin:
        M = dict([t for t in zip(dkeys,M)])

    if np.any(Lerr):
        #dM=1.0857362047581294*Lerr/L #-1.0857362047581294 = -2.5/ln(10)
        dM = _fluxerr_to_magerr(Lerr,L)
        return M,dM
    else:
        return M
    
def color_to_flux_ratio(color,band1,band2):
    """
    transforms a color as band1mag-band2mag into a flux ratio flux1/flux2
    """
    return mag_to_lum(color,band1)/mag_to_lum(0,band2)

def flux_ratio_to_color(ratio,band1,band2):
    """
    transforms a flux ratio band1flux/band2flux to a color band1mag-band2mag
    """
    return lum_to_mag(ratio,band1) - lum_to_mag(1,band2)
    
def intensities_to_sig(Is,In,exposure=1,area=1):
    """
    converts photon count intensities (i.e. photon cm^-2 sr^-1) to signficance
    values from poisson statistics.
    
    Is is intensity of the signal
    In is intensity of the background
    """
    return exposure*area*Is*(exposure*area*(Is+In))**-0.5
    
def cosmo_surface_brightness_correction(Sobs,z,mag=True):
    """
    computes
    
    mag determines if mag/as^2 is used or if False, flux
    """
    if mag:
        return Sobs-10*np.log10(1+z)
    else:
        return Sobs*(1+z)**4
    
def kcorrect(mags,zs,magerr=None,filterlist=['U','B','V','R','I']):
    """
    Uses the Blanton et al. 2003 k-correction
    
    requires pidly (http://astronomy.sussex.ac.uk/~anthonys/pidly/) and IDL with
    kcorrect installed
    
    input magnitudes should be of dimension (nfilter,nobj), as should magerr
    zs should be sequence of length nobj
    if magerr is None, 
    
    returns absmag,kcorrections,chi2s
    """
    from .constants import H0
    import pidly
    idl=pidly.IDL()
    idl('.com kcorrect')
    #TODO: figure out if it worked
    
    mags = np.array(mags,copy=False)
    zs = np.array(zs,copy=False).ravel()
    if magerr is None:
        magerr = np.ones_like(mags)
    else:
        magerr = np.array(magerr,copy=False)
    
    if mags.shape[0] != len(filterlist) or magerr.shape[0] != len(filterlist):
        raise ValueError("number of filters and magnitude shapes don't match")
    if mags.shape[1] != zs.size or magerr.shape[1] != zs.size:
        raise ValueError("number of redshifts doesn't match magnitude shapes")
        
    fd = {'U':'bessell_U.par',
          'B':'bessell_B.par',
          'V':'bessell_V.par',
          'R':'bessell_R.par',
          'I':'bessell_I.par',
          'u':'sdss_u0.par',
          'g':'sdss_g0.par',
          'r':'sdss_r0.par',
          'i':'sdss_i0.par',
          'z':'sdss_z0.par'}    
    
    filterlist = [fd.get(fl,fl) for fl in filterlist]
    
    try:
        idl.mags = mags
        idl.magerr = magerr
        idl.zs = zs
        idl.flst = filterlist
        idl('kcorrect,mags,magerr,zs,kcorr,absmag=absmag,chi2=chi2,filterlist=flst,/magnitude,/stddev')
        #idl.pro('kcorrect',mags,magerr,zs,'kcorr',absmag='absmag',chi2='chi2',filterlist=filterlist,magnitude=True,stddev=True,)
        
        kcorr = idl.kcorr
        absmag = idl.absmag +5*np.log10(H0/100)
        chi2 = idl.chi2    
        
        return absmag,kcorr,chi2
        
    finally:
        idl.close()
        
def color_to_Teff(c,colorbands='g-r'):
    """
    estimate color to effective temperature - currently only g-r is available
    from McGurk 10 or B-V from Sekiguchi & Fukugita 00
    """    
    if colorbands == 'g-r':
        if not (-.2 <= c <= .9):
            raise ValueError('g-r can only be between -.2 and .9')
        logteff = np.polyval([.0283,.0488,-.316,3.882],c)
        return 10**logteff
    elif colorbands == 'B-V':
        teffsf = [7078,6621,6214,5852,5532,5248,4996,4770,4567,4382,4209,4044,3882]
        BmVsf = [.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5]
        return np.interp(c,BmVsf,teffsf)
    else:
        raise ValueError('unrecognized color')
    
def teff_to_color(teff,colorbands='g-r'):
    """
    estimate color to effective temperature - currently only g-r is available
    from McGurk 10 or B-V from Flower 96
    """    
    if colorbands == 'g-r':
        .9, -.2
        ps = np.poly1d([.0283,.0488,-.316,3.882-np.log10(teff)])
        roots = ps.r
        goodroot = None
        for r in roots:
            if r.imag==0 and -.2 <= r.real <= .9:
                if goodroot is not None:
                    raise ValueError('multiple matches found - invalid teff')
                goodroot = r.real
        return goodroot
    elif colorbands == 'B-V':
        teffsf = [7078,6621,6214,5852,5532,5248,4996,4770,4567,4382,4209,4044,3882]
        BmVsf = [.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5]
        return np.interp(teff,teffsf[::-1],BmVsf[::-1])
    else:
        raise ValueError('unrecognized color')
    
#<---------------------Load built-in data-------------------------------------->
#TODO: check UBVRI/ugriz S function - energy or quantal?

def __load_UBVRI():
    from .io import _get_package_data
    bandlines = _get_package_data('UBVRIbands.dat').split('\n')
    
    src = bandlines.pop(0).replace('#Source:','').strip()
    bandlines.pop(0)
    
    d={}
    for ln in bandlines:
        if ln.strip() == '':
            pass
        elif 'Band' in ln:
            k = ln.replace('Band','').strip()
            xl = []
            Rl = []
            d[k]=(xl,Rl)
        else:
            t = ln.split()
            xl.append(t[0])
            Rl.append(t[1])
            
    d = dict([(k,ArrayBand(np.array(v[0],dtype='f8'),np.array(v[1],dtype='f8'),name=k)) for k,v in d.iteritems()])
    for v in d.itervalues():
        v.source = src
    return d

def __load_JHK():
    from .io import _get_package_data
    bandlines = _get_package_data('JHKbands.dat').split('\n')
    
    src = bandlines.pop(0).replace('#2MASS J,H, and K bands:','').strip()
    
    d={}
    for ln in bandlines:
        if ln.strip() == '':
            pass
        elif 'Band' in ln:
            k = ln.replace('Band','').strip()
            xl = []
            Rl = []
            d[k]=(xl,Rl)
        else:
            t = ln.split()
            xl.append(t[0])
            Rl.append(t[1])
            
    d = dict([(k,ArrayBand(np.array(v[0],dtype='f8'),np.array(v[1],dtype='f8'),name=k,unit='microns')) for k,v in d.iteritems()])
    for v in d.itervalues():
        v.source = src
    return d

def __load_ugriz():
    from .io import _get_package_data
    bandlines = _get_package_data('ugrizbands.dat').split('\n')
    
    bandlines.pop(0)
    psrc = bandlines.pop(0).replace('#primes:','').strip()
    src = bandlines.pop(0).replace('#SDSS ugriz:','').strip()
    bandlines.pop(0)
    
    d={}
    for ln in bandlines:
        if ln.strip() == '':
            pass
        elif 'Band' in ln:
            k = ln.replace('Band','').strip()
            xl = []
            Rl = []
            d[k]=(xl,Rl)
        else:
            t = ln.split()
            xl.append(t[0])
            Rl.append(t[1])
            
    d = dict([(k,ArrayBand(np.array(v[0],dtype='f8'),np.array(v[1],dtype='f8'),name=k)) for k,v in d.iteritems()])
    for k,v in d.iteritems():
        if "'" in k:
            v.source = psrc
        else:
            v.source = src
    dp = {}
    for k in d.keys():
        if "'" in k: 
            dp[k]=d.pop(k)
    return d,dp

def __load_washington():
    from .io import _get_package_data
    bandlines = _get_package_data('washingtonbands.dat').split('\n')
    
    src = bandlines.pop(0).replace('#2MASS J,H, and K bands:','').strip()
    
    d={}
    for ln in bandlines:
        if ln.strip() == '':
            pass
        elif 'Band' in ln:
            k = ln.replace('Band','').strip()
            xl = []
            Rl = []
            d[k]=(xl,Rl)
        else:
            t = ln.split()
            xl.append(t[0])
            Rl.append(t[1])
            
    d = dict([(k,ArrayBand(np.array(v[0],dtype='f8'),np.array(v[1],dtype='f8'),name=k,unit='nm')) for k,v in d.iteritems()])
    for v in d.itervalues():
        v.source = src
    return d

def __load_human_eye():
    from .io import _get_package_data
    bandlines = _get_package_data('eyeresponse.dat').split('\n')
    
    src = bandlines.pop(0).replace('#from','').strip()
    d={'cone_s':[],'cone_m':[],'cone_l':[]}
    for ln in bandlines:
        if ln.strip() == '':
            pass
        else:
            t = ln.strip().split(',')
            d['cone_l'].append((t[0],t[1]))
            d['cone_m'].append((t[0],t[2]))
            if t[3]!='':
                d['cone_s'].append((t[0],t[3]))
                
    
    d = dict([(k,np.array(v,dtype='f8')) for k,v in d.iteritems()])
    #logarithmic response data - must do 10**data, also, convert x-axis to angstroms
    d = dict([(k,ArrayBand(10*v[:,0],10**v[:,1],name=k)) for k,v in d.iteritems()])
    for v in d.itervalues():
        v.source = src
    return d


#register all the built-in bands
bands = DataObjectRegistry('bands',Band)
str_to_bands = bands.getObjects

bands.register(__load_human_eye(),'eye')
d,dp = __load_ugriz()
bands.register(d,'ugriz')
bands.register(d,'SDSS')
bands.register(dp,'ugriz_prime')
del d,dp
bands.register(__load_UBVRI(),['UBVRI','UBVRcIc'])
bands.register(__load_JHK(),'JHK')
bands.register(__load_washington(),'washington')

#default all to AB mags
set_zeropoint_system('AB','all')

#set the Spectrum eye sensitivity assuming a T=5800 blackbody
from .spec import Spectrum
def __generateRGBSensTuple():
        from .models import BlackbodyModel
        from .phot import bands

        x=np.linspace(3000,9000,1024)
        bm = BlackbodyModel(T=5800)
        spec = Spectrum(x,bm(x))
        eyefluxes = np.array(spec.rgbEyeColor())
        return tuple(eyefluxes)
Spectrum._rgbsensitivity = __generateRGBSensTuple()


class _BwlAdapter(dict): 
    def __getitem__(self,key):
        if key not in self:
            return bands[key].cen
        return dict.__getitem__(self,key)
bandwl = _BwlAdapter()
#photometric band centers - B&M ... deprecated, use bands[band].cen instead
bandwl.update({'U':3650,'B':4450,'V':5510,'R':6580,'I':8060,'u':3520,'g':4800,'r':6250,'i':7690,'z':9110})

del ABCMeta,abstractmethod,abstractproperty,Spectrum,pi,division #clean up namespace