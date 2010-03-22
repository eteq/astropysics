#Copyright (c)2008 Erik Tollerud (etolleru@uci.edu) 
"""

====
spec
====

The :mod:`spec` module contains classes and funtions focused on plotting and
analysis of arbitrary spectra and SEDs, as well as related utility functions.

Defaults tend to be oriented towards optical, but everything should still be
valid in other bands.

.. todo:: examples/tutorials


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.spec
   :parts: 1

Module API
----------

"""

#TODO: RV corrections for heliocentric and LSR 
#TODO: better line-list storage techniques
#TODO: more vetting of zfind
#TODO: reduction framework that accounts for errors with ccd noise models

from __future__ import division,with_statement
from .constants import pi
import numpy as np
try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type

#Spectrum related io module functions
from .io import load_deimos_spectrum,load_all_deimos_spectra,load_wcs_spectrum

class HasSpecUnits(object):
    """
    This class is a mixin superclass for objects that have spectral-like units
    (e.g. flux spectral density :math:`f_\\lambda` or :math:`f_\\nu`, luminosity
    spectral density, etc.). It adds a property :attr:`unit` that can be set
    with a string to change the units of the object.
    
    The following steps be done by subclasses:
    
    1. override method :meth:`HasSpecUnits._applyUnits` - see method doc for
       details.
    2. call :meth:`HasSpecUnits.__init__`  with the `unit` argument providing
       the units of the initial input data.
    """
    __metaclass__ = ABCMeta
    
    def __init__(self,unit):
        self._phystype,self._unit,self._xscaling = self.strToUnit(unit) 
    
    @abstractmethod
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        """
        This method must be overridden by subclasses.  The overriding method
        should use the xftrans OR xfinplace (not both) to convert
        the old x-axis and flux axis into new x axes and flux values
        
        xtrans and xitrans are of the form newx = xtrans(oldx)
        and oldx = xitrans(newx)
        
        xftrans is a function with the signature newx,newf = xftrans(oldx,oldf)
        x is expected to be 1D, and f is expected to be either 1D or must have
        last dimension equal to len(x)
        
        xfinplace(x,f) returns None but uses in-place operators to 
        adjust the values of x and f
        """
        raise NotImplementedError
    
    @staticmethod
    def strToUnit(typestr):
        """
        converts a unit string into a standardized form
        
        returns phystype,unit,scaling
        (scaling relative to cgs)
        """
        if not isinstance(typestr,basestring):
            raise TypeError('invalid type provided as unit')
        
        u = typestr.lower()
        if u == 'wl' or u == 'lambda' or u == 'ang' or u == 'angstroms' or u == 'wavelength' or u == 'wavelength-angstrom':
            val =  'wavelength-angstrom'
            scaling=1e-8
        elif u == 'nm' or u == 'wavelength-nm':
            val =  'wavelength-nm'
            scaling=1e-7
        elif u == 'microns' or u == 'um' or u == 'wavelength-micron':
            val =  'wavelength-micron'
            scaling=1e-4
        elif u == 'm' or u == 'wavelength-m':
            val =  'wavelength-m'
            scaling=1e2
        elif u == 'cm' or u == 'wavelength-cm':
            val =  'wavelength-cm'
            scaling=1
        elif u == 'f' or u == 'nu' or u == 'hz' or u == 'frequency' or u == 'frequency-hz':
            val =  'frequency-Hz'
            scaling=1
        elif u == 'thz' or u == 'frequency-THz':
            val =  'frequency-THz'
            scaling=1e12
        elif u == 'e' or u == 'en' or u == 'energy' or u == 'energy-eV':
            from .constants import ergperev
            val =  'energy-eV'
            scaling=ergperev
        elif u == 'erg' or u == 'energy-erg':
            val =  'energy-erg'  
            scaling=1
        elif u == 'J' or u == 'energy-J':
            val =  'energy-J'   
            scaling=1e-7
        else:
            raise ValueError('unrecognized unit')
        
        phystype,unit = val.split('-')
        return phystype,unit,scaling
        
    def _convertType(self,oldtype,newtype,oldscale,newscale):
        from .constants import c,h
        
        x,flux = self._unitGetX()*oldscale,self._unitGetY()*oldscale # convert to cgs
    
        if newtype == 'energy': #TODO:check
            if oldtype == 'energy': 
                fluxscale = 1
            elif oldtype == 'wavelength':
                x = h*c/x
                fluxscale = x*x/c/h #convert to fnu and divide by h
                
            elif oldtype == 'frequency':
                x = h*x
                fluxscale = 1/h
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'wavelength':
            if oldtype == 'energy': #TODO:check
                x = h*c/x
                fluxscale = x*x*h/c #convert to fnu and then to  by h
            elif oldtype == 'wavelength':
                fluxscale = 1
            elif oldtype == 'frequency':
                x = c/x
                fluxscale = x*x/c #flambda = fnu*nu^2/c
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'frequency':
            if oldtype == 'energy': #TODO:check
                x = x/h
                fluxscale = h
            elif oldtype == 'wavelength':
                x = c/x
                fluxscale = x*x/c #fnu = flambda*lambda^2/c
            elif oldtype == 'frequency':
                fluxscale = 1
            else:
                raise ValueError('unrecognized oldtype')
        else:
            raise ValueError('unrecognized newtype')
        
        return x/newscale,flux*fluxscale/newscale
    
    def __convertUnitType(self,oldx,oldf,oldtype,newtype):
        """
        takes in x and f and returns a new x and f
        """
        from .constants import h,c
        
        if newtype == 'wavelength':
            if oldtype == 'wavelength':
                newf = oldf
                newx = oldx
            elif oldtype == 'frequency':
                newf = oldx*oldx*oldf/c
                newx = c/oldx
            elif oldtype == 'energy':
                newf = oldx*oldx*oldf/h/c
                newx = h*c/oldx
            else:
                raise RuntimeError('Invalid newtype - this should be an unreachable code point')
        elif newtype == 'frequency':
            if oldtype == 'wavelength':
                newf = oldx*oldx*oldf/c
                newx = c/oldx
            elif oldtype == 'frequency':
                newf = oldf
                newx = oldx
            elif oldtype == 'energy':
                newf = h*oldf
                newx = oldx/h
            else:
                raise RuntimeError('Invalid newtype - this should be an unreachable code point')
        elif newtype == 'energy':
            if oldtype == 'wavelength':
                newf = oldx*oldx*oldf/h/c
                newx = h*c/oldx
            elif oldtype == 'frequency':
                newf = oldf/h
                newx = h*oldx
            elif oldtype == 'energy':
                newf = oldf
                newx = oldx
            else:
                raise RuntimeError('Invalid newtype - this should be an unreachable code point')
        else:
            raise RuntimeError('Invalid newtype - this should be an unreachable code point')
        
        return newx,newf
    
    def __xtrans(self,x):
        newx,newf = self.__convertUnitType(x*self.__oldscale,0,self.__oldtype,self.__newtype)
        return newx/self.__newscale
    
    def __xitrans(self,x):
        newx,newf = self.__convertUnitType(x*self.__newscale,0,self.__newtype,self.__oldtype)
        return newx/self.__oldscale
    
    def __xftrans(self,x,f):
        newx,newf = self.__convertUnitType(x*self.__oldscale,f/self.__oldscale,self.__oldtype,self.__newtype)
        #TODO: figure out why in-place doesnt work
        #newx *= self.__rescalefactor
        #newf /= self.__rescalefactor
        return newx/self.__newscale,newf*self.__newscale
    
    def __xfinplace(self,x,f):
        newx,newf = self.__convertUnitType(x*self.__oldscale,f/self.__oldscale,self.__oldtype,self.__newtype)
        #newx *= self.__rescalefactor
        #newf /= self.__rescalefactor
        x[:] = newx/self.__newscale
        f[:] = newf*self.__newscale
    
    def _getUnit(self):
        return self._phystype+'-'+self._unit
    def _setUnit(self,typestr):
        newtype,newunit,newscaling = self.strToUnit(typestr)
        if not (newunit == self._unit and newscaling == self._xscaling):
            
            oldsettings = self._xscaling,self._phystype,self._unit
            self.__newscale = newscaling
            self.__oldscale = self._xscaling
            self.__newtype = newtype
            self.__oldtype = self._phystype
            
            self._phystype = newtype
            self._xscaling = newscaling
            self._unit = newunit
                
            try:    
                self._applyUnits(self.__xtrans,self.__xitrans,self.__xftrans,self.__xfinplace)
            except:
                self._xscaling,self._phystype,self._unit = oldsettings
                raise
            
    unit = property(_getUnit,_setUnit)

class Spectrum(HasSpecUnits):
    """
    Represents a flux/luminosity/intensity
    
    Note that operations are performed in-place, and properties retrieve the 
    same versions that are changed (except ivar)
    """
    def __init__(self,x,flux,err=None,ivar=None,unit='wl',name='',copy=True,sort=True):
        """
        sets the x-axis values and the flux.  Optionally, an error can be
        given by supplying an error or an inverse variance (bot not both)
        
        copy determines if the inputs will be copied if they are already arrays
        
        sort determines if the array values will be sorted by x
        """
        
        x = np.array(x,copy=copy)
        flux = np.array(flux,copy=copy)
        if x.shape != flux.shape:
            raise ValueError("x and flux don't match shapes")
        
        if ivar is None and err is None:
            err = np.zeros_like(flux)
        elif ivar is not None:
            if np.isscalar(ivar):
                err = ivar**-0.5*np.ones_like(flux)
            else:
                err = np.array(ivar,copy=False)**-0.5
            if err.shape != flux.shape:
                raise ValueError("ivar and flux don't match shapes")
            
        elif err is not None:
            if np.isscalar(err):
                err = err*np.ones_like(flux)
            else:
                err=np.array(err,copy=copy)
                err[err<0]=-1*err[err<0]
            if err.shape != flux.shape:
                raise ValueError("err and flux don't match shapes")
            
        else:
            raise ValueError("can't set both err and ivar at the same time")
        
        HasSpecUnits.__init__(self,unit)
        
        if sort:
            sorti = np.argsort(x)
            x = x[sorti]
            flux = flux[sorti]
            err = err[sorti]
        
        self._x = x
        self._flux = flux
        self._err = err
        
        self.continuum = None
        
        self.name = name
        self.z = 0 #redshift
        self._zqual = -1 #redshift quality
        
        self._features = []
    
    def __getstate__(self):
        #state = super(HasSpecUnits,self).__getstate__()
        state = self.__dict__
        #necessary because spylot sometimes replaces this with a feature
        if not type(state['_features']) is list:
            state = dict(state) #make a new dictionary so as not to override the current one
            state['_features'] = list(state['_features'])
        return state
    
    def save(self,fn,**kwargs):
        """
        Save this Spectrum as a pickled object with the provided file name.  
        kwargs are passed into utils.fpickle
        """
        from .utils import fpickle
        fpickle(self,fn,**kwargs)
    
    @staticmethod
    def load(fn):
        """
        Load a saved Spectrum from the given file
        """
        from .utils import funpickle
        obj = funpickle(fn,0)
        if obj.__class__.__name__ != 'Spectrum':
            raise TypeError('file does not contain a Spectrum')
        return obj
    
    _zqmap = {-1:'unknown',0:'none',1:'bad',2:'average',3:'good',4:'excellent'}
    _zqmapi = dict([(v,k) for k,v in _zqmap.iteritems()])
    def _getZqual(self):
        return self._zqual
    def _setZqual(self,val):
        if isinstance(val,basestring):
            self.zqual = self._zqmapi[val]
        else:
            val = int(val)
            if val in self._zqmap.keys():
                self._zqual = val
            else:
                raise KeyError('invalid zqual value')
    zqual = property(_getZqual,_setZqual,doc=None)
    
    def getZQualStr(self):
        return self._zqmap[self._zqual]
        
    def copy(self):
        """
        Generates a deep copy of this Spectrum
        """
        from copy import deepcopy
        return deepcopy(self)
    
    #units support
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        if hasattr(self,'_contop'):
            raise ValueError('continuum operation applied - revert before changing units')
        if self.continuum is None:
            pass
        elif callable(self.continuum):
            self.continuum = None
        else:
            self.continuum = xftrans(self._x,self.continuum)[1]
        
        err = xftrans(self._x,self._err)[1]
        #x,flux = xftrans(self._x,self._flux)
        xfinplace(self._x,self._flux)
        self._err[:] = err
        #self._x[:] = x
        #self._flux[:] = flux
    
    
    #------------------------Properties--------------------------------->
    @property
    def npix(self):
        return self._x.shape[-1]
    
    @property
    def shape(self):
        return self._x.shape
    
    def _getFlux(self):
        return self._flux
    def _setFlux(self,flux):
        flux = np.array(flux,copy=False)
        if flux.shape != self._flux.shape:
            raise ValueError("new flux doesn't match old flux shape")
        self._flux[:] = np.array(flux)
    flux = property(_getFlux,_setFlux)
    
    def _getNFlux(self):
        from .constants import h
        
        oldtype = self.unit
        try:
            self.unit = 'hz'
            return self._flux/h/self._x
        finally:
            self.type = oldtype
    def _setNFlux(self,nflux):
        nflux = np.array(nflux)
        if nflux.shape != self._flux.shape:
            raise ValueError("new flux doesn't match flux shape")
        raise NotImplementedError
        self._flux = 42
    nflux = property(_getNFlux,_setNFlux)
    
    def _getErr(self):
        return self._err
    def _setErr(self,err):
        if err.shape != self._err.shape:
            raise ValueError("new err doesn't match old err shape")
        self._err[:]=np.array(err)
    err = property(_getErr,_setErr)
    
    def _getIvar(self):
        return 1/self._err/self._err
    def _setIvar(self,ivar):
        if ivar.shape != self._flux.shape:
            raise ValueError("newivar doesn't match flux shape")
        self._err=np.array(ivar)**-0.5
    ivar = property(_getIvar,_setIvar)
    
    def _getX(self):
        return self._x
    def _setX(self,x):
        x = np.array(x)
        if x.shape != self._flux.shape:
            raise ValueError("new x doesn't match flux shape")
        self._x = np.array(x)
    x = property(_getX,_setX,doc='x-axis as measured')
    
    def _getX0(self):
        return self._x/(self.z+1)
    def _setX0(self,x):
        self.x = x
        self._x*=(self.z+1)
    x0 = property(_getX0,_setX0,doc='x-axis in rest frame')
    
    def _getFeatures(self):
        return tuple(self._features)
    def _setFeatures(self,val):
        for v in val:
            if not isinstance(v,SpectralFeature):
                raise TypeError('features must all be SpectralFeature objects')
        for i in range(len(self._features)):
            del self.features[0]
        self._features.extend(val)
    features = property(_getFeatures,_setFeatures,doc='The spectral features in this spectrum')
    
    
    
    #<----------------------Tests/Info--------------------------------->
    def __len__(self):
        return len(self._x)
    
    def __eq__(self,other):
        try:
            xeq = self._x == other._x
            feq = self._flux == other._flux
            eeq = self._err == other._err
            return np.all(xeq & feq & eeq) 
        except AttributeError:
            return False
    
    def getUnitFlux(self,units,err=False):
        """
        returns x and flux of this spectrum in a new unit system without 
        changing the selected unit
        
        err can be False, True or 'ivar'
        
        if err is False, returns x,flux
        if err is True, returns x,flux,err
        if err is 'ivar', returns x,flux,ivar
        """
        oldunit = self.unit
        try:
            self.unit = units
            if err == 'ivar':
                res = self.x.copy(),self.flux.copy(),self.ivar.copy()
            elif err:
                res = self.x.copy(),self.flux.copy(),self.err.copy()
            else:
                res = self.x.copy(),self.flux.copy()
        finally:
            self.unit = oldunit
        
        return res
    
    def getPhotonFlux(self):
        """
        returns flux adjusted to represent photons s^-1 cm^-2 xunit^-1
        """
        from .constants import h,c
        x = self.x*self._xscaling
        if 'wavelength' in self.unit:
            factor = x/h/c
        elif 'frequency' in self.unit:
            factor = 1/h/x
        elif 'energy' in self.unit:
            factor = 1/x
        else:
            raise ValueError('Unrecognized unit for photon conversion')
        
        return factor*self.flux
            
    def getDx(self,mean=True):
        """
        get the spacing of the x-axis, which is always 1 element shorter than x
        
        if mean, returns the mean of the spacing
        """
        dx = np.convolve(self._x,(1,-1),mode='valid')
        if mean:
            return dx.mean()
        else:
            return dx
        
    def getDlogx(self,mean=True,logbase=10):
        """
        get the logarithmic spacing of the x-axis, which is always 1 element
        shorter than x
        
        if mean, returns the mean of the spacing
        """
        x = np.log(self._x)/np.log(logbase)
        dlogx = np.convolve(x,(1,-1),mode='valid')
        if mean:
            return dlogx.mean()
        else:
            return dlogx
        
    def isXMatched(self,other,tol=1e-10):
        """
        tests if the x-axis of this Spectrum matches that of another Spectrum
        or equal length array, with an average deviation less than tol
        """
        from operator import isSequenceType
        if isSequenceType(other):
            ox = other
        else:
            ox = other.x
            
        try:
            return np.std(self.x - ox) < tol
        except (TypeError,ValueError):
            return False
    
    def isLinear(self,eps=1e-10):
        """
        Determines if the x-spacing is linear (e.g. evenly spaced so that 
        dx ~ constant with a standard deviation of eps)
        """
        return np.std(self.getDx(False)) < eps
        
    def isLogarithmic(self,eps=1e-10):
        """
        Determines if the x-spacing is logarithmic (e.g. evenly spaced in 
        logarithmic bins so that dx ~ x to tolerance of eps)
        """
        return np.std(self.getDlogx(False)) < eps
    
    #<----------------------Operations---------------------------->
    
    def smooth(self,width=1,filtertype='gaussian',replace=True):
        """
        smooths the flux in this object by a filter of the given filtertype 
        (can be either 'gaussian' or 'boxcar'/'uniform')
        
        if replace is True, the flux in this object is replaced by the smoothed
        flux and the error is smoothed in the same fashion
        
        width is in pixels, either as sigma for gaussian or half-width 
        for boxcar
        
        returns smoothedflux,smoothederr
        """
        import scipy.ndimage as ndi 
        
        if filtertype is None:
            if width > 0:
                fitertype = 'gaussian'
            else:
                filtertype = 'boxcar'
                width = -1*width
        
        if filtertype == 'gaussian':
            filter = ndi.gaussian_filter1d
        elif filtertype == 'boxcar' or type == 'uniform':
            filter = ndi.uniform_filter1d
            width = 2*width
        else:
            raise ValueError('unrecognized filter type %s'%filtertype)
        
        smoothedflux = filter(self._flux,width)
        smoothederr = filter(self._err,width)
        
        if replace:
            self.flux = smoothedflux
            self.err = smoothederr
        
        return smoothedflux,smoothederr
    
    def resample(self,newx,interpolation='linear',replace=True,**kwargs):
        """
        this interpolates the flux to populate the supplied x-axis
        
        kwargs go into the interpolation routine as described below
        
        interpolations can be:
        'linear': simple linear interpolation
        'spline': a k-order spline with smoothing factor s is used, where s and 
        k are set by kwargs.  if the 'save' kwarg is True, the spline is saved
        and will be used for subsequent resamplings.  if 'clear' is True, the 
        existing spline will be cleared and a new spline will be calculated.
        note that default spline has smoothing=0, which interpolates through
        every point
        
        WARNING: this does not treat the errors properly yet - currently just interpolates
        
        returns newx,newflux,newerr
        """
        if interpolation == 'linear':
            newflux = np.interp(newx,self._x,self._flux)
            #TODO: fix errors
            newerr = np.interp(newx,self._x,self._err)
        elif 'spline' in interpolation:
            from scipy.interpolate import UnivariateSpline
            
            s = kwargs.pop('s',0)
            k = kwargs.pop('k',3)
            save = kwargs.pop('save',False)
            clear = kwargs.pop('clear',False)
            
            fspline = espline = None
            if hasattr(self,'_spline') and not (clear or save):
                fspline,espline = self._spline
                
                    
            if fspline is None:
                fspline = UnivariateSpline(self._x,self._flux,k=k,s=s)
            if espline is None:
                espline = UnivariateSpline(self._x,self._flux,k=k,s=s)
                
            if save:
                self._spline = (fspline,espline)
            
            newflux = fspline(newx)
            #TODO: fix errors
            newerr = espline(newx)
        else:
            raise ValueError('unrecognized interpolation technique')
        
        if len(kwargs)>0:
            raise TypeError('unexpected keyword %s'%kwargs.keys()[0])
        
        if replace:
            if newx.shape != self._x.shape:
                raise ValueError("can't replace if new x-axis shape doesn't match old")
            
            self._x[:] = newx 
            self._flux[:] = newflux
            self._err[:] = newerr   
        return newx,newflux,newerr
    
    def linearize(self,lower=None,upper=None,**kwargs):
        """
        convinience function for resampling to an equally-spaced linear x-axis
        
        if lower or upper are None, the upper and lower x values are used
        kwargs go into Spectrum.resample
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        
        newx = np.linspace(lower,upper,self.npix)
        return self.resample(newx,**kwargs)
    
    def logify(self,lower=None,upper=None,**kwargs):
        """
        convinience function for resampling to an x-axis that is evenly spaced
        in logarithmic bins.  Note that lower and upper are the x-axis values 
        themselves, NOT log(xvalue)
        
        if lower or upper are None, the upper and lower x values are used
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        
        newx = np.logspace(np.log10(lower),np.log10(upper),self.npix)
        return self.resample(newx,**kwargs)
    
    def computeMag(self,bands,**kwargs):
        """
        this computes the magnitude of the spectrum in the provided bands
        
        bands can be a sequence of band objects or strings which will be mapped
        to bands through the phot.bands registry
        
        kwargs are passed into phot.Band.computeMag
        """
        kwargs['__domags']=True
        return self.computeFlux(bands,**kwargs)
        
    
    def computeFlux(self,bands,**kwargs):
        """
        this computes the flux of the spectrum in the provided bands
        
        bands can be a sequence of band objects or strings which will be mapped
        to bands through the phot.bands registry
        
        kwargs are passed into phot.Band.computeFlux
        """
        from operator import isMappingType
        from .phot import str_to_bands
        
#        if isinstance(bands,basestring) or isinstance(bands,phot.Band):
#            bands = [bands]
#            scalarout = True
#        else:
#            scalarout = False
        
#        bands = [phot.bands[b] if isinstance(b,basestring) else b for b in bands]
#        bl = []
#        for b in bands:
#            if isMappingType(b):
#                bl.extend(b.values())
#            else:
#                bl.append(b)
#        bands = bl
        
        scalarout = isinstance(bands,basestring) or isinstance(bands,phot.Band)
        bands = str_to_bands(bands)
        
        if kwargs.pop('__domags',False):
            res = [b.computeMag(self,**kwargs) for b in bands]
        else:
            res = [b.computeFlux(self,**kwargs) for b in bands]
        
        if scalarout and len(res) == 1:
            return res[0]
        else:
            return res
        
        
    _rgbsensitivity = (1,1,1) #this is adjusted to a 5800 K blackbody after the Spectrum class is created
    def rgbEyeColor(self):
        """
        this uses the 'eye' group in phot.bands to convert a spectrum to
        an (r,g,b) tuple
        """
        from .phot import bands
        spec = self
        
        eyed = bands['eye']
        if len(eyed) != 3:
            raise ValueError('eye bands are not length 3')
        eyefluxes = np.array([b.computeFlux(spec) for b in sorted(eyed.values())])
        eyefluxes = eyefluxes[::-1] #b,g,r -> r,g,b
        
        #now normalize by eye sensitivities -- default is computed by assuming
        #a T=5800 blackbody gives (1,1,1)
        eyefluxes /= self._rgbsensitivity
        
        maxe = eyefluxes.max()
        eyefluxes /= maxe
        
        return tuple(eyefluxes)
    
    def fitContinuum(self,model=None,weighted=False,evaluate=False,
                          interactive=False,**kwargs):
        """
        this method computes a continuum fit to the spectrum using a model
        from astropysics.models (list_models will give all options) or
        an callable with a fitData(x,y) function
        
        if model is None, the existing model will be used if present, 
        or if there is None, the default is 'uniformknotspline'.  Otherwise,
        it may be any acceptable model (see `models.get_model`)
        
        kwargs are passed into the constructor for the model
        
        if weighted, the inverse variance will be used as weights to the 
        continuum fit
        
        if interactive is True, the fitgui interface will be displayed to 
        tune the continuum fit
        
        the fitted model is assigned to self.continuum or evaluated at the
        spectrum points if evaluate is True and the results are set to 
        self.continuum
        """
        
        if model is None and self.continuum is None:
            model = 'uniformknotspline'
        
        #for the default, choose a reasonable number of knots
        if model == 'uniformknotspline' and 'nknots' not in kwargs:
            kwargs['nknots'] = 4
        
        if model is None and self.continuum is not None:
            model = self.continuum
            if not interactive:
                model.fitData(self.x,self.flux,weights=(self.ivar if weighted else None))
        else:
            if isinstance(model,basestring):
                from .models import get_model
                model = get_model(model)(**kwargs)
        
            if not (callable(model) and hasattr(model,'fitData')):
                raise ValueError('provided model object cannot fit data')
            
            model.fitData(self.x,self.flux,weights=(self.ivar if weighted else None))
        
        if interactive:
            from .gui.fitgui import FitGui
            
            if interactive == 'reuse' and hasattr(self,'_contfg'):
                fg = self._contfg
            else:
                fg = FitGui(self.x,self.flux,model=model,weights=(self.ivar if weighted else None))
                fg.plot.plots['data'][0].marker = 'dot'
                fg.plot.plots['data'][0].marker_size = 2
                fg.plot.plots['model'][0].line_style = 'solid'
                
            if fg.configure_traits(kind='livemodal'):
                model = fg.tmodel.model
            else:
                model = None
                
            if interactive == 'reuse':
                self._contfg = fg
            elif hasattr(self,'_contfg'):
                del self._contfg
                
        if model is not None:    
            self.continuum = model(self.x) if evaluate else model
        
    def subtractContinuum(self):
        """
        Subtract the continuum from the flux
        """
        if hasattr(self,'_contop'):
            raise ValueError('%s already performed on continuum'%self._contop)
        
        if self.continuum is None:
            raise ValueError('no continuum defined')
        elif callable(self.continuum):
            cont = self.continuum(self.x)
        else:
            cont = self.continuum
            
        self.flux = self.flux - cont
        self._contop = 'subtraction'
            
    def normalizeByContinuum(self):
        """
        Divide by the flux by the continuum
        """
        if hasattr(self,'_contop'):
            raise ValueError('%s already performed on continuum'%self._contop)
        
        if self.continuum is None:
            raise ValueError('no continuum defined')
        elif callable(self.continuum):
            cont = self.continuum(self.x)
        else:
            cont = self.continuum
            
        self.flux = self.flux/cont
        self._contop = 'normalize'
        
    def rejectOutliersFromContinuum(self,sig=3,iters=1,center='median',savecont=False):
        """
        rejects outliers and returns the resulting continuum. see 
        `utils.sigma_clip` for arguments
        
        returns a pair of maksed arrays xmasked,contmasked
        
        if savecont is True, the outlier-rejected value will be saved as
        the new continuum
        """
        from .utils import sigma_clip
        
        if self.continuum is None:
            raise ValueError('no continuum defined')
        elif callable(self.continuum):
            cont = self.continuum(self.x)
        else:
            cont = self.continuum
            
        contma = sigma_clip(cont,sig=sig,iters=iters,center=center,maout='copy')
        xma = np.ma.MaskedArray(self.x,contma.mask,copy=True)
        
        if savecont:
            self.continuum = contma
        
        return xma,contma
            
    def revertContinuum(self):
        """
        Revert to flux before continuum subtraction
        """
        if self.continuum is None:
            raise ValueError('no continuum defined')
        elif callable(self.continuum):
            cont = self.continuum(self.x)
        else:
            cont = self.continuum
            
        if hasattr(self,'_contop'):
            if self._contop == 'subtraction':
                self.flux = self.flux+cont
            elif self._contop == 'normalize':
                self.flux = self.flux*cont
            else:
                raise RuntimeError('invalid continuum operation')
            del self._contop
        else:
            raise ValueError('no continuum action performed')
        
    def addFeatureRange(self,lower,upper,continuum='fromspec',identity=None):
        sf = SpectralFeature(unit=self.unit,extent=(lower,upper),continuum=continuum)
        
        sf.computeFeatureData(self,edges='interp')
        if identity is not None:
            sf.identify(identity)
            
        self._features.append(sf)
    
    def addFeatureLocation(self,loc,smooth=None,window=200,**kwargs):
        """
        add a spectral feature at the provided x-location.  Edges are inferred
        from either the continuum zero-point (if continuum is present) or 
        where the derivative changes sign
        
        If window is specified, it is the number of pixels to check for 
        the edge.
        
        If smooth is not None, it specifies a smoothing size to apply - 
        positive for gaussian sigma or negative for boxcar half-width.
        
        kwargs are passed into `addFeatureRange`
        """
        
        xi = np.interp(loc,self.x,np.arange(self.x.size))
        if window is not None:
            wl,wu = int(np.floor(xi-window/2)),int(np.ceil(xi+window/2))
            xwi = (wu-wl)/2.0
            
            if wl < 0:
                xwi += wl
                wl = 0
                                
            
            m = slice(wl,wu)
        else:
            m = slice(None)
            xwi = (self.x.size-1)/2.0
            
        x,y = self.x[m],self.flux[m]
        if smooth is not None:
            if smooth > 0:
                from scipy.ndimage import gaussian_filter1d as filter
            else:
                from scipy.ndimage import unifom_filter1d as filter
                smooth = -2*smooth
            y = filter(y,smooth)
        
        if self.continuum is None:
            dysg = np.sign(np.convolve(y,[0,-1,1],mode='same'))
            trans = np.convolve(dysg,[-1,1,0],mode='same')
            
            transi = np.where(trans!=0)[0]
            transorti = np.argsort(np.abs(transi-xwi))
            sgn = trans[transorti[0]]
            loweri = transi[transorti][(transi[transorti] < xwi) & (trans[transorti]!=sgn)][0]
            upperi = transi[transorti][(transi[transorti] > xwi) & (trans[transorti]!=sgn)][0]
            
            cont = 'fromspec'
        else:
            if callable(self.continuum):
                cont = self.continuum(x)
            else:
                cont = self.continuum[m]
            y = y - cont
            
            transi = np.where(np.sign(y)!=0)[0]
            transorti = np.argsort(np.abs(transi-xwi))
            loweri = transi[transorti][transi[transorti] < xwi][0]
            upperi = transi[transorti][transi[transorti] > xwi][0]
            
        lower,upper = x[loweri],x[upperi]
        kwargs['continuum'] = cont
        
        return self.addFeatureRange(lower,upper,**kwargs)
    
    def removeFeatureLocation(self,loc):
        """
        removes the feature with a center nearest to the requested location
        
        raises IndexError if there are no features
        """
        cens = np.array([f.center for f in self._features],dtype=float)
        seps = np.abs(cens-loc)
        self.removeSpectralFeature(np.argsort(seps)[0])
            
    
    def removeSpectralFeature(self,iorname):
        """
        remove the requested spectral feature, either by index
        or by line name
        """
        if isinstance(iorname,basestring):
            i = iorname.index([sf.name for sf in self._features])
        elif isinstance(iorname,int):
            i = iorname
        else:
            raise TypeError()
        del self._features[i]
        
    def plot(self,fmt=None,ploterrs=.1,plotcontinuum=True,smoothing=None,
                  step=True,clf=True,colors=('b','g','r','k'),restframe=True,
                  **kwargs):
        """
        Use :mod:`matplotlib` to plot the :class:`Spectrum` object. The
        resulting plot shows the flux, error (if present), and continuum (if
        present).
        
        If `step` is True, the plot will be a step plot instead of a line plot.
        
        `colors` should be a 3-tuple that applies to
        (spectrum,error,invaliderror,continuum) and kwargs go into spectrum and
        error plots.
        
        If `restframe` is True, the x-axis is offset to the rest frame.
        
        If `ploterrs` or `plotcontinuum` is a number, the plot will be scaled so
        that the mean value matches the mean of the spectrum times the numeric
        value. If either are True, the scaling will match the actual value. If
        False, the plots will not be shown.
        
        kwargs are passed into either the :func:`matplotlib.pyplot.plot` or
        :func:`matplotlib.pyplot.step` function.
        """
        
        import matplotlib.pyplot as plt
        
        if step:
            kwargs.setdefault('where','mid')
        
        if smoothing:
            x,(y,e) = self.x0 if restframe else self.x,self.smooth(smoothing,replace=False)
        else:
            x,y,e = self.x0 if restframe else self.x,self.flux,self.err
            
        if len(x)==3:
            dx1 = x[1]-x[0]
            dx2 = x[2]-x[1]
            x = np.array((x[0]-dx1/2,x[0]+dx1/2,x[1]+dx2/2,x[2]+dx2/2))
            y = np.array((y[0],y[0],y[1],y[2]))
            e = np.array((e[0],e[0],e[1],e[2]))
        elif len(x)==2:
            dx = x[1]-x[0]
            x = np.array((x[0]-dx/2,x[0]+dx/2,x[1]+dx/2))
            y = np.array((y[0],y[0],y[1]))
            e = np.array((e[0],e[0],e[1]))
        elif len(x)==1:
            x = np.array((0,2*x[0]))
            y = np.array((y[0],y[0]))
            e = np.array((e[0],e[0]))
            
        if clf:
            plt.clf()
            
        kwargs['c'] = colors[0]
        if fmt is None:
            if step:
                res = [plt.step(x,y,**kwargs)]
            else:
                res = [plt.plot(x,y,**kwargs)]
        else:
            if step:
                res = [plt.step(x,y,fmt,**kwargs)]
            else:
                res = [plt.plot(x,y,fmt,**kwargs)]
            
        if ploterrs and np.any(e):
            from operator import isMappingType
            
            m = (e < np.max(y)*2) & np.isfinite(e)
            
            if ploterrs is True:
                scale = 1
            elif np.isscalar(ploterrs):
                scale = float(ploterrs)*np.mean(y)/np.mean(e[m])
            
            if not isMappingType(ploterrs):
                ploterrs = {}
            ploterrs.setdefault('ls','-')
            
            if np.sum(m) > 0:
                kwargs['c'] = colors[1]
                if step:
                    res.append(plt.step(x[m],scale*e[m],**kwargs))
                else:
                    res.append(plt.plot(x[m],scale*e[m],**kwargs))
            if np.sum(~m) > 0:
                if step:
                    res.append(plt.step(x[~m],scale*np.mean(e[m] if np.sum(m)>0 else y)*np.ones(sum(~m)),'*',mew=0,color=colors[2]))
                else:
                    res.append(plt.plot(x[~m],scale*np.mean(e[m] if np.sum(m)>0 else y)*np.ones(sum(~m)),'*',mew=0,color=colors[2]))
                
        if plotcontinuum and self.continuum is not None:
            if callable(self.continuum):
                cont = self.continuum(self.x)
            else:
                cont = self.continuum
            
            if plotcontinuum is True:
                scale = 1
            elif np.isscalar(plotcontinuum):
                scale = float(plotcontinuum)*np.mean(y)/np.mean(cont)
                
            kwargs['c'] = colors[3]
            kwargs['ls'] =  '--'
            if step:
                res.append(plt.step(self.x,scale*cont,**kwargs))
            else:
                res.append(plt.plot(self.x,scale*cont,**kwargs))
                
                
        plt.xlim(np.min(x),np.max(x))
        
        xl=self.unit
        xl=xl.replace('wavelength','\\lambda')
        xl=xl.replace('frequency','\\nu')
        xl=xl.replace('energy','E')
        xl=xl.replace('angstrom','\\AA')
        xl=xl.replace('micron','\\mu m')
        xl=tuple(xl.split('-'))
        plt.xlabel('$%s/{\\rm %s}$'%xl)
        
        plt.ylabel('$ {\\rm Flux}/({\\rm erg}\\, {\\rm s}^{-1}\\, {\\rm cm}^{-2} {\\rm %s}^{-1})$'%xl[1])
            
        return res
    
    def plotSpylot(self,spylotinstance=None,show=True):
        """
        Displays the spectrum using the astropysics.gui.spylot spectrum plotter.
        
        if `spylotinstance` is None, a new instance of Spylot will be generated
        with this Spectrum as the only Spectrum.  If it is a Spylot instance,
        this Spectrum will be added to the end of the list on that instance.
        
        If `show` is true and `spylotinstance` is None, the GUI will be shown
        using edit_traits (kwargs can be passed in if show is a dictionary),
        and if `spylotinstance` is an instance, the current Spectrum in that 
        instance will be switched to this Spectrum.  Otherwise, no GUI changes
        will occur.
        
        returns the Spylot instance
        """
        from operator import isMappingType
        from .gui.spylot import Spylot
        
        if spylotinstance is None:
            spylotinstance = Spylot([self])
            if show:
                if not isMappingType(show):
                    show = {}
                spylotinstance.edit_traits(**show)
            
        elif isinstance(spylotinstance,Spylot):
            spylotinstance.specs.append(self)
            
            if show:
                spylotinstance.currspeci = len(spylotinstance.specs)-1
            else:
                spylotinstance._specs_changed() #TODO: remove this when spylot properly responds
            
        else:
            raise TypeError('spylotinstance is not a Spylot object or None')
        
        return spylotinstance


class FunctionSpectrum(Spectrum):
    """
    This is a :class:`Spectrum` generated by functional forms rather than using
    fixed samples. Setting x will determine where the function is evaluated (and
    do the evaluation).
    """
    def __init__(self,xi,fluxf,errf=None,ivarf=None,unit='wl'):
        if ivarf is not None and errf is not None:    
            raise ValueError("can't set both err and ivar at the same time")
        
        flux = fluxf(xi)
        self._fluxf = fluxf
        if ivarf is not None:
            errf = lambda x:ivarf(x)**-0.5
        err = errf(xi) if errf is not None else np.zeros_like(xi)
        self._errf = errf
        super(FunctionSpectrum,self).__init__(np.array(xi,copy=True),flux,err,unit=unit)
        
        self._xitranses = []
        self._xftranses = []
        self._unithist = []
        
    def _setX(self,x):
        self._x = x = np.array(x)
        
        for xit in self._xitranses:
            x = xit(x)
            
        flux = self._fluxf(x)
        err = self._errf(x) if self._errf is not None else np.zeros_like(x)
        
        for xft in self._xftranses:
            err = xft(x,err)[1]
            x,flux = xft(x,flux)
            
        self._flux = flux
        self._err = err    
    x = property(Spectrum._getX,_setX)
    
    @property
    def flux(self):
        return self._flux
    
    @property
    def ivar(self):
        return 1/self._err/self._err
    
    @property
    def err(self):
        return self._err
    
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        super(FunctionSpectrum,self)._applyUnits(xtrans,xitrans,xftrans,xfinplace)
        if self.unit in self._unithist:
            i = self._unithist.index(self.unit)+1
            self._xitranses = self._xitranses[:i]
            self._xftranses = self._xftranses[:i]
            self._unithist = self._unithist[:i]
        else:
            self._xitranses.append(xitrans)
            self._xftranses.append(xftrans)
            self._unithist.append(self.unit)
        
#Make a special Spectrum to use for AB Magnitude calibrations        
ABSpec = FunctionSpectrum(np.linspace(5e13,3e15,1024),lambda x:np.ones_like(x)*10**(48.6/-2.5),unit='hz')
        
class ModelSpectrum(FunctionSpectrum):
    """
    A Spectrum that follows the provided `astropysics.models.FunctionModel1D`
    """
    def __init__(self,xi,model,noisetype=None,noisefunc=None,unit='wl'):     
        """
        `noisefunc` should be of the form nf(x,flux)
    
        `noisetype` determines the noise model - can be:
        
        * None
            no noise
        * 'poisson'
            poisson noise where the noisefunc determines the poissonian scale factor
        * 'normal' or 'gauss'
            normally-distributed noise where the noisefunc determines sigma
        """
        super(ModelSpectrum,self).__init__(xi,model,errf=None,ivarf=None,unit=unit)
        self.noisefunc = noisefunc
        self.noisetype = noisetype
        self._mod = model
        
    def _setX(self,x):
        super(ModelSpectrum,self)._setX(x)
        if self._noise is not None:
            self._noise(x,self._flux)
#        if self._noisetype is not None:
#            nf = self._noisefunc(x,self._flux)
#            if self._noisetype == 'normal':
#                offset = np.random.normal(scale=nf,size=x.size)
#                self._flux += offset
#            elif self._noisetype == 'poisson':
#                frac = np.random.poisson(nf,size=x.size)/nf
#                self._flux *= frac
#            else:
#                raise ValueError('internally invalid noise type')
    x = property(Spectrum._getX,_setX)
        
    @property
    def model(self):
        return self._mod
    
    def _getNoisefunc(self):
        return self._noisefunc
    def _setNoisefunc(self,val):
        if val is None:
            self._noisefunc = None
        else:
            if not callable(val):
                if np.isscalar(val):
                    scalarval = val
                    val = lambda x,f:scalarval*np.ones(x.size)
                else:
                    raise ValueError('noisefunc is not a callable or a scalar')
            
            oldfunc = self._noisefunc
            try:
                self._noisefunc = val
                if self._noisetype is not None:
                    self._setX(self.x)
            except:
                self._noisefunc = oldfunc
                raise
                
            
    noisefunc = property(_getNoisefunc,_setNoisefunc)
    
    def _getNoisetype(self):
        return self._noisetype
    def _setNoisetype(self,val):
        if val is None:
            self._noisetype = None
            self._errf = None
            self._noise = None
        else:
            if self._noisefunc is None:
                raise ValueError("can't apply noise model without a noise function")   
            
            if val == 'poisson':
                val = 'poisson'
                newerrf = lambda x:self._flux*self._noisefunc(x,self._flux)**-0.5
                def newnoise(x,flux):
                    nf = self._noisefunc(x,flux)
                    flux *= np.random.poisson(nf,size=x.size)/nf
            elif 'gauss' in val or 'normal' == val:
                val = 'normal'
                newerrf = lambda x:self._noisefunc(x,self._flux)
                def newnoise(x,flux):
                    nf = self._noisefunc(x,flux)
                    flux += np.random.normal(scale=nf,size=x.size)
            else:
                raise ValueError('unrecognized noise %s'%val)
            
            oldtype = self._noisetype
            olderrf = self._errf
            oldnoise = self._noise
            try:
                self._noisetype = val
                self._noise = newnoise
                self._errf = newerrf
                self._setX(self.x)
            except:
                self._noisetype = oldtype
                self._errf = olderrf
                self._noise = oldnoise
                raise
    noisetype = property(_getNoisetype,_setNoisetype)
    
    def fitToSpectrum(self,spec,**kwargs):
        """
        Fit the model parameters to match this Spectrum to the supplied 
        Spectrum
        
        kwargs are passed into the fitData method of the model
        """
        return self._mod.fitData(spec._x,spec._flux,**kwargs)
    def plotModels(self,*args,**kwargs):
        """
        calls the model's plot method with args and kwargs
        """
        return self._mod.plot(*args,**kwargs)

class SpectralFeature(HasSpecUnits):
    """
    This class represents a Spectral Feature/line in a Spectrum.
    
    All SpecralFeatures have a rest value,observed value (possibly with
    error), flux value (possibly with error), and equivalent width (possibly
    with error)
    
    Note that equivalent width is always expected to be in angstroms
    """
    def __init__(self,extent,unit='wavelength',continuum='fromspec'):
        HasSpecUnits.__init__(self,unit)
        
        self.extent = extent
        self.continuum = continuum
        self.rest = 1
        self._center = None
        self.centererr = None
        self.flux = 0
        self.fluxerr = 0
        self.ew = None
        self.ewerr = None
        self.known = None
        self.model = None
    
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        self.rest = xtrans(self.rest)
            
        if self.center is not None:
            oldobs = self.center
            self.center = newobs = xtrans(self.center)
            if self.centererr is not None:
                op = oldobs+self.centererr
                om = oldobs-self.centererr
                op,om = xtrans(op),xtrans(om)
                self.centererr = (op+om)/2
        if self.extent is not None:
            self.extent = (xtrans(self.extent[0]),xtrans(self.extent[1]))
            
    def _getCenter(self):
        if self._center is None:
            if self.extent is None:
                return None
            else:
                return np.sum(self.extent)/len(self.extent) 
        else:
            return self._center
    def _setCenter(self,val):
        self._center = val
    center = property(_getCenter,_setCenter,doc='central wavelength/frrequency/energy of this feature')
    
                
    def computeFeatureData(self,spec,fit=False,interactive=False,edges='interp'):  
        """
        computes the data (center,flux, equivalent width, etc.) for 
        this feature
        
        if fit is true, it specifies a model that should be used
        to fit the feature.  
        
        if interactive is True, a fitgui will be displayed to fit the
        data.  if fit is False, interactive will be ignored
        
        edges determine the behavior between the last pixel and the extent 
        edge:
        
        * 'tweak': adjust the extent to the nearest pixel and recompute
        * 'interp': interpolate the flux where the edges 
        * None/False/'': Do nothing
        """
        from .models import ParametricModel,LinearModel,get_model_instance
        from .utils import centroid
        
        if edges and edges not in ('interp','tweak'):
            raise ValueError('invalid value for edges parameter %s'%edges)
        
        if self.extent is None:
            raise ValueError('no extent available for line')
        x1,x2 = self.extent
        if x1 > x2:
            x1,x2 = x2,x1
            
        if self.continuum == 'fromspec' and spec.continuum is not None:
            cont = spec.continuum
        elif self.continuum == 'fromspec':
            cont = None
        else:
            cont = self.continuum
            
        xi1 = spec.x.searchsorted(x1)
        xi2 = spec.x.searchsorted(x2)
        
        if edges == 'tweak':
            self.extent = (spec.x[xi1],spec.x[xi2])
            return self.computeFeatureData(spec,fit,interactive,None)
        
        x = spec.x[xi1:xi2]
        y = spec.flux[xi1:xi2]
        err = spec.err[xi1:xi2]
        
        self.center = centroid(y,x)[0]
        self.centererr = None #TODO:figure out
        
        if isinstance(cont,tuple) and len(cont)==2:
            linmod = LinearModel()
            linmod.fitData((x1,x2),cont)
            cont = linmod
        
        if cont is None:
            cont = np.zeros(x.size)
        elif isinstance(cont,ParametricModel):
            cont = cont(x)
        else:
            cont = np.array(cont,copy=False)[xi1:xi2]
        conterr = 0 #TODO:figure out
        
        if fit:
            if fit is True:
                if self.model is None:
                    model = 'gaussian'
                else:
                    model = self.model
            else:
                model = get_model_instance(fit)
                
            if interactive:
                from .gui import fitgui
                fg = fitgui.FitGui(x,y-cont,weights=1/err,model=model)
                fg.configure_traits(kind='modal')
            else:
                model.fitData(x,y-cont,weights=1/err)
                
            lineflux = model.integrate(x1,x2)
            linefluxerr = 0 #TODO:figure out
            self.model = model
        else:
            # direct computation
            #TODO:fix for edges
            if xi1 == 0:
                dxi1 = 0
            else:
                dxi1 = xi1-1
            if xi2 == 0:
                dxi2 = 0
            else:
                dxi2 = xi2+1
            
            dx = np.convolve(spec.x[dxi1:dxi2],[0.5,0,-0.5],mode='valid') #TODO:figure out why sign is flipped
            lineflux = np.sum((y-cont)*dx)
            linefluxerr = 0 #TODO:figure out
            
            edges = False #TODO:remove when fixed
            if edges == 'interp':
                if xi1 == 0:
                    dx1 = spec.x[0] - x1
                    flux1 = spec.flux[0]*dxl
                else:
                    dx1 = x1 - spec.x[xi1]
                    flux1 = spec.flux[xi1]*dx1
                
                if xi2 == spec.x.size:
                    dx2 = x2 - spec.x[-1]
                    flux2 = spec.flux[-1]*dx2
                else:
                    dx2 = spec.x[xi2] - x2
                    flux2 = spec.flux[xi2]*dx2
                
                lineflux += flux1 + flux2
                
            self.model = None
            
        self.flux = lineflux
        self.fluxerr = linefluxerr
        
        mcont,mconterr = np.mean(cont),np.mean(conterr) #continuum should be array by here
        self.ew = lineflux/mcont 
        self.ewerr = linefluxerr/mcont - lineflux*mconterr*mcont**-2
    
    def identify(self,nameorknown):
        if isinstance(nameorknown,basestring):
            raise NotImplementedError('TODO:add stored line list')
            self.known = nameorknown
        elif isinstance(nameorknown,KnownFeature):
            self.known = nameorknown
        else:
            raise ValueError('unrecognized identification')
        
    def _getIdname(self):
        if self.known is None:
            return 'unknown'
        else:
            return self.known.name
    def _setIdname(self,val):
        if isinstance(val,basestring):
            self.identify(val)
        else:
            raise TypeError('feature identity name must be a string')
    idname = property(_getIdname,_setIdname,doc='name string of the identity of this feature')
    
    
class KnownFeature(HasSpecUnits):
    """
    Represents a spectral feature with known location to be matched
    to an unknown feature (or included in line lists).
    
    The strength attribute can be positive to imply emmission, negative to
    mean absorption, or None to leave unspecified
    """
    
    def __init__(self,loc,name='',strength=None,unit='wavelength'):
        super(KnownFeature,self).__init__(unit)
        
        self.loc = loc
        if name == '':
            name = str(loc)
        self.name = name
        
        self.strength = None
        
    def __str__(self):
        ls = str(int(round(self.loc)))
        if ls not in self.name:
            return '%s(%s)'%(self.name,ls)
        else:
            return self.name
            
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        if self.strength is None:
            self.loc = xtrans(self.loc)
        else:
            self.loc,self.strength = xftrans(self.loc,self.strength)
            
    def getUnitLoc(self,unit):
        """
        gets the location in the specified units (without changing the current
        unit)
        """
        if unit == self.unit:
            return self.loc
        else:
            oldunit = self.unit
            try:
                self.unit = unit
                return self.loc
            finally:
                self.unit = oldunit
            
def load_line_list(lst,unit='wavelength',tol=None,ondup='warn',sort=True):
    """
    Generates a spectral line list.
    
    lst can be:
    
    * either 'galaxy' or 'stellar' to load built-in line lists
      of common optical features
    * a filename for a file with one spectral feature per line - either 
      one column of the spectral location, two columns (location,name), or 
      three columns (location,name,strength)
    * a sequence of KnownFeature objects.
    
    tol is the seperation between locations allowed until two features
    are considered duplicates (or None to require an exact match)
    
    ondup specifies the action to take if a duplicate is encountered: 
    
    * 'raise':raise an IndexError
    * 'warn': issue a warning
    * 'remove':removes all but last in a set of duplicates
    * None: ignore
    
    The list will be sorted in increasing (decreasing) order if sort is 
    positive (negative), or no sorting will be done if it is 0/False
    
    returns a list of KnownFeature objects 
    """
    from warnings import warn
    from .io import _get_package_data
    
    if isinstance(lst,basestring):
        name = lst
        if name in ('galaxy','stellar'):
            #load builtin line list
            f = _get_package_data(name+'_lines.dat').split('\n')
            fopen = False
        else:
            f = open(name)
            fopen = True
            
        try:
            kfs = []
            for i,l in enumerate(f):
                llstrip = l.lstrip()
                if llstrip != '' and llstrip[0]!='#':
                    ls = l.split()
                    if len(ls)==1:
                        loc = float(ls[1])
                        name = ''
                        strength = None
                    elif len(ls)==2:
                        loc = float(ls[0])
                        name = ls[1]
                        strength = None
                    elif len(ls) == 3:
                        loc = float(ls[0])
                        name = ls[1]
                        strength = ls[2]
                    else:
                        raise ValueError('could not parse line #'+str(i))
                                           
                    kfs.append(KnownFeature(loc,name,strength,unit)) 
        finally:
            if fopen:
                f.close()
    else:
        kfs = lst
        
    if tol is None:
        tol = 0
    dupset = []
    ignores = []
    for i,kf in enumerate(kfs):
        if i not in ignores:
            locarr = np.array([kfi.loc for kfi in kfs])
            dups=[j for j,dup in  enumerate(np.abs(locarr-kf.loc) <= tol) if dup]
            if len(dups)>1:
                dups.append(i)
                dupset.append(set(dups))
                ignores.extend(dupset[-1])
    
    torem = []
    for ds in dupset:
        if ondup == 'raise':
            raise IndexError('Lines '+str(ds)[4:-1]+' are duplicates')
        elif ondup == 'warn':
            
            warn('Lines '+str(ds)[4:-1]+' are duplicates')
        elif ondup == 'remove':
            torem.extend(ds)
            del torem[-1]
        elif not ondup:
            pass
        else:
            raise ValueError('invalid ondup')
        
    for i in sorted(torem,reverse=True):
        del kfs[i]
        
    if sort:
        if int(sort) > 0: #increasing order
            kfs.sort(key=lambda o:o.loc)
        else: #decreasing
            kfs.sort(key=lambda o:-o.loc)
        
    return kfs
    
    
#<------------------------Spectrum-related functions--------------------------->

def align_spectra(specs,ressample='super',interpolation='linear',copy=False):
    """
    resample the spectra in the sequence specs so that they all have the same 
    x-axis.  
    
    ressample can be 'super' lower-resolution spectra to the highest or 'sub' 
    to interolate the highest-resolution spectrum to the lowest
    alternateively, 'logsuper' or 'logsub' will use the logarithmic resolution
    to determine
    
    copy makes new copies of the Spectrum objects 
    
    returns specs, or new copies if copy is True
    """
    from operator import isSequenceType
    if not isSequenceType(specs):
        raise ValueError('specs must be a sequence')
    
    for s in specs:
        if s.__class__.__name__ != 'Spectrum':
            raise ValueError(str(s)+' is not a spectrum')
    
    if copy:
        from copy import deepcopy
        specs= [deepcopy(s) for s in specs]
    
    if 'super' in ressample:
        super = True
    elif 'sub' in ressample:
        super = False
    else:
        raise ValueError('unrecognized ressample value')
    
    logres = 'log' in ressample
        
    
    if logres:
        reses=[s.getDlogx() for s in specs]
    else:
        reses=[s.getDx() for s in specs]
    mins=[np.min(s) for s in specs]
    maxs=[np.max(s) for s in specs]
    
    if super:
        templi = np.where(reses == max(reses))[0][0]
    else:
        templi = np.where(reses == max(reses))[0][0]
    
    x = specs[templi].x
    
    
    for s in specs:
        s.resample(x,interpolation)
    
    return specs
#<---------------------spectral utility functions------------------------------>

def air_to_vacuum(airwl,nouvconv=True):
    """
    Returns vacuum wavelength of the provided air wavelength array or scalar.
    Good to ~ .0005 angstroms.

    If nouvconv is True, does nothing for air wavelength < 2000 angstroms.
    
    Input must be in angstroms.
    
    Adapted from idlutils airtovac.pro, based on the IAU standard 
    for conversion in Morton (1991 Ap.J. Suppl. 77, 119)
    """
    airwl = np.array(airwl,copy=False,dtype=float,ndmin=1)
    isscal = airwl.shape == tuple()
    if isscal:
        airwl = airwl.ravel()
    
    #wavenumber squared
    sig2 = (1e4/airwl)**2
    
    convfact = 1. + 6.4328e-5 + 2.94981e-2/(146. - sig2) +  2.5540e-4/( 41. - sig2)
    newwl = airwl.copy() 
    if nouvconv:
        convmask = newwl>=2000
        newwl[convmask] *= convfact[convmask]
    else:
        newwl[:] *= convfact
    return newwl[0] if isscal else newwl

def vacuum_to_air(vacwl,nouvconv=True):
    """
    Returns air wavelength of the provided vacuum wavelength array or scalar.
    Good to ~ .0005 angstroms.
    
    If nouvconv is True, does nothing for air wavelength < 2000 angstroms.
    
    Input must be in angstroms.
    
    Adapted from idlutils vactoair.pro.
    """
    vacwl = np.array(vacwl,copy=False,dtype=float)
    isscal = vacwl.shape == tuple()
    if isscal:
        vacwl = vacwl.ravel()
    
    #wavenumber squared
    wave2 = vacwl*vacwl
    
    convfact = 1.0 + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2*wave2)
    newwl = vacwl.copy()/convfact
    
    if nouvconv:    
        #revert back those for which air < 2000
        noconvmask = newwl<2000
        newwl[noconvmask] *= convfact[noconvmask] 
    
    return newwl[0] if isscal else newwl

def zfind(specobj,templates,lags=(0,200),checkspec=True,checktemplates=True,verbose=True,interpolation = None):
    """
    computes the best fit by linear least-squares fitting of templates to the 
    spectrum for each possible pixel offset.  Weighted fits will be done if 
    the spectrum ivars are different.
    
    lags can either be a sequence of lags or a 2-tuple specifying the lower and 
    upper possible lags
    
    specobj must be a Spectrum object or a sequence of (flux,[x],[ivar]) or 
    flux.  If it is not logarithmically spaced, it will be interpolated
    
    templates can be either a sequence of Spectrum objects or an array with at
    least one dimension matching the pixel dimension.  
    (Note that templates longer than the spectrum will not use information off
    the edges)
    
    interpolation is the technique for interpolating for sub-pixel lags.  If 
    None, no interpolation is used, so lags must be integers (but this method is
    much faster)
    
    returns besti,lags,zs,coeffs,xs,fitfluxes,rchi2s
    """
    if interpolation is not None:
        raise NotImplementedError('interpolation not implemented yet')
    
    from operator import isSequenceType
    if not isinstance(specobj,Spectrum) and isSequenceType(specobj):
        if len(specobj) > 3:
            specobj = Spectrum(np.logspace(0,1,len(specobj)),specobj)
        else:
            flux = specobj[0]
            if len(specobj) < 3:
                ivar = None
            else:
                ivar = specobj[2]
            if len(specobj) < 2:
                x = np.arange(len(flux))
            else:
                x = specobj[1]
            specobj = Spectrum(x,flux,ivar=ivar)
        
    x = specobj.x
    if checkspec:
        if not specobj.isLogarithmic():
            newx = np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),len(x))
            specobj = specobj.copy()
            specobj.resample(newx)
            x = specobj.x
    flux = specobj.flux
    npix = specobj.shape[0]
    ivar = specobj.ivar
    
    if checktemplates:
        templates=[t if isinstance(t,Spectrum) else Spectrum(x,t) for t in templates]
        for i,t in enumerate(templates):
            if not t.isXMatched(x):
                if verbose:
                    print 'template',i,'does not match spectrum -- resampling'
                t.resample(x)
        templates = [t.flux for t in templates]
    templates=np.array(templates,ndmin=1)
    
    if templates.shape[0] == npix:
        tm = np.asmatrix(templates)
    elif templates.shape[1] == npix:
        tm = np.asmatrix(templates).T
        
    y = np.asmatrix(flux).T
    ws = np.asmatrix(ivar).T
    
    
    if type(lags) is tuple and len(lags) == 2:
        lags = np.arange(*lags)
    llags = lags[lags<0]
    ulags = lags[lags>0]
    
    #generate slice objects to match offsets to lags
    #TODO:index directly - will this be much slower than rolling?
    ls,slices=[],[]
    for l in llags:
        ls.append(l)
        slices.append((np.s_[-l:,:],np.s_[:l,:]))
    if 0 in lags:
        ls.append(0)
        slices.append((np.s_[:,:],np.s_[:,:]))
    for l in ulags:
        ls.append(l)
        slices.append((np.s_[:-l,:],np.s_[l:,:]))
        
        
    cs,dsq,fitfluxes=[],[],[]
    
    #don't do weighting if all of the errors are identical -- matrix becomes singular
    useweights = np.any(ivar-ivar[0]) and not np.all(~np.isfinite(ivar)) 
    for l,s in zip(ls,slices):
        if verbose:
            print 'doing lag',l
        A = tm[s[0]]
        v = y[s[1]]
        w = ivar[s[1][0]]
        
        if useweights:
            try:
                AT = np.multiply(A.T,w)
                cs.append(np.linalg.inv(AT*A)*AT*v)
            except np.linalg.LinAlgError,e:
                if verbose:
                    print 'Error inverting matrix in lag',l,':',e
                cs.append(np.linalg.pinv(A)*v)
        else:
            cs.append(np.linalg.pinv(A)*v)
            #v->4096,a->3096
            #TODO: faster inversion schemes?
        
        fitspec = A*cs[-1]
        fitfluxes.append(fitspec)
        fitdiff = (v-fitspec).A
        dsq.append(np.sum(fitdiff*fitdiff))
    ls=np.array(ls)
    cs=np.array(cs)
    dsq=np.array(dsq)
    dofs = np.sum(ivar!=0) - np.abs(ls)
    rchi2s = dsq/dofs
    
    mins = np.where(rchi2s==min(rchi2s))[0]
    if len(mins) == 1:
        besti = mins[0]
    else:
        besti = mins
        
    xs = [x[s[1][0]] for s in slices]
    fitfluxes = [f.A[:,0] for f in fitfluxes]
    zs = np.mean(lag_to_z(x,ls),1)
    
    try:
        from collections import namedtuple
        tinit = namedtuple('zfind_out','besti lags zs coeffs xs fitfluxes rchi2s')
    except ImportError: #support for pre-2.6 - use ordinary tuples
        tinit = lambda *args:args
    return tinit(besti,ls,zs,cs,xs,fitfluxes,rchi2s)

def lag_to_z(x,lag,xunit='ang',avgbad=True):
    """
    this converts an integer pixel lag for a given x-axis into a 
    redshift
    
    currently the only supported unit for the x-axis is angstroms
    
    avgbad replaces undeterminable velocities with the average of 
    the others, else they are set to 0
    """
    x = np.array(x,copy=False)
    
    if np.isscalar(lag):
        scalarout = True
        lag = [lag]
    else:
        scalarout = False
        
    zs = np.ndarray((len(lag),x.size))
    for i,l in enumerate(lag):
        z = np.roll(x,-l)/x-1
        
        if l>=0:
            if avgbad:
                z[-l:] = np.mean(z[:-l])
            else:
                z[-l:] = 0
        else:
            if avgbad:
                z[:-l] = np.mean(z[-l:])
            else:
                z[:-l] = 0
            
        zs[i] = z
    
    return zs[0] if scalarout else zs


def build_model_spectrum(n=1024,T=5800,range=None,peak=1,z=0,noise=None,lines=None,linesampling=1,name='Model'):
    """
    Constructs a basic spectral object composed of a blackbody with temperature
    T, a set of lines, and noise.
    
    `peak` sets the peak of the blackbody curve, z is the redshift.
    
    `noise` can be:
    
    * positive
        additive gaussian noise with the specified amplitude
    * negative
        multiplicative uniform noise with the specified amplitude
    * True
        poisson
    * False/None
        no noise
    
    `lines` can be None or a sequence, where each element of the sequence is:
    
    * a 3-sequence (location,amplitude,width) where if width is positive it
      gives a gaussian sigma, or negative gives a lorentzian FWHM
    * a 4-sequence (location,amplitude,sigma,FWHM) to use a voigt profile
    * a  FunctionModel1D for the line flux
    
    linesampling is passed into FunctionModel1D.pixelize or the models will
    be called directly if it is 1
    """
    from . import models
    bm = models.BlackbodyModel(T=T)
    bm.peak = peak
    
    if range is None:
        range = bm.rangehint
    x = np.linspace(range[0],range[1],n)
    x0 = x/(z+1)
    
    flux = bm(x0)
    
    if lines:
        for l in lines:
            if isinstance(l,models.FunctionModel1D):
                lmod = l
            elif len(l) == 3:
                loc,A,w = l
                if w>0:
                    lmod = models.GaussianModel(mu=loc,sig=w)
                    lmod.peak = A
                else:
                    lmod = models.LorentzianModel(mu=loc,gamma=-w)
                    lmod.peak = A
            elif len(l) ==4 :
                loc,A,sig,w = l
                lmod = models.VoigtModel(mu=loc,sig=sig,gamma=w)
                lmod.peak = A
            else:
                raise ValueError('unrecognized line specifier %s'%l)
            if linesampling == 1:
                flux += lmod(x0)
            else:
                flux += lmod.pixelize(x0,sampling=linesampling)
    if noise:
        if noise is True:
            flux = np.random.poisson(flux)
            err = flux**0.5
        else:
            noise = float(noise)
            if noise > 0:
                err = np.random.randn(flux.size)*noise
            else:
                err = flux*(1+(2*np.random.rand(flux.size)-1)*-noise) - flux
            flux = flux + err
    else:
        err = None
                
    return Spectrum(x,flux,err,name=name)



del ABCMeta,abstractmethod,abstractproperty #clean up namespace
