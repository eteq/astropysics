#Copyright (c)2008 Erik Tollerud (etolleru@uci.edu) 
"""
This module is for observed or synthetic spectra (Spectrum object) and 
related operations.

Methods/functions tend to be oriented towards optical, but will work in other
bands.
"""

from __future__ import division,with_statement
from math import pi
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
from .io import load_deimos_spectrum,load_all_deimos_spectra,load_spylot_spectrum

class HasSpecUnits(object):
    """
    Use this class as a mixin superclass for objects that have spectrum-like 
    units (e.g. flux spectral density f_lambda/f_nu, luminosity spectral 
    density, etc.).  It adds a property 'unit' that can be set 
    
    The following must be done by subclasses:
    *override method _applyUnits(self,xtrans,xitrans,xftrans,xfinplace) - see method doc for details
    *call HasSpecUnits.__init__(self,unit) in the class 
    initializer where unit is the units of the initial input data
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
        and oldx = xtrans(newx)
        
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
        converts a unit string into a standardized for
        
        returns phystype,unit,scaling
        (scaling relative to cgs)
        """
        if not isinstance(typestr,basestring):
            raise TypeError('invalid type provided as unit')
        
        u = typestr.lower()
        if u == 'wl' or u == 'lambda' or u == 'ang' or u == 'angstroms' or u == 'wavelength-angstrom':
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
        elif u == 'f' or u == 'nu' or u == 'hz' or u == 'frequency-hz':
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
    def __init__(self,x,flux,err=None,ivar=None,unit='wl',copy=True):
        """
        sets the x-axis values and the flux.  Optionally, an error can be
        given by supplying an error or an inverse variance (bot not both)
        
        copy determines if the inputs will be copied if they are already arrays
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
        
        self._x = x
        self._flux = flux
        self._err = err
        
        self.continuum = None
        
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
    x = property(_getX,_setX)
    
    
    
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
        
        width is in pixels
        
        returns smoothedflux,smoothederr
        """
        import scipy.ndimage as ndi 
        
        if filtertype == 'gaussian':
            filter = ndi.gaussian_filter1d
        elif filtertype == 'boxcar' or type == 'uniform':
            filter = ndi.uniform_filter1d
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
    
    def fitContinuum(self,model='uniformknotspline',weighted=False,evaluate=False,
                          interactive=False,**kwargs):
        """
        this method computes a continuum fit to the spectrum using a model
        from astropysics.models (list_models will give all options) or
        an callable with a fitData(x,y) function
        
        kwargs are passed into the constructor for the model
        
        if weighted, the inverse variance will be used as weights to the 
        continuum fit
        
        if interactive is True, the fitgui interface will be displayed to 
        tune the continuum fit
        
        the fitted model is assigned to self.continuum or evaluated at the
        spectrum points if evaluate is True and the results are set to 
        self.continuum
        """
        #for the default, choose a reasonable number of knots
        if model == 'uniformknotspline' and 'nknots' not in kwargs:
            kwargs['nknots'] = 4
        
        if isinstance(model,basestring):
            from .models import get_model
            model = get_model(model)(**kwargs)
        
        if not (callable(model) and hasattr(model,'fitData')):
            raise ValueError('provided model object cannot fit data')
        
        model.fitData(self.x,self.flux,weights=(self.ivar if weighted else None))
        
        if interactive:
            from .gui import fit_data,FitGui
            #model = fit_data(self.x,self.flux,model,weights=(self.ivar if weighted else None))
            fg = FitGui(self.x,self.flux,model,weights=(self.ivar if weighted else None))
            fg.plot.plots['data'][0].marker = 'dot'
            fg.plot.plots['data'][0].marker_size = 2
            fg.plot.plots['model'][0].line_style = 'solid'
            fg.configure_traits(kind='livemodal')
            model = fg.tmodel.model
            
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
        
    def rejectOutliersFromContinuum(self,sig=3):
        raise NotImplementedError
            
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
        
    def plot(self,fmt=None,ploterrs=.1,plotcontinuum=True,smoothing=None,clf=True,colors=('b','g','r','k'),**kwargs):
        """
        uses matplotlib to plot the Spectrum object
        
        if the linestyle is 'steps' (the default), the spectrum will
        be offset so that the steps are centered on the pixels
        
        colors should be a 3-tuple that applies to 
        (spectrum,error,invaliderror,continuum) and kwargs go into spectrum
        and error plots
        
        if ploterrs or plotcontinuum is a number, the plot will be 
        scaled so that the mean value matches the mean  of the spectrum 
        times the numeric value
        if True, the scaling will match the actual value
        if False, the plots will not be shown
        """
        
        import matplotlib.pyplot as plt
        
        if smoothing:
            x,(y,e) = self.x,self.smooth(smoothing,replace=False)
        else:
            x,y,e = self.x,self.flux,self.err
        
        if 'ls' not in kwargs and 'linestyle' not in kwargs:
            kwargs['ls']='steps'
            
        if kwargs['ls']=='steps' and len(x)>=4:
            x = np.hstack((1.5*x[0]-x[1]/2,np.convolve(x,[0.5,0.5],mode='valid'),1.5*x[-1]-x[-2]/2))
            y = np.hstack((y[0],y))
            e = np.hstack((e[0],e))
        elif len(x)==3:
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
            res = [plt.plot(x,y,**kwargs)]
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
                res.append(plt.plot(x[m],scale*e[m],**kwargs))
            if np.sum(~m) > 0:
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


class FunctionSpectrum(Spectrum):
    """
    This is a Spectrum generated by functional forms rather than
    using fixed samples.  Setting x will determine where the function is
    evaluated (and do the evaluation)
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
    A Spectrum that follows the provided astropysics.models.FunctionModel1D .
    
    noisefunc should be of the form nf(x,flux)
    
    Attributes
    ------------
    `noisetype` determines the noise model - can be:
    *None: no noise
    *'poisson': poisson noise where the noisefunc determines the poissonian
    scale factor
    *'normal' or 'gauss': normally-distributed noise where the noisefunc 
    determines sigma
    """
    def __init__(self,xi,model,noisetype=None,noisefunc=None,unit='wl'):     
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
    This class represents a Spectral Feature.
    
    All SpecralFeatures have a rest value,observed value (possibly with
    error), flux value (possibly with error), and equivalent width (possibly
    with error)
    
    Note that equivalent width is always expected to be in angstroms
    """
    def __init__(self,unit='wavelength'):
        HasSpecUnits.__init__(self,unit)
        
        self.rest = 1
        self.observed = None
        self.observederr = None
        self.flux = 0
        self.fluxerr = 0
        self.ew = None
        self.ewerr = None
    
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        self.rest = xtrans(self.rest)
            
        if observed is not None:
            oldobs = self.observed
            self.observed = newobs = xtrans(self.observed)
            if self.observederr is not None:
                op = oldobs+self.observederr
                om = oldobs-self.observederr
                op,om = xtrans(op),xtrans(om)
                self.observederr = (op+om)/2
    
    @property
    def name(self):
        raise NotImplementedError
    
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
    
    returns besti,lags,coeffs,xs,fitfluxes,rchi2s
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
    templates=np.atleast_1d(templates)
    
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
    


del ABCMeta,abstractmethod,abstractproperty #clean up namespace