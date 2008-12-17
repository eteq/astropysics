from __future__ import division
from math import pi
import numpy as np

try:
    import spylot
except ImportError:
    print 'Spylot not found - some spec package functions may not work correctly'
    
class Flux(object):
    def __init__(self,x,flux,ivar=None,type='wl'):
        
        x,flux = np.array(x),np.array(flux)
        if x.shape != flux.shape:
            raise ValueError("x and flux don't match shapes")
        
        if ivar is None:
            ivar = np.zeros_like(flux)
        elif ivar.shape != flux.shape:
            raise ValueError("ivar and flux don't match shapes")
        
        self._type,self._unit,self._scaling = self._strToType(type)
        
        self._x = x
        self._flux = flux
        self._err = ivar**-0.5
        
    def _strToType(self,typestr):
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
        elif u == 'e' or u == 'energy-eV':
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
            raise ValueError('unrecognized type')
        
        tyep,unit = val.split('-')
        return type,unit,scaling
        
    def _convertType(self,oldtype,newtype,oldscale,newscale):
        from .constants import c,h
        
        x,flux,err = self._x/oldscale,self._flux*oldscale,self._err*oldscale # convert to cgs
    
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
        
        return x*newscale,flux*fluxscale/newscale,err*fluxscale/newscale
        
    def _getFlux(self):
        return self._flux.copy()
    def _setFlux(self,flux):
        flux = np.array(flux)
        if flux.shape != self._x.shape:
            raise ValueError("new flux doesn't match x shape")
        self._flux = np.array(flux)
    flux = property(_getFlux,_setFlux)
    
    def _getNFlux(self):
        raise NotImplementedError
    def _setNFlux(self,nflux):
        nflux = np.array(nflux)
        if nflux.shape != self._x.shape:
            raise ValueError("new flux doesn't match x shape")
        raise NotImplementedError
        self._flux = 42
    nflux = property(_getNFlux,_setNFlux)
    
    def _getErr(self):
        return self._err.copy()
    def _setErr(self,err):
        if err.shape != self._flux.shape:
            raise ValueError("new err doesn't match flux shape")
        self._err=np.array(err)
    err = property(_getErr,_setErr)
    
    def _getIvar(self):
        return 1/self._err/self._err
    def _setIvar(self,ivar):
        if ivar.shape != self._flux.shape:
            raise ValueError("newivar doesn't match flux shape")
        self._err=np.array(ivar)**-0.5
    ivar = property(_getIvar,_setIvar)
    
    def _getX(self):
        return self._x.copy()
    def _setX(self,x):
        x = np.array(x)
        if x.shape != self._flux.shape:
            raise ValueError("new x doesn't match flux shape")
        self._x = np.array(x)
    x = property(_getX,_setX)
    
    def _getUnit(self):
        return self._type+'-'+self._unit
    def _setUnit(self,typestr):
        newtype,newunit,newscaling = self._strToType(typestr)
        
        x,flux,err = self._convertType(self._type,newtype,self._scaling,newscaling)
        
        self._x = x
        self._flux = flux
        self._err = err
        self._type = newtype
        self._unit = newunit
        self._scaling = newscaling
    unit = property(_getUnit,_setUnit)
    
    #<----------------------Tests/Info--------------------------------->
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
        
    def getDlogx(self,mean=True):
        """
        get the logarithmic spacing of the x-axis, which is always 1 element
        shorter than x
        
        if mean, returns the mean of the spacing
        """
        x = self._x
        dlogx = np.convolve(x,(1,-1),mode='valid')/x[1:]
        if mean:
            return dlogx.mean()
        else:
            return dlogx
    
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
    
    def smooth(self,width=1,filtertype='gaussian',replace=False):
        """
        smooths the flux in this object by a filter of the given filtertype 
        (can be either 'gaussian' or 'boxcar')
        
        if replace is True, the flux in this object is replaced by the smoothed
        flux and the error is smoothed in the same fashion
        
        returns smoothedflux
        """
        import scipy.ndimage as ndi 
        
        if filtertype == 'gaussian':
            filter = ndi.gaussian_filter1d
        elif filtertype == 'boxcar' or type == 'uniform':
            filter = ndi.uniform_filter1d
        else:
            raise ValueError('unrecognized filter type %s'%filtertype)
        
        smoothed = filter(self._flux,width)
        
        if replace:
            self.flux = smoothed
            self.err = filter(self._err,width)
        
        return smoothed
    
    def linearize(self,lower=None,upper=None,n=None):
        """
        convinience function for resampling to an equally-spaced linear x-axis
        
        if lower or upper are None, the upper and lower x values are used
        if n is None, the same number of samples are used
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        if n is None:
            n = self._x.shape[0]
        
        newx = np.linspace(lower,upper,n)
        return self.resample(newx)
    
    def logify(self,lower=None,upper=None,n=None):
        """
        convinience function for resampling to an x-axis that is equally spaced
        in logarithmic bins.  Note that lower and upper are the x-axis values 
        themselves, NOT log(xvalue)
        
        if lower or upper are None, the upper and lower x values are used
        if n is None, the same number of samples are used
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        if n is None:
            n = self._x.shape[0]
        
        newx = np.logspace(np.log10(lower),np.log10(upper),n)
        return self.resample(newx)
    
    def resample(self,newx,interpolation='linear'):
        """
        this interpolates the flux to populate the supplies x-axis
        
        currently the only available interpolation is 'linear'
        
        WARNING: this does not treat the errors properly yet - currently just interpolates
        """
        newflux = np.interp(newx,self._x,self._flux)
        #TODO: fix errors
        newerr = np.interp(newx,self._x,self._err)
        
        self._x = newx 
        self._flux = newflux
        self._err = newerr       
        
    
    def plot(self,fmt='-',clf=True,log=False,**kwargs):
        import matplotlib.pyplot as plt
        
        if clf:
            plt.clf()
        
        if log:
            plot = plt.semilogx
        else:
            plot = plt.plot
        return plot(self.x,self.flux,fmt,**kwargs)

def load_deimos_spectrum(fn,plot=True,extraction='horne',retdata=False,smoothing=None):
    """
    extraction type can 'horne' or 'boxcar'
    
    returns Flux object with ivar, [bdata,rdata]
    """
    import pyfits
    if 'spec1d' not in fn or 'fits' not in fn:
        raise ValueError('loaded file must be a 1d spectrum from deep pipeline')
    
    if extraction == 'horne':
        extname = 'HORNE'
    elif extraction == 'boxcar':
        extname = 'BXSPF'
    else:
        raise ValueError('unrecgnized extraction type %s'%extraction)
    
    f=pyfits.open(fn)
    try:
        extd = dict([(f[i].header['EXTNAME'],i) for i in range(1,len(f))])
        bi,ri = extd[extname+'-B'],extd[extname+'-R']
        bd,rd=f[bi].data,f[ri].data
        x=np.concatenate((bd.LAMBDA[0],rd.LAMBDA[0]))
        flux=np.concatenate((bd.SPEC[0],rd.SPEC[0]))
        ivar=np.concatenate((bd.IVAR[0],rd.IVAR[0]))
        sky=np.concatenate((bd.SKYSPEC[0],rd.SKYSPEC[0]))
        
        changei = len(bd.LAMBDA[0])
        
        fobj = Flux(x,flux,ivar)
        fobj.sky = sky
        
        if smoothing:
            fobj.smooth(smoothing,replace=True)
        
        if plot:
#            from operator import isMappingType
#            if not isMappingType(plot):
#                plot={}
            from matplotlib import pyplot as plt
            plt.plot(fobj.x[:changei],fobj.flux[:changei],'-b')
            plt.plot(fobj.x[changei:],fobj.flux[changei:],'-r')
        
        if retdata:
            return fobj,bd,rd
        else:
            return fobj
    finally:
        f.close()
        
def load_deimos_templates(fns):
    """
    This will generate a dictionary of Flux objects from the specified templates 
    fn can be a single file, a list of files, or a unix-style pattern
    """
    import pyfits
    from operator import isSequenceType
    from warnings import warn
    
    if type(fns) == str:
        from glob import glob
        fns = glob(fns)
        
    if not isSequenceType(fns):
        raise ValueError('improper filename format')
    
    tempd={}
    xd={}
    for fn in fns:
        f = pyfits.open(fn)
        try:
            h=f[0].header
            d=f[0].data
            if len(d.shape) == 1:
                d = (d,)
            if len(d) == 1:
                if h.has_key('NAME0'):
                    k=h['NAME0'].strip()
                elif h.has_key('NAME'):
                    k=h['NAME'].strip()
                elif h.has_key('OBJECT'):
                    k=h['OBJECT'].strip()
                else:
                    k = fn.strip()
                    warn('could not infer spectrum name - using filename %s as key name'%k)
                    
                xd[k] = 10**(h['COEFF0']+h['COEFF1']*np.arange(d[0].shape[-1]))
                tempd[k] = d[0].copy()
            else:
                x = 10**(h['COEFF0']+h['COEFF1']*np.arange(d.shape[-1]))
                for i,spec in enumerate(d):
                    if  h.has_key('NAME%i'%i):
                        k=h['NAME%i'%i].strip()
                    elif h.has_key('OBJECT') and h.has_key('EIGEN%i'%i):
                        k = '%s%i'%(h['OBJECT'].strip(),h['EIGEN%i'%i])
                    elif h.has_key('NAME'):
                        k = '%s-%i'%(h['NAME'].strip(),i)
                    elif h.has_key('OBJECT'):
                        k = '%s-%i'%(h['OBJECT'].strip(),i)
                    else:
                        k = '%s-%i'%(fn.strip(),i)
                        warn('could not infer spectrum name - using filename %s as key name'%k)
                        
                    xd[k] = x
                    tempd[k]=spec.copy()
                    
        finally:
            f.close()
    
    return dict(((k,Flux(xd[k],tempd[k])) for k in tempd))
    
    
         