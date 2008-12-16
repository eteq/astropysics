from __future__ import division
from math import pi
import numpy as np

try:
    import spylot
except ImportError:
    print 'Spylot not found - some spec package functions may not work correctly'
    
class Flux(object):
    def __init__(self,x,flux,type='wl'):
        
        x,flux = np.array(x),np.array(flux)
        if x.shape != flux.shape:
            raise ValueError("x and flux don't match shapes")
        
        
        self._type,self._unit,self._scaling = self._strToType(type)
        
        self._x = x
        self._flux = flux
        
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
        elif inunit == 'm' or u == 'wavelength-m':
            val =  'wavelength-m'
            scaling=1e2
        elif inunit == 'cm' or u == 'wavelength-cm':
            val =  'wavelength-cm'
            scaling=1
        elif inunit == 'f' or inunit == 'nu' or u == 'hz' or u == 'frequency-hz':
            val =  'frequency-Hz'
            scaling=1
        elif u == 'thz' or inunit == 'frequency-THz':
            val =  'frequency-THz'
            scaling=1e12
        elif u == 'e' or inunit == 'energy-eV':
            from .constants import ergperev
            val =  'energy-eV'
            scaling=ergperev
        elif u == 'erg' or inunit == 'energy-erg':
            val =  'energy-erg'  
            scaling=1
        elif inunit == 'J' or inunit == 'energy-J':
            val =  'energy-J'   
            scaling=1e-7
        else:
            raise ValueError('unrecognized type')
        
        tyep,unit = val.split('-')
        return type,unit,scaling
        
    def _convertType(self,oldtype,newtype,oldscale,newscale):
        from .constants import c,h
        
        x,flux = self._x/oldscale,self._flux*oldscale # convert to cgs
    
        if newtype == 'energy': #TODO:check
            if oldtype == 'energy': 
                pass
            elif oldtype == 'wavelength':
                x = h*c/x
                flux = x*x*flux/c/h #convert to fnu and divide by h
                
            elif oldtype == 'frequency':
                x = h*x
                flux = flux/h
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'wavelength':
            if oldtype == 'energy': #TODO:check
                x = h*c/x
                flux = x*x*flux*h/c #convert to fnu and then to  by h
            elif oldtype == 'wavelength':
                pass
            elif oldtype == 'frequency':
                x = c/x
                flux = x*x*flux/c #flambda = fnu*nu^2/c
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'frequency':
            if oldtype == 'energy': #TODO:check
                x = x/h
                flux = flux*h
            elif oldtype == 'wavelength':
                x = c/x
                flux = x*x*flux/c #fnu = flambda*lambda^2/c
            elif oldtype == 'frequency':
                pass
            else:
                raise ValueError('unrecognized oldtype')
        else:
            raise ValueError('unrecognized newtype')
        
        return x*newscale,flux/newscale
        
    def _getFlux(self):
        return self._flux.copy()
    def _setFlux(self,flux):
        flux = np.array(flux)
        if flux.shape != self._x.shape:
            raise ValueError("new flux doesn't match x shape")
        self._flux = flux
    flux = property(_getFlux,_setFlux)
    
    def _getNFlux(self):
        raise NotImplementedError
    def _setNFlux(self,nflux):
        nflux = np.array(nflux)
        if nflux.shape != self._x.shape:
            raise ValueError("new flux doesn't match x shape")
        raise NotImplementedError
        self._flux = something
    nflux = property(_getNFlux,_setNFlux)
    
    def _getX(self):
        return self._x.copy()
    def _setX(self,x):
        x = np.array(x)
        if x.shape != self._flux.shape:
            raise ValueError("new x doesn't match flux shape")
        self._x = x
    x = property(_getX,_setX)
    
    def _getUnit(self):
        return self._type+'-'+self._unit
    def _setUnit(self,typestr):
        newtype,newunit,newscaling = self._strToType(typestr)
        
        x,flux = self._convertType(self._type,newtype,self._scaling,newscaling)
        
        self._x = x
        self._flux = flux
        self._type = newtype
        self._unit = newunit
        self._scaling = newscaling
    unit = property(_getUnit,_setUnit)
    
    def smooth(self,width=1,filtertype='gaussian',replace=False):
        """
        smooths the flux in this object by a filter of the given filtertype 
        (can be either 'gaussian' or 'boxcar')
        
        if replace is True, the flux in this object is replaced by the smoothed flux
        
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
        
        return smoothed
    
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
        
        fobj = Flux(x,flux)
        fobj.ivar = ivar
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