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
        
        self._x = x/self._scaling
        self._flux = flux/self._scaling
        
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
        
    def _convertType(self,oldtype,newtype):
        from .constants import c,h
    
        if newtype == 'energy':
            if oldtype == 'energy':
                pass
            elif oldtype == 'wavelength':
                raise NotImplementedError
            elif oldtype == 'frequency':
                raise NotImplementedError
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'wavelength':
            if oldtype == 'energy':
                raise NotImplementedError
            elif oldtype == 'wavelength':
                pass
            elif oldtype == 'frequency':
                self._x = c/self._x
                self._flux = self._x*self._x*self._flux/c #flambda = fnu*nu^2/c
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'frequency':
            if oldtype == 'energy':
                raise NotImplementedError
            elif oldtype == 'wavelength':
                self._x = c/self._x
                self._flux = self._x*self._x*self._flux/c #fnu = flambda*lambda^2/c
            elif oldtype == 'frequency':
                pass
            else:
                raise ValueError('unrecognized oldtype')
        else:
            raise ValueError('unrecognized newtype')
        
    def _getFlux(self):
        return self._flux.copy()*self._scaling
    def _setFlux(self,flux):
        flux = np.array(flux)
        if flux.shape != self._x.shape:
            raise ValueError("new flux doesn't match x shape")
        self._flux = flux/self._scaling
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
        return self._x.copy()*self._scaling
    def _setX(self,x):
        x = np.array(x)/self._scaling
        if x.shape != self._flux.shape:
            raise ValueError("new x doesn't match flux shape")
        self._x = x
    x = property(_getX,_setX)
    
    def _getUnit(self):
        return self._type+'-'+self._unit
    def _setUnit(self,typestr):
        newtype,newunit,newscaling = self._strToType(typestr)
        
        self._convertType(self._type,newtype)
        
        self._type = newtype
        self._unit = newunit
        self._scaling = newscaling
    unit = property(_getUnit,_setUnit)