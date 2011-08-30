#Copyright 2011 Erik Tollerud
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

#These are all the builtins that should by in astropysics instead of pymodelfit

"""

==================================================
models -- astronomy-specific models for pymodelfit
==================================================

This module makes use of the :mod:`pymodelfit` package to define a number of
models inteded for data-fitting that are astronomy-specific. Note that
:mod:`pymodelfit` was originally written as part of astropysics, and hence it is
tightly integrated (including the Fit GUI, where relevant) to other parts of
astropysics, often using the models in this module.

.. note::

    Everything in the :mod:`pymodelfit` package is imported to this package.
    This is rather bad form, but is currently done because other parts of
    astropysics is written for when pymodelfit was part of astropysics.  This 
    may change in the future, though, so it is recommended that user code always
    use ``from pymodelfit import ...`` instead of 
    ``from astropysics.models import ...``.

Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.models
   :parts: 1
   
Module API
^^^^^^^^^^

"""
from __future__ import division,with_statement
import numpy as np
from pymodelfit.core import *
from pymodelfit.builtins import *
from .spec import HasSpecUnits as _HasSpecUnits

    
class BlackbodyModel(FunctionModel1DAuto,_HasSpecUnits):
    """
    A Planck blackbody radiation model.  

    Output/y-axis is taken to to be specific intensity.
    """
    from .constants import h,c,kb
    
    def __init__(self,unit='wl'):
        _HasSpecUnits.__init__(self,unit)
        self.unit = unit #need to explicitly set the unit to initialize the correct f
        self.stephanBoltzmannLaw = self._instanceSBLaw
        
    def f(self,x,A=1,T=5800):
        x = x*self._xscaling
        if self._phystype == 'wavelength': 
            val = self._flambda(x,A,T)
        elif self._phystype == 'frequency':
            val = self._fnu(x,A,T)
        elif self._phystype == 'energy':
            val = self._fen(x,A,T)
        else:
            raise ValueError('unrecognized physical unit type!')
        return val*self._xscaling
    
    def _flambda(self,x,A=1,T=5800):
        h,c,k=self.h,self.c,self.kb
        return A*2*h*c*c*x**-5/(np.exp(h*c/(k*T*x))-1)
    
    def _fnu(self,x,A=1,T=5800):
        h,c,k=self.h,self.c,self.kb
        return A*2*h/c/c*x**3/(np.exp(h*x/(k*T))-1)
    
    def _fen(self,x,A=1,T=5800):
        return self._fnu(x,A,T)/self.h
    
    def _applyUnits(self,xtrans,xitrans,xftrans,xfinplace):
        pass #do nothing because the checking is done in the f-function
#        if self._phystype == 'wavelength': #TODO:check
#            self.f = self._flambda
#        elif self._phystype == 'frequency':
#            self.f = self._fnu
#        elif self._phystype == 'energy':
#            self.f = self._fen
#        else:
#            raise ValueError('unrecognized physical unit type!')
    @property
    def xaxisname(self):
        if self._phystype == 'wavelength': 
            return 'lambda'
        elif self._phystype == 'frequency':
            return 'nu'
        elif self._phystype == 'energy':
            return 'E'
        
    yaxisname = 'I'
    
    @property
    def rangehint(self):
        cen = self.wienDisplacementLaw(None)
        return(cen/3,cen*3)
    
    
    def setIntensity(self):
        """
        Sets A so that the output is specific intensity/surface brightness.
        """
        self.A = 1
    
    def setFlux(self,radius,distance):
        """
        Sets A so that the output is the flux at the specified distance from
        a spherical blackbody with the specified radius.
        
        :param radius: Radius of the blackbody in cm
        :type radius: float
        :param distance: distance to the blackbody in cm
        :type distance: float
        """
        from .phot import intensity_to_flux
        self.A = intensity_to_flux(radius,distance)
        
    def getFlux(self,x,radius=None,distance=None):
        """
        Returns the flux density at the requested wavelength for a blackbody 
        of the given radius at a specified distance.
        
        :param x: x-value of the model
        :type x: float
        :param radius: radius of the blackbody in cm
        :type radius: float
        :param distance: distance to the blackbody in cm
        :type distance: float
        
        :returns: flux/wl at the specified distance from the blackbody
        
        """
        if distance is None:
            if radius is None:
                pass
            else:
                distance = self.getFluxDistance(radius)
                self.setFlux(radius,distance)
        else:
            if radius is None:
                radius = self.getFluxRadius(distance)
                self.setFlux(radius,distance)
        
            else:
                self.setFlux(radius,distance)
        
        return self(x)
        
    def getFluxRadius(self,distance):
        """
        Determines the radius of a spherical blackbody at the specified distance
        assuming the flux is given by the model at the given temperature.
        """
        return (self.A*distance*distance/pi)**0.5
     
    def getFluxDistance(self,radius):
        """
        Determines the distance to a spherical blackbody of the specified radius
        assuming the flux is given by the model at the given temperature.
        """
        return (pi*radius*radius/self.A)**0.5
    
    def _getPeak(self):
        h,k = self.h,self.kb
        if 'wavelength' in self.unit:
            b = .28977685 #cm * K
            peakval = b/self.T/self._xscaling
        elif 'frequency' in self.unit:
            a=2.821439 #constant from optimizing BB function
            peakval=a/h*k*self.T/self._xscaling
        elif 'energy' in self.unit:
            raise NotImplementedError
        else:
            raise RuntimeError('Should never see this - bug in BB code')
        return self(peakval)
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
    
    peak=property(_getPeak,_setPeak)
    
    def wienDisplacementLaw(self,peakval):
        """
        Uses the Wien Displacement Law to calculate the temperature given a peak
        *input* location (wavelength, frequency, etc) or compute the peak
        location given the current temperature.
        
        :param peakval: 
            The peak value of the model to use to compute the new temperature.
        :type peakval: float or None
        
        :returns: 
            The temperature for the provided peak location or the peak location
            for the current temperature if `peakval` is None.
        """
        h,k = self.h,self.kb
        if self._phystype == 'wavelength':
            b = .28977685 #cm * K
            if peakval is None:
                out = b/self.T/self._xscaling
            else:
                out = b/peakval/self._xscaling
        elif self._phystype == 'frequency':
            a=2.821439 #constant from optimizing BB function
            if peakval is None:
                out = a*k*self.T/h/self._xscaling
            else:
                out = peakval*h/a/k/self._xscaling
        elif self._phystype == 'energy':
            a=2.821439 #constant from optimizing BB function
            if peakval is None:
                out = a*self.T/h/self._xscaling
            else:
                out = peakval*h/a/self._xscaling
        else:
            raise RuntimeError('Should never see this - bug in BB code')
        
        return out
    
    def _instanceSBLaw(self,T=None,area=1):
        if T is not None:
            self.T = T
        return BlackbodyModel.stephanBoltzmannLaw(self.T,area)*self._enscale
    
    @staticmethod
    def stephanBoltzmannLaw(T,area=1):
        """
        assumes cgs units
        """
            
        h,c,kb=BlackbodyModel.h,BlackbodyModel.c,BlackbodyModel.kb
        sigma = 2*pi**5*kb**4*h**-3*c**-2/15
        return area*sigma*T**4
    
    
class BurkertModel(FunctionModel1DAuto):
    r"""
    Burkert (1996) profile defined as:
    
    .. math::
        \frac{\rho_0 r_0^3}{(r+r_0)(r^2+r_0^2)}
        
    where :attr:`rho0` is the central density.
    
    """
    
    xaxisname = 'r'
    yaxisname = 'rho'
    
    def f(self,r,rho0=1,r0=1):
        return rho0*r0**3/((r+r0)*(r**2+r0**2))
    
    @property
    def rangehint(self):
        return 0,self.r0*3
    
    
class EinastoModel(FunctionModel1DAuto):
    """
    Einasto profile given by
    
    .. math::
        \\ln(\\rho(r)/A) = (-2/\\alpha)((r/r_s)^{\\alpha} - 1) . 
        
    :attr:`A` is the density where the log-slope is -2, and :attr:`rs` is the
    corresponding radius. The `alpha` parameter defaults to .17 as suggested by
    Navarro et al. 2010 
    
    .. note::
        The Einasto profile is is mathematically identical to the Sersic profile
        (:class:`SersicModel`), although the parameterization is different. By
        Convention, Sersic generally refers to a 2D surface brightness/surface 
        density profile, while Einasto is usually treated as a 3D density 
        profile.
    
    """
    xaxisname = 'r'
    yaxisname = 'rho'
    
    def f(self,r,A=1,rs=1,alpha=.17):
        return A*np.exp((-2/alpha)*((r/rs)**alpha - 1))
    
class ExponentialDiskModel(FunctionModel2DScalarAuto):
    """
    A disk with an exponential profile along both vertical and horizontal axes.
    The first coordinate is the horizontal/in-disk coordinate (scaled by `l`)
    wile the second is `z`. i.e.
    
    .. math::
        A e^{-|s/l|} e^{-|z/h|}
        
    for `pa` = 0 ; non-0 `pa` rotates the profile counter-clockwise by `pa`
    radians.
    """
    _fcoordsys='cartesian'
    def f(self,inarr,A=1,l=2,h=1,pa=0):
        s,z = inarr
        
        sinpa,cospa = np.sin(pa),np.cos(pa)
        sr = cospa*s+sinpa*z
        zr = -sinpa*s+cospa*z
        return A*np.exp(-np.abs(sr/l)-np.abs(zr/h))
    
    @property
    def rangehint(self):
        bigscale = max(self.h,self.l)
        return (-2*bigscale,2*bigscale,-2*bigscale,2*bigscale)
    
    
class InclinedDiskModel(FunctionModel2DScalarDeformedRadial):
    """
    Inclined Disk model -- identical to
    :class:`FunctionModel2DScalarDeformedRadial` but with inclination
    (:attr:`inc` and :attr:`incdeg`) and position angle (:attr:`pa` and
    :attr:`padeg`) in place of axis-ratios.
    """
    
    def __init__(self,inc=0,pa=0,degrees=True,**kwargs):
        super(InclinedDiskModel,self).__init__('sersic',**kwargs)
        
        self.n = 1
        
        if degrees:
            self.incdeg = inc
            self.padeg = pa
        else:
            self.inc = inc
            self.pa = pa
        
class RoundBulgeModel(FunctionModel2DScalarSeperable):
    """
    A bulge modeled as a radially symmetric sersic profile (by default the
    sersic index is kept fixed - to remove this behavior, set the fixedpars
    attribute to None)
    
    By default, the 'n' parameter does not vary when fitting.
    """
    
    fixedpars = ('n',)    
    def __init__(self,Ae=1,re=1,n=4):
        super(RoundBulgeModel,self).__init__('sersic')
        self.Ae = Ae
        self.re = re
        self.n = n
        
class ExponentialSechSqDiskModel(FunctionModel2DScalarAuto):
    """
    A disk that is exponential along the horizontal/in-disk (first) coordinate,
    and follows a :math:`{\\rm sech}^2(z)` profile along the vertical (second)
    coordinate. i.e.
    
    .. math::
        A e^{-|s/l|} {\\rm sech}^2 (z/h)
        
    for `pa` = 0 ; non-0 `pa` rotates the profile counter-clockwise by `pa`
    radians.
    """
    
    _fcoordsys='cartesian'
    def f(self,inarr,A=1,l=2,h=1,pa=0):
        s,z = inarr
        
        if pa == 0:
            sr,zr = s,z
        else:
            sinpa,cospa = np.sin(pa),np.cos(pa)
            sr = cospa*s+sinpa*z
            zr = -sinpa*s+cospa*z
            
        return A*np.exp(-np.abs(sr/l))*np.cosh(zr/h)**-2
    
    @property
    def rangehint(self):
        return (0,3*self.l,-3*self.h,3*self.h)
    
class GaussHermiteModel(FunctionModel1DAuto):
    """
    Model notation adapted from van der Marel et al 94
    
    hj3 are h3,h4,h5,... (e.g. N=len(hj3)+2 )
    """
    
    paramnames = 'h'
    paramoffsets = 3
    
    def f(self,v,A=1,v0=0,sig=1,*hj3):
        hj3arr = np.array(hj3,copy=False)
        hj3arr = hj3arr.reshape((hj3arr.size,1))
        w = (v-v0)/sig
        alpha = np.exp(-w**2/2)*(2*pi)**-0.5
        return A*alpha/sig*(1+np.sum(hj3arr*self._Hjs(w,len(hj3)+3,exclude=(0,1,2)),axis=0)) #sum start @ 3
    
    _Hpolys = None
    def _Hjs(self,w,N,exclude=None):
        """
        generates hermite polynomial arrays and evaluates them if w is not None
        """
        
        if self._Hpolys is None or N != len(self._Hpolys):
            from scipy.special import hermite
            self._Hpolys = [hermite(i) for i in range(N)]
            
        if w is not None:
            warr = np.array(w,copy=False).ravel()
            if exclude is None:
                return np.array([H(warr) for H in self._Hpolys ])
            else:
                return np.array([H(warr) for i,H in enumerate(self._Hpolys) if i not in exclude])
    
    @property
    def rangehint(self):
        return self.v0-self.sig*4,self.v0+self.sig*4
    
    def gaussHermiteMoment(self,l,f=None,lower=-np.inf,upper=np.inf):
        """
        compute the requested moment on the supplied function, or this
        object if f is None
        """
        from scipy.integrate import quad
        
        if int(l)!=l:
            raise ValueError('moment specifier must be an integer')
        l = int(l)
        
        self._Hjs(None,len(self.params)) 
        def gHJac(v,A,v0,sig,*hj3):
            w = (v-v0)/sig
            alpha = np.exp(-w**2/2)*(2*pi)**-0.5
            return alpha*self._Hpolys[l](w)
        
        if f is None:
            f = self
        
        intnorm = quad(lambda x:self._Hpolys[l](x)*self._Hpolys[l](x)*np.exp(-x*x),lower,upper)[0]/(2*pi)
        #return self.integrate(lower,upper,jac=gHJac)/self.A/intnorm
        return quad(lambda x:f(x)*gHJac(x,*self.parvals),lower,upper)[0]/self.A/intnorm
    
class HernquistModel(FunctionModel1DAuto):
    """
    Hernquist (1990) profile defined as:
    
    .. math::
        \\frac{A r_0}{2 \\pi r (r+r_0)^3}
    
    where :attr:`A` is the total mass enclosed. Note that :attr:`r0` does not
    enclose half the mass - the radius enclosing half the mass is 
    :math:`r_h = \\frac{r_0}{\\sqrt{2}-1}` .
    
    .. note::
        This form is equivalent to an :class:`AlphaBetaGammaModel` with
        :math:`(\\alpha,\\beta,\\gamma) = (1,4,1)` , but with a slightly
        different overall normalization.
    
    """
    
    xaxisname = 'r'
    yaxisname = 'rho'
    
    def f(self,r,A=1,r0=1):
        return (A/2/pi)*(r0/r)*(r+r0)**-3
    
    @property
    def rangehint(self):
        return 0,self.r0*3
    

class JaffeModel(FunctionModel1DAuto):
    """
    Jaffe (1983) profile as defined in Binney & Tremaine as:
    
    .. math::
        \\frac{A}{4 \\pi} \\frac{r_j}{r^2 (r+r_j)^2)}
    
    where :attr:`A` is the total mass enclosed. :attr:`rj` is the radius that
    encloses half the mass.
    
    .. note::
        This form is equivalent to an :class:`AlphaBetaGammaModel` with
        :math:`(\\alpha,\\beta,\\gamma) = (1,4,2)` , but with a slightly
        different overall normalization. 
    
    """
    
    xaxisname = 'r'
    yaxisname = 'rho'
    
    def f(self,r,A=1,rj=1):
        return (A/4/pi)*rj*r**-2*(r+rj)**-2
    
    @property
    def rangehint(self):
        return 0,self.rj*3
    
class NFWModel(FunctionModel1DAuto):
    """
    A Navarro, Frenk, and White 1996 profile -- the canonical
    :math:`\\Lambda {\\rm CDM}` dark matter halo profile:
    
    .. math::
        \\rho(r) = \\frac{\\rho_0}{r/r_c (1+r/r_c)^2} .
    
    Where unspecified, units are r in kpc and rho in Msun pc^-3, although units can
    always be arbitrarily scaled by rescaling `rho0`.
    
    .. note::
        This form is equivalent to an :class:`AlphaBetaGammaModel` with
        :math:`(\\alpha,\\beta,\\gamma) = (1,3,1)` , but with a slightly
        different overall normalization. This class also has a number of
        helper methods.
    
    """
    xaxisname = 'r'
    yaxisname = 'rho'
        
    def f(self,r,rho0=1,rc=1):
        x = r/rc
        return rho0/(x*(1+x)**2)
        #return rho0*rc*rc*rc*((r+rc)**(-2))*(r**-1)
    
    @property
    def rangehint(self):
        return self.rc/1000,1000*self.rc
    
    def toAlphaBetaGamma(self):
        """
        returns a new :class:`AlphaBetaGammaModel` based on this model's
        parameters.
        """
        mod = AlphaBetaGammaModel(A=self.rho0,rs=self.rc,alpha=1,beta=3,gamma=1)
        if self.data is None:
            mod.data = None
        else:
            mod.data = (self.data[0].copy(),
                        self.data[1].copy(),
                        None if self.data[2] is None else self.data[2].copy())
        return mod
    
    def integrateSpherical(self,lower,upper,method=None,**kwargs):
        """
        NFW Has an analytic form for the spherical integral - if the lower 
        is not 0 or or if `method` is not None, this function will
        fall back to the numerical version.
        """    
        if method is not None or np.any(lower!=0):
            return FunctionModel1D.integrate(self,lower,upper,method,**kwargs)
            
        x = upper/self.rc
        return 4*pi*self.rho0*self.rc**3*(np.log(1+x)-x/(1+x))
        
    def getV(self,r,**kwargs):
        r"""
        Computes the circular velocity of the halo at a given radius.  i.e.
        :math:`v = \sqrt{\frac{G M(<r)}{r}}`.
        
        :param r: 
            The radius (or radii) at which to compute the circular velocity, in
            kpc.
            
        Additional keyword arguments are passed into :meth:`integrateSpherical`.
        
        :returns: 
            The circular velocity at `r` (if `r` is a scalar, this will return a
            scalar, otherwise an array)
        
        """
        from .constants import GMspc
        
        r = np.array(r,copy=False)
        
        call = self.getCall()
        try:
            M = self.integrateSpherical(0,r,**kwargs)
        finally:
            if call is None:
                self.setCall(None)
            else:
                self.setCall(*call)
        
        return ((GMspc/1000)*M/r)**0.5 #GMspc/1000 converts from pc to kpc
    
    #numerically computed derivative of spherical integral divided by x
    #_xmin = scipy.optimize.fmin(finv,1,xtol=1e-16,ftol=1e-16,disp=1)[0]    
    _xmin = 2.16258158706460983485655366960326457351
    def getVmax(self):
        """
        Computes the maximum circular velocity of this profile.
        
        :returns: 
            (vmax,rvmax) where vmax is the maximum circular velocity in km/s and
            rvmax is the radius at which the circular velocity is maximized.
        
        """
        
        rmin = self.rc*NFWModel._xmin
        return self.getV(rmin),rmin
        
    def setC(self,c,Rvir=None,Mvir=None,z=0):
        """
        Sets the model parameters to match a given concentration given virial
        radius and mass.
        
        :param c: concentration = rvir/rc
        :type c: scalar
        :param Rvir:
            virial radius to assume in kpc - if None, will be inferred from
            `Mvir`, `c` and the virial overdensity (see `deltavir`).
        :type Rvir: scalar or None
        :param Mvir: 
            virial mass to assume in Msun - if None, will be inferred from
            `Rvir`, `c`, and the virial overdensity (see `deltavir`), or if
            `Rvir` is None, from the current virial mass.
        :type Mvir: scalar or None
        :param z: 
            redshift at which to compute virial overdensity if `Rvir` or `Mvir`
            are None. If both are given, this is ignored.
        :type z: scalar
        
        """
        from .constants import get_cosmology
#        if Rvir is None and Mvir is None:
#            raise ValueError('need to specify Rvir or Mvir')
#        elif Rvir is None:
#            Rvir = self.Mvir_to_Rvir(Mvir,z=0)
#        elif Mvir is None:
#            Mvir = self.Rvir_to_Mvir(Rvir,z=0)
            
        if Rvir is None:
            if Mvir is None:
                Mvir = self.getMv()
                
            rhov = get_cosmology().rhoC(z,'cosmological')*1e-18*self.deltavir(z)
            Rvir = 1e-3*(3*Mvir/(4*pi*rhov))**(1/3)
            
        elif Mvir is None:
            rhov = get_cosmology().rhoC(z,'cosmological')*1e-18*self.deltavir(z)
            Mvir = (4*pi*(Rvir*1e3)**3/3)*rhov
        else: #both are specified, implying a particular deltavir
            self._c = c
            
        self.rc = Rvir/c
        self.rho0 = 1
        a0 = self.integrateSpherical(0,Rvir)
        self.rho0 = Mvir/a0
        
    def getC(self,z=0):
        """
        Returns the concentration of this profile(rvir/rc) for a given redshift.
        
        
        .. note::
            If :meth:`setC` has been called with an explicitly provided Rvir and
            Mvir, the concentration will be saved and this will return the saved
            concentration instead of one computed for a given redshift. If this
            is not desired, delete the :attr:`_c` attribute after calling
            :meth:`setC` .
        
        """
        if hasattr(self,'_c'):
            return self._c
        else:
            return self.getRv(z)/self.rc
        
    def getRhoMean(self,r):
        """
        Computes the mean density within the requested radius, i.e.
        
        .. math::
            \\rho_m = M(<r)/(4/3 pi r^3)
            
        density in units of Msun pc^-3 assuming `r` in kpc
        """
        vol = 4*pi*(r*1e3)**3/3
        return self.integrateSpherical(0,r)/vol
    
    def deltavir(self,z=0):
        """
        Returns the virial overdensity at a given redshift.  
        
        This method should be overridden if a particular definition of virial
        radius is desired. E.g. to use r200 instead of rvir, do::
        
            nfwmodel.deltavir = lambda z:200
        
        
        :param z: redshift at which to compute virial overdensity
        :type z: scalar
        
        :returns: virial overdensity        
        
        """
        from .constants import get_cosmology
        
        return get_cosmology().deltavir(z)
            
    def getRv(self,z=0):
        """
        Get the virial radius at a given redshift, with `deltavir` choosing the
        type of virial radius - if 'fromcosmology' it will be computed from the
        cosmology, otherwise, it should be a multiple of the critical density
        
        units in kpc for mass in Msun
        """
        from scipy.optimize import newton
        from .constants import get_cosmology,Ms
        
        cosmo = get_cosmology()
        overden = self.deltavir(z)
        
        try:
            rhov = self.deltavir(z)*cosmo.rhoC(z,'cosmological')*1e-18 
            # *1e-18  does Mpc^-3->pc^-3
        except:
            raise ValueError('current cosmology does not support critical density')
        
        oldcall = self.getCall()
        try:
            self.setCall('getRhoMean')
            return self.inv(rhov,self.rc)
        finally:
            if oldcall is None:
                self.setCall(None)
            else:
                self.setCall(*oldcall)
                
        
    
    def getMv(self,z=0):
        """
        Gets the mass within the virial radius (for a particular redshift)
        """
        rv = self.getRv(z)
        return self.integrateSpherical(0,rv)
    
    
    #Below are static scalings that are approximations 
    
    @staticmethod
    def Rvir_to_Mvir(Rvir,z=0):
        """
        from Bullock '01
        M_sun,kpc
        
        .. warning::
            These scalings are approximations.
            
        :param float Rvir: Virial Radius in kpc
        :param float z: Redshift
        
        :returns: Virial mass in Msun
            
        """
        #Mvir = 10^12 Msun/h [Omega_0 Delta(z)/97.2] [Rvir(1+z)/203.4 kpc/h]^3
        from .constants import get_cosmology
        c = get_cosmology()
        
        return 1e12/c.h*(c.omega*c.deltavir(z)/97.2)*(Rvir*c.h*(1+z)/203.4)**3
    
    @staticmethod
    def Mvir_to_Rvir(Mvir,z=0):
        """
        from Bullock '01
        M_sun,kpc
        
        .. warning::
            These scalings are approximations.
            
        :param float Mvir: Virial Mass in Msun
        :param float z: Redshift
        
        :returns: Virial radius in kpc
        
        """
        #Rvir = 203.4 kpc/h [Omega_0 Delta(z)/97.2]^-1/3 [Mvir/10^12 h^-1 Msun]^1/3 (1+z)^-1 ~= 300 kpc [Mvir/10^12 h^-1 Msun]^1/3 (1+z)^-1
        from .constants import get_cosmology
        c = get_cosmology()
        
        return 203.4/c.h*(c.omega*c.deltavir(z)/97.2)**(-1/3)*(Mvir*c.h/1e12)**(1/3)/(1+z)
    
    @staticmethod
    def Mvir_to_Vvir(Mvir,z=0):
        """
        from Bullock '01
        km/s,M_sun
        
        .. warning::
            These scalings are approximations.
            
        :param float Mvir: Virial Mass in Msun
        :param float z: Redshift
        
        :returns: Virial Velocity in km/s
        
        """
        #Vvir = 143.8 km/s [Omega_0 Delta(z)/97.2]^1/6 [Mvir/10^12 Msun/h]^1/3 (1+z)^1/2
        from .constants import get_cosmology
        c = get_cosmology()
        
        return 143.8*(c.omega*c.deltavir(z)/97.2)**(1/6)*(Mvir*c.h/1e12)**(1/3)*(1+z)**0.5
    
    @staticmethod
    def Vvir_to_Mvir(Vvir,z=0):
        """
        from Bullock '01
        km/s,M_sun
        
        .. warning::
            These scalings are approximations.
            
        :param float Vvir: Virial velocity in km/s
        :param float z: Redshift
        
        :returns: Virial mass in Msun
            
        """
        from .constants import get_cosmology
        c = get_cosmology()
        
        return ((c.omega*c.deltavir(z)/97.2)**-0.5*(1+z)**-1.5*1e12*(Vvir/143.8)**3)/c.h
    
    @staticmethod
    def Vvir_to_Vmax(Vvir):
        """
        from Bullock '01
        good from 80-1200 km/s in vvir
        
        .. warning::
            These scalings are approximations.
            
        :param float Vvir: Virial velocity in km/s
        
        :returns: Virial velocity in km/s
            
        """
        return (Vvir/0.468)**(1/1.1)
    
    @staticmethod
    def Vmax_to_Vvir(Vmax):
        """
        from Maller & Bullock '04
        good from 80-1200 km/s in vvir
        
        .. warning::
            These scalings are approximations.
            
        :param float Vmax: Maximum circular Velocity in km/s 
        :param float z: Redshift
        
        :returns: Virial Velocity in km/s
            
        """
        return (0.468)*Vmax**1.1
    
    @staticmethod
    def Mvir_to_Cvir(Mvir,z=0):
        """
        from Maller & Bullock '04
        M_sun
        
        .. warning::
            These scalings are approximations.
            
        :param float Mvir: Virial Mass in Msun
        :param float z: Redshift
        
        :returns: Concentration parameter (rvir/rc)
            
        """
        from .constants import H0
        
        cvir = 9.6*(Mvir*(H0/72.)/1e13)**-.13*(1+z)**-1
        return cvir if cvir>4 else 4.0
    
    @staticmethod
    def Cvir_to_Mvir(Cvir,z=0):
        """
        from Maller & Bullock '04
        M_sun
        
        .. warning::
            These scalings are approximations.
            
        :param float Cvir: Concentration parameter rvir/rc
        :param float z: Redshift
        
        :returns: Virial mass in Msun
            
        """
        return 1e13/(c.h/.72)*((Cvir/9.6)*(1+z))**(-1.0/.13)
    
    @staticmethod
    def Mvir_to_Vmax(Mvir,z=0):
        """
        Convert virial mass to vmax using the scalings here
        
        .. warning::
            These scalings are approximations.
        
        :param float Mvir: Virial Mass in Msun
        :param float z: Redshift
        
        :returns: Maximum circular velocity in km/s
            
        """
        return NFWModel.Vvir_to_Vmax(NFWModel.Mvir_to_Vvir(Mvir,z))
    
    @staticmethod
    def Vmax_to_Mvir(Vmax,z=0):
        """
        Convert vmax to virial mass using the scalings here
        
        .. warning::
            These scalings are approximations.
            
        :param float Vmax: Maximum circular velocity in km/s
        :param float z: Redshift
        
        :returns: Virial mass in Msun
            
        """
        return NFWModel.Vvir_to_Mvir(NFWModel.Vmax_to_Vvir(Vmax),z)
    
    @staticmethod
    def Vmax_to_Rvir(Vmax,z=0):
        """
        Convert vmax to virial radius using the scalings here.
        
        .. warning::
            These scalings are approximations.
            
        :param float Vmax: Maximum circular velocity in km/s
        :param float z: Redshift
        
        :returns: Virial radius in kpc
            
        """
        return NFWModel.Mvir_to_Rvir(NFWModel.Vmax_to_Mvir(Vmax),z)
    
    @staticmethod
    def RvirMvir_to_Vvir(Rvir,Mvir):
        """
        Computes virial velocity from virial radius and virial mass.
        
        :param float Rvir: Virial Radius in kpc
        :param float Mvir: Virial Mass in Msun
        
        :returns: Virial Velocity in km/s
        """
        from .constants import G
        
        return (G*Mvir/Rvir)**0.5
    
    @staticmethod
    def Vvir_to_Tvir(Vvir,mu=.59):
        r"""
        Compute the Virial temperature for a given virial radius.
        
        From Barkana & Loeb 2001:
        :math:`T_{\rm vir} = \frac{\mu m_p V_c^2}{2 k_B}`
        
        :math:`\mu` is the mean molecular weight, which is 0.59 for fully
        ionized primordial gas, .61 for fully ionized hydrogen but singly
        ionized helium, and 1.22 neutral. :math:`m_p` is the proton mass,
        :math:`k_B` is the Boltzmann constant, and :math:`V_c` is the virial
        velocity.
        
        
        :param float Vvir: Virial velocity in km/s
        :param float mu: Mean molecular weight
        
        :returns: Virial Temperatur in Kelvin
            
        
        """
        from constants import mp,kb
        
        return  mu*mp*Vvir**2/(2*kb)*1e10 #(km/s)^2 -> (cm/2)^2
    
    @staticmethod
    def create_Mvir(Mvir,z=0):
        """
        Generates a new NFWModel with the given Mvir using the
        :meth:`Mvir_to_Cvir` function to get the scaling.
        
        .. warning::
            These scalings are approximations.
            
        :param float Mvir: Virial Mass in Msun
        :param float z: Redshift at which model is valid
        
        :returns: 
            A :class:`NFWModel` object with the provided Mvir and concentration
            set by the approximate scalings here.
        
        """
        m = NFWModel()
        Cvir = m.Mvir_to_Cvir(Mvir,z)
        Rvir = m.Mvir_to_Rvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
    @staticmethod
    def create_Rvir(Rvir,z=0):
        """
        Generates a new NFWModel with the given Mvir using the Mvir_to_Cvir and
        :meth:`Rvir_to_Mvir` functions to get the scaling.
        
        .. warning::
            These scalings are approximations.
            
        :param float Rvir: Virial Radius in kpc
        :param float z: Redshift at which the model is valid 
        
        :returns: 
            A :class:`NFWModel` object with the provided Rvir and concentration
            set by the approximate scalings here.
        
        """
        m = NFWModel()
        Mvir = m.Rvir_to_Mvir(Rvir,z)
        Cvir = m.Mvir_to_Cvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
    @staticmethod
    def create_Cvir(Cvir,z=0):
        """
        Generates a new NFWModel with the given Mvir using the
        :meth:`Cvir_to_Mvir` function to get the scaling.
        
        .. warning::
            These scalings are approximations.
            
        :param float Cvir: Concentration parameter (rvir/rc)
        :param float z: Redshift at which the model is valid
        
        :returns: 
            A :class:`NFWModel` object with the provided Cvir and normalization
            set by the approximate scalings here.
        
        """
        m = NFWModel()
        Mvir = m.Cvir_to_Mvir(Cvir,z)
        Rvir = m.Mvir_to_Rvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
    @staticmethod
    def create_VmaxRvmax(vmax,rvmax):
        """
        Generates a new NFWModel with the given Vmax and R_Vmax.
        
        :param float vmax: maximum circular velocity in km/s.
        :param float rvmax: radius of maximum circular velocity in kpc.
        
        See Bullock et al. 2001 for reflated reference.
        
        :returns: 
            A :class:`NFWModel` object matching the supplied `vmax` and `rvmax`
        
        """
        from scipy.optimize import fmin
        
        #generate an approximate "best guess" from the scalings at z=0
        m = NFWModel.create_Mvir(NFWModel.Vmax_to_Mvir(vmax),z=0)
        
        #from Bullock+ 01
        m.rc = rvmax/2.16
        
        def toopt(v,vmaxwanted,model):
            model.rho0 = v[0]
            return (model.getVmax()[0]-vmaxwanted)**2
        
        fmin(toopt,(m.rho0,),(vmax,m),disp=0)
        #rho0 is now set correctly
        
        return m
    
class NFWProjectedModel(FunctionModel1DAuto):
    
    def f(self,R,rc=1,sig0=1):
        x = R/rc
        xsqm1 = x*x - 1
        
        Cinv = np.arccos(1/x)
        nanmsk = np.isnan(Cinv)
        Cinv[nanmsk] = np.arccosh(1/x[nanmsk])
        
        xterm = (1-np.abs(xsqm1)**-0.5*Cinv)/xsqm1
        xterm[x==1] = 1/3
                
        return sig0*xterm/(2*pi*rc**2)
        #sig0=Mv*g
        
    def integrateCircular(self,lower,upper,method=None,**kwargs):
        """
        NFW Has an analytic form for the spherical integral - if the lower 
        is not 0 or or if the keyword 'numerical' is True, this function will
        fall back to FunctionModel1D.integrateCircular 
        """        
        if method is not None or np.any(lower!=0):
            return FunctionModel1D.integrate(self,lower,upper,method,**kwargs)
        
        x = upper/self.rc
        xterm = None
        
        if x > 1:
            Cinv = np.arccos(1/x)
        elif x < 1:
            Cinv = np.arccosh(1/x)
        else:
            xterm = 1-np.log(2)
        
        if xterm is None:
            xterm = Cinv*np.abs(x*x - 1)**-0.5 + np.log(x/2)
                
        return self.sig0*xterm
        
    @property
    def rangehint(self):
        return self.rc/1000,1000*self.rc
    
class MoffatModel(FunctionModel1DAuto):
    """
    Moffat function given by:
    
    .. math::
        A \\frac{(\\beta-1}{\\pi \\alpha^2} \\left(1+\\left(\\frac{r}{\\alpha}\\right)^2\\right)^{-\\beta}
    """
    def f(self,r,A=1,alpha=1,beta=4.765):
        roa=r/alpha
        return A*(beta-1)/(pi*alpha**2)*(1+roa*roa)**-beta
    
    def _getFWHM(self):
        return self.alpha*2*(2**(1/self.beta)-1)**0.5
    def _setFWHM(self,val):
        self.alpha = val*(4*(2**(1/self.beta)-1))**-0.5
    FWHM = property(_getFWHM,_setFWHM,doc='Full Width at Half Maximum for this model - setting changes alpha for fixed beta')
    
    
    @property
    def rangehint(self):
        return(-self.alpha,self.alpha)
        


class PlummerModel(FunctionModel1DAuto):
    """
    Plummer model of the form
    
    .. math::
        \\frac{3M}{4 \\pi r_p^3} (1+(r/r_p)^2)^{-5/2}
    
    """
    
    xaxisname = 'r'
    
    def f(self,r,rp=1.,M=1.):
        return 3*M/(4.*pi*rp**3)*(1+(r/rp)**2)**-2.5
    
    @property
    def rangehint(self):
        return 0,self.rp*2

class King2DrModel(FunctionModel1DAuto):  
    """
    2D (projected/surface brightness) King model of the form:
    
    .. math::
        f(r) = A r_c^2 (\\frac{1}{\\sqrt{r^2+r_c^2}} - \\frac{1}{\\sqrt{r_t^2+r_c^2}})^2
    
    .. seealso:: King (1966) AJ Vol. 71, p. 64
    """
      
    xaxisname = 'r'
    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        return A*rcsq*((r*r+rcsq)**-0.5 - (rt*rt+rcsq)**-0.5)**2
    
    @property
    def rangehint(self):
        return 0,self.rt
    
class King3DrModel(FunctionModel1DAuto):
    """
    3D (deprojected) King model of the form:
    
    .. math::
        f(r) = A/(z^2 \\pi r_c) ((1+(r_t/r_c)^2)^{-3/2}) \\arccos(z)/(z-\\sqrt{1-z^2})
    
    """
    
    xaxisname = 'r'
    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        z = ((r*r+rcsq)**0.5) * ((rt*rt+rcsq)**-0.5)
        res = (A/z/z/pi/rc)*((1+rt*rt/rcsq)**-1.5)*(np.arccos(z)/z-(1-z*z)**0.5)
        res[r>=rt] = 0
        return res
    
    @property
    def rangehint(self):
        return 0,self.rt

class SchechterMagModel(FunctionModel1DAuto):
    """
    The Schechter Function, commonly used to fit the luminosity function of
    galaxies, in magnitude form:
    
    .. math::
        \\Phi(M) = \\phi^* \\frac{2}{5} \\ln(10) \\left[10^{\\frac{2}{5} (M_*-M)}\\right]^{\\alpha+1} e^{-10^{\\frac{2}{5} (M_*-M)}}
    
    .. seealso:: :class:`SchechterLumModel`, Schechter 1976, ApJ 203, 297 
    """
    
    xaxisname = 'M'
    yaxisname = 'phi'
    
    from math import log
    _frontfactor = log(10)*0.4
    del log #hide this so as not to clutter the namespace
    
    def f(self,M,Mstar=-20.2,alpha=-1,phistar=1.0):
        from math import log #single-variable version
        x=10**(0.4*(Mstar-M))
        return SchechterMagModel._frontfactor*phistar*(x**(1+alpha))*np.exp(-x)

    def derivative(self,M,dx=None):
        """
        Compute Schechter derivative for magnitude form. if `dx` is not None,
        will fallback on the numerical version.
        """
        if dx is not None:
            return FunctionModel1D.derivative(self,M,dx)
        
        a = self.alpha
        x = 10**(0.4*(self.Mstar-M))
        return -SchechterMagModel._frontfactor**2*self.phistar*np.exp(-x)*x**a*(a-x+1)
    
    def integrate(self,lower,upper,method=None,**kwargs):
        """
        Analytically compute Schechter integral for magnitude form using
        incomplete gamma functions. If `method` is not None, numerical
        integration will be used. The gamma functions break down for alpha<=-1,
        so numerical is used if that is the case.
        """
        if self.alpha<=-1:
            method = True #use default numerical method, because gamma functions fail for alpha<=-2
        if method is not None:
            return FunctionModel1D.integrate(self,lower,upper,method,**kwargs)
        
        from scipy.special import gamma,gammainc,gammaincc
    
        s = self.alpha+1
        u = 10**(0.4*(self.Mstar-upper))
        l = 10**(0.4*(self.Mstar-lower))
        
        if upper==np.inf and lower<=0:
            I = gamma(s)
        elif upper==np.inf:
            I = gammaincc(s,l)*gamma(s)
        elif lower==0:
            I = gammainc(s,u)*gamma(s)
        else:
            I = (gammainc(s,u) - gammainc(s,l))*gamma(s)

        return -self.phistar*I
    
    @property
    def rangehint(self):
        return self.Mstar-3,self.Mstar+3
    
class SchechterLumModel(FunctionModel1DAuto):
    """
    The Schechter Function, commonly used to fit the luminosity function of 
    galaxies, in luminosity form:
    
    .. math::
        \\phi(L) = \\frac{\\phi^*}{L_*} \\left(\\frac{L}{L_*}\\right)^\\alpha e^{-L/L_*}   
    
    .. seealso:: :class:`SchechterMagModel`, Schechter 1976, ApJ 203, 297 

    """
    
    xaxisname = 'L'
    yaxisname = 'phi'
    
    def f(self,L,Lstar=1e10,alpha=-1.0,phistar=1.0):
        x = L/Lstar
        return phistar*(x**alpha)*np.exp(-x)/Lstar
    
    def derivative(self,L,dx=None):
        """
        Compute Schechter derivative.  if `dx` is not None, will fallback on the
        numerical version.
        """
        if dx is not None:
            return FunctionModel1D.derivative(self,L,dx)
        
        a = self.alpha
        x = L/self.Lstar
        return self.phistar*self.Lstar**-2*(a-x)*np.exp(-x)*x**(a-1)
    
    def integrate(self,lower,upper,method=None,**kwargs):
        """
        Analytically Compute Schechter integral using incomplete gamma
        functions. If `method` is not None, numerical integration will be used.
        The gamma functions break down for alpha<=-1, so numerical is used if
        that is the case.
        """
        if self.alpha<=-1:
            method = True #use default numerical method, because gamma functions fail for alpha<=-1
        if method is not None:
            return FunctionModel1D.integrate(self,lower,upper,method,**kwargs)
        
        
        from scipy.special import gamma,gammainc,gammaincc
            
        s = self.alpha+1
        u = upper/self.Lstar
        l = lower/self.Lstar
        
        if upper==np.inf and lower<=0:
            I = gamma(s)
        elif upper==np.inf:
            I = gammaincc(s,l)*gamma(s)
        elif lower==0:
            I = gammainc(s,u)*gamma(s)
        else:
            I = (gammainc(s,u) - gammainc(s,l))*gamma(s)

        return self.phistar*I
    
    @property
    def rangehint(self):
        return self.Lstar/3,self.Lstar*3
        


class SersicModel(FunctionModel1DAuto):
    """
    Sersic surface brightness profile:
    
    .. math::
        A_e e^{-b_n[(R/R_e)^{1/n}-1]}
    
    Ae is the value at the effective radius re
    
    .. note::
        The Sersic profile is is mathematically identical to the Einasto profile
        (:class:`EinastoModel`), although the parameterization is different. By
        Convention, Sersic generally refers to a 2D surface brightness/surface 
        density profile, while Einasto is usually treated as a 3D density 
        profile.
        
    """
    xaxisname = 'r'
    
    def f(self,r,Ae=1,re=1,n=4):
        #return EinastoModel.f(self,r,A,rs,1/n)
        #return A*np.exp(-(r/rs)**(1/n))
        return Ae*np.exp(-self.getBn()*((r/re)**(1.0/n)-1))
    
    @property
    def rangehint(self):
        return 0,2*self.re
    
    def _getA0(self):
        return self.f(0,self.Ae,self.re,self.n)
    def _setA0(self,val):
        self.Ae *= val/self.f(0,self.Ae,self.re,self.n)
    A0 = property(_getA0,_setA0,doc='value at r=0')
    
    _exactbn = True
    _bncache = {}
    def getBn(self):
        """
        Computes :math:`b_n` for the current sersic index. If the
        :attr:`SersicModel.exactbn` attribute is True, this is calculated using
        incomplete gamma functions (:func:`bn_exact`), otherwise, it is
        estimated based on MacArthur, Courteau, and Holtzman 2003
        (:func:`bn_estimate`).
        
        """
        n = self.n
        if n in SersicModel._bncache:
            return SersicModel._bncache[n]
        else:
            if SersicModel._exactbn:
                bn = SersicModel.bn_exact(n)
            else:
                bn = SersicModel.bn_estimate(n)
            SersicModel._bncache[n] = bn
            return bn
        
    @staticmethod
    def exactBn(val=None):
        """
        Sets whether the exact Bn calculation is used, or the estimate (see
        :meth:`getBn`).
        
        :param val: 
            If None, nothing will be set, otherwise, if True, the exact Bn
            computation will be used, or if False, the estimate will be used.
            
        :returns: 
            True if the exact computation is being used, False if the estimate. 
        """
        if val is not None:
            SersicModel._exactbn = bool(val)
            SersicModel._bncache = {}
        return SersicModel._exactbn
    
    @staticmethod
    def bn_exact(n):
        """
        Computes :math:`b_n` exactly for the current sersic index, using
        incomplete gamma functions.
        """
        from scipy.special import gammaincinv
        
        n = float(n) #sometimes 0d array gets passed in
        return gammaincinv(2*n,0.5)
    
    _bnpoly1=np.poly1d([-2194697/30690717750,131/1148175,46/25515,4/405,-1/3])
    _bnpoly2=np.poly1d([13.43,-19.67,10.95,-0.8902,0.01945])
    @staticmethod
    def bn_estimate(n):
        """
        bn is used to get the half-light radius.  If n is 
        None, the current n parameter will be used
        
        The form is a fit from MacArthur, Courteau, and Holtzman 2003 
        and is claimed to be good to ~O(10^-5)
        """
        n = float(n) #sometimes 0d array gets passed in
        
        return (2*n+SersicModel._bnpoly1(1/n)) if n>0.36 else SersicModel._bnpoly2(n)
        
    def sbfit(self,r,sb,zpt=0,**kwargs):
        """
        fit surface brightness using the SersicModel
        
        r is the radial value,sb is surface brightness, zpt is the zero point
        of the magnitude scale, and kwargs go into fitdata
        """
        flux = 10**((zpt-sb)/2.5)
        return self.fitData(r,flux,**kwargs)
        
    def sbplot(self,lower=None,upper=None,data=None,n=100,zpt=0,clf=True):
        """
        plots the surface brightness for this flux-based SersicModel.  arguments
        are like fitDat
        """
        from matplotlib import pyplot as plt
        
        if data is None and (lower is None or upper is None):
            raise ValueError('need data for limits or lower/upper')
        if data is not None:
            if upper is None:
                upper = np.max(data[0])
            if lower is None:
                lower = np.min(data[0])
        
        if clf:
            plt.clf()
        
        x = np.linspace(lower,upper,n)
        plt.plot(x,zpt-2.5*np.log10(self(x)))
        if data:
            skwargs={'c':'r'}
            plt.scatter(*data,**skwargs)
        
        plt.ylim(*reversed(plt.ylim()))
        plt.xlim(lower,upper)
    
class DeVaucouleursModel(SersicModel):
    """
    Sersic model with n=4.
    """
    
    xaxisname = 'r'
    
    def f(self,r,Ae=1,re=1):
        return SersicModel.f(self,r,Ae,re,4)
    
#register everything in this module
from inspect import isclass
for o in locals().values():
    if isclass(o) and not o.__name__.startswith('_') and issubclass(o,ParametricModel):
        if 'FunctionModel' not in o.__name__ and 'CompositeModel' not in o.__name__:
            register_model(o)

del ABCMeta,abstractmethod,abstractproperty,isclass,o #clean up namespace
