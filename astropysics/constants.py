#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""

=========
constants
=========

The :mod:`constants` module contains attributes storing physical constants and
conversion factors. Most of these are at the package level and should be
imported as::

    from astropysics.constants import c,G

It also includes classes representing various cosmologies that are used to
derive relevant cosmological parameters. The current default is the
:class:`WMAP7Cosmology`, based on the LCDM cosmology with parameters favored by
`WMAP7 <http://lambda.gsfc.nasa.gov/product/map/dr4/parameters.cfm>`_

.. todo:: examples


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.constants
   :parts: 1

Module API
----------

"""

from __future__ import division,with_statement
from math import pi
import numpy as np

#<--------------------------------Constants------------------------------------>
#all in cgs including gaussian esu for charge
unit_system='cgs'

G = 6.673e-8 #gravitational constant
mp = 1.67262171e-24 #proton mass
me = 9.1093897e-28 #electron mass
e = 4.8032068e-10 #electron charge
Ms = 1.9891e33 #solar mass
Mj = 1.8986e30 #jupiter mass
Me = 5.9742e27 #earth mass
Rs = 6.96e10 #solar radius
Rj = 7.1492e9 #jupiter radius
Re = 6.371e8 #earth radius
Lsun = 3.839e33 #solar luminosity
kb = 1.3807e-16 #boltzmann's constant
Rb = 8.314472e7 #gas constant
c = 2.99792458e10 #speed of light - exact
h = 6.626068E-27 #planck's constant
hbar = h/2/pi
g0 = 980.665 #mean earth gravitational acceleration at sea level

#<-------------------------------Conversions----------------------------------->
ergperev = 1.60217646e-12
secperday = 86400 #IAU standard
secperyr = 365.25*secperday#3.15576926e7
secpergyr = secperyr*1e9
cmperpc = 3.08568025e18
pcpercm = 1.0/cmperpc
lyperpc = 3.26
pcperly = 1.0/lyperpc
cmperau = 1.49597887e13
aupercm = 1.0/cmperau


#TODO:remove this once these are no longer in use in favor of spec.HasSpectrumUnits?
def flambda_to_fnu_l(flambda,lamb):
    return flambda*lamb*lamb/c

def fnu_to_flambda_l(fnu,lamb):
    return fnu*c/lamb/lamb

def flambda_to_fnu_n(flambda,nu):
    return flambda*c/nu/nu

def fnu_to_flambda_n(fnu,nu):
    return fnu*nu*nu/c


#<--------------------------------Cosmology------------------------------------>
class Cosmology(object):
    """
    This class represents a cosmology and should be subclassed
    
    all cosmologies should have a hubble constant (H0) in km/s/Mpc
    
    Cosmologies should also include a sequence called "_params_" with a list of
    strings that specify the names of the cosmological parameters associated to
    be exported to the constants module.
    
    Error bars for parameters can be specified as "<paramname>_err" as a tuple
    (lowererr,uppererr)
    """
    _params_=('H0',)
    _autoupdate = False
    
    H0 = 0
    
    def __init__(self,*args,**kwargs):
        ps = self._params_
        
    h = property(lambda self:self.H0/100.0)
    h7 = h70 = property(lambda self:self.H0/70.0)

    __params_cache=None
    def _getParams(self):
        if self.__params_cache is None:
            import inspect
            pars = [cls._params_ for cls in inspect.getmro(self.__class__) if 
                    hasattr(cls,'_params_')]
            s = set()
            for p in pars:
                s.update(p)
            self.__params_cache = tuple(s)
        return self.__params_cache
    params = property(_getParams)
    
    def getParamWithError(self,parname):
        """
        this retreives the requested parameter and returns it as a tuple of the
        form
        
        (param,lowererr,uppererr)
        """
        from operator import isSequenceType
        
        if parname not in self.params:
            raise ValueError('invalid parameter name %s'%parname)
        
        val = getattr(self,parname)
        if hasattr(self,parname+'_err'):
            err = getattr(self,parname+'_err')
        else:
            err = (0,0)
        
        if isSequenceType(err):
            return (val,err[0],err[1])
        else:
            return (val,err,err)
    
    def _exportParams(self):
        pd=dict([(p,getattr(self,p)) for p in self.params])
        globals().update(pd)
        
    def _removeParams(self):
        from warnings import warn
        d=globals()
        for p in self.params:
            out=d.pop(p,None)
            if out is None:
                warn('Cosmological parameter %s not present despite being current cosmology'%p)
        
    def __setattr__(self,name,value):
        object.__setattr__(self, name, value)
        if self._autoupdate:
            self._exportParams()
    
class FRWCosmology(Cosmology):
    """
    A cosmology based on the FRW metric  with a global density, a matter 
    density, and a radiation density, and a comological constant as 
    specified at z=0
    
    default values are approximately LambdaCDM
    """
    _params_=('omega','omegaR','omegaM','omegaL')
    
    H0=72
    omega = property(lambda self:self.omegaR+self.omegaM+self.omegaL)
    omegaR=0 #radiation density
    omegaM=0.3 #matter density
    omegaL=0.7 #dark energy density
    
    @property
    def omegaK(self):
        return 1-self.omegaR-self.omegaM-self.omegaL
    
    def H(self,z):
        z=np.array(z)
        M,L,R=self.omegaM,self.omegaL,self.omegaR
        K=1-M-L-R
        a=1/(1+z)
        return self.H0*(R*a**-4 + M*a**-3 + L + K*a**-2)**0.5
        
    def computeOmegaRz(self,z):
        """
        compute radiation density at arbitrary redshift
        """
        a=1/(1+z)
        return (self.rhoC(0)/self.rhoC(z))*self.omegaR*a**-4
        
    def computeOmegaMz(self,z):
        """
        compute matter density at arbitrary redshift
        """
        a=1/(1+z)
        return (self.rhoC(0)/self.rhoC(z))*self.omegaM*a**-3
    
    def computeOmegaLz(self,z):
        """
        compute cosmological constant density at arbitrary redshift
        """
        return self.omegaL*(self.rhoC(0)/self.rhoC(z))
    
    def computeOmegaKz(self,z,units='cgs'):
        """
        compute curvature density at arbitrary redshift
        """
        a=1/(1+z)
        return (self.rhoC(0)/self.rhoC(z))*self.omegaK*a**-2
    
    def rhoC(self,z=0,units='cgs'):
        """
        critical density at a given redshift
        
        units can be 'cgs' or 'cosmological' (Mpc,Msun)
        TODO:check
        """
        H = self.H(z)*1e5*1e-6*pcpercm #km/s->cm/s, Mpc->cm
        cgsres = 3*H*H/(8*pi*G)
        
        if units == 'cgs':
            return cgsres
        elif units == 'cosmological':
            return cgsres/Ms*(1e-6*pcpercm)**-3#g->Msun,cm^-3->Mpc^-3
        else:
            raise ValueError('unrecognized units')
        
    
    def rho(self,z=0,units='cgs'):
        """
        mean density in this cosmology
        
        units can be 'cgs' or 'cosmological' (Mpc,Msun)
        """
        return self.omega*self.rhoC(z,units)
    
    def deltavir(self,z=0):
        """
        virial overdensity as paramaterized in Bryan&Norman 98 for omega=1
        """
        if self.omegaK !=0:
            raise NotImplementedError("can't compute for omega!=1")
        om = self.computeOmegaMz(z)
        x= om - 1
        return (18*pi**2+82*x-39*x**2)/om

class WMAP7Cosmology(FRWCosmology):
    """
    WMAP7-only (http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm)
    """
    _params_ = ('t0','sigma8','omegaB','omegaC','ns')
    t0 = 13.71#Gyr
    t0_err = .13
    sigma8 = .801
    sigma8_err = .030
    ns = .963
    ns_err = .014
    H0 = 71.0
    H0_err = 2.5
    omegaB = 0.044
    omegaB_err = .0028
    omegaC = 0.222
    omegaC_err = 0.026
    omegaL = 0.734
    omegaL_err = .029
    omegaM = property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)
    
class WMAP7BAOH0Cosmology(FRWCosmology):
    """
    WMAP7+BAO+H0 (http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7_bao_h0.cfm)
    """
    _params_ = ('t0','sigma8','omegaB','omegaC','ns')
    t0 = 13.78#Gyr
    t0_err = .11
    sigma8 = 0.809
    sigma8_err = .024 
    ns = .963
    ns_err = .012
    H0 = 70.4
    H0_err = (1.4,1.3)
    omegaB = 0.045
    omegaB_err = .0016
    omegaC = 0.227
    omegaC_err = .014
    omegaL = 0.728
    omegaL_err = (.016,.015)
    omegaM = property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)
    
class WMAP5Cosmology(FRWCosmology):
    """
    WMAP5-only (http://lambda.gsfc.nasa.gov/product/map/dr3/parameters_summary.cfm)
    """
    _params_=('t0','sigma8','omegaB','omegaC','ns')
    t0 = 13.69 #Gyr
    t0_err = .13
    sigma8 = .796
    sigma8_err = .036
    ns = 0.963
    ns_err = (.015,.014)
    H0 = 71.9
    H0_err = (2.7,2.6)
    omegaB = 0.044
    omegaB_err = .003
    omegaC = 0.214
    omegaC_err = .027
    omegaL = 0.742
    omegaL_err = .030
    omegaM = property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)
    
class WMAP5BAOSNCosmology(FRWCosmology):
    """
    WMAP5+BAO+SN (http://lambda.gsfc.nasa.gov/product/map/dr3/parameters_summary.cfm)
    """
    _params_=('t0','sigma8','omegaB','omegaC','ns')
    t0 = 13.73 #Gyr
    t0_err = .12
    sigma8 = .817
    sigma8_err = .026
    ns = 0.960
    ns_err = (.013,.014)
    H0 = 70.1
    H0_err = 1.3
    omegaB = 0.046
    omegaB_err = .0015
    omegaC = 0.233
    omegaC_err = .013
    omegaL = 0.721
    omegaL_err = .015
    omegaM=property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)
    
class WMAP3Cosmology(FRWCosmology):
    """
    WMAP3 only (http://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_wmap.cfm)
    """
    #_params_=('t0','sigma8')
    _params_=('sigma8','omegaB','omegaC','ns')
    #t0=13.69 #Gyr
    sigma8 = .761
    sigma8_err = (.048,.049)
    ns = 0.958
    ns_err = .016
    H0 = 73.2
    H0_err = (3.2,3.1)
    omegaB = 0.044
    omegaB_err = .0014
    omegaC = 0.224
    omegaC_err = .014
    omegaL = 0.732
    omegaL_err = .034
    omegaM = property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)
    
class WMAP3AllCosmology(FRWCosmology):
    """
    WMAP3+all (http://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_all.cfm)
    """
    #_params_=('t0','sigma8')
    _params_=('sigma8','omegaB','omegaC','ns')
    #t0=13.69 #Gyr
    sigma8 = .776
    sigma8_err = (.032,.031)
    ns = 0.947
    ns_err = .015
    H0 = 70.4
    H0_err = (1.6,1.5)
    omegaB = 0.042+.002/3 #omega=1
    omegaB_err = .0013 
    omegaC = 0.197+.002/3 #omega=1
    omegaC_err = (.0079,.0077)
    omegaL = 0.759+.002/3 #omega=1
    omegaL_err = .018
    omegaM = property(lambda self:self.omegaB+self.omegaC)
    omegaM_err = property(lambda self:self.omegaB_err+self.omegaC_err)


__current_cosmology=WMAP7BAOH0Cosmology() #default value
__current_cosmology._exportParams()
__cosmo_registry={}

def register_cosmology(cosmocls,name=None):
    """
    Add the provided subclass of Cosmology to the cosmology registry
    
    if name is None, the name will be inferred from the class name, otherwise
    
    """
    if not name:
        name = cosmocls.__name__
    name = name.lower().replace('cosmology','')
    try:
        if not issubclass(cosmocls,Cosmology):
            raise TypeError("Supplied object is not a subclass of Cosmology")
    except TypeError:
        raise TypeError("Supplied object to register is not a class")
    
    __cosmo_registry[name]=cosmocls
    
#register all Cosmologies in this module
for o in locals().values():
    if type(o)==type and issubclass(o,Cosmology) and o != Cosmology:
        register_cosmology(o)

def choose_cosmology(nameorobj,autoupdate=True,args=None,kwargs=None):
    """
    Change the currently active cosmology and export its cosmological parameters
    into the package namespace.
    
    :param nameorobj: the new cosmology to use
    :type nameorobj: string or :class:`Cosmology` object
    :param autoupdate: 
        If True, the cosmology object will automatically propogate changes to
        its parameters up to the module variables. Otherwise,
        :func:`update_cosmology` must be called explicitly to have this behavior
        occur.
    :type autoupdate: bool
    :param args: 
        If `nameorobj` is a string, these are passed in as positional arguments
        to the object initializer (if not None). Otherwise it is ignored.
    :type args: sequence or None
    :param args: 
        If `nameorobj` is a string, these are passed in as keyword arguments to
        the object initializer (if not None). Otherwise it is ignored.
    :type args: dictionary or None
    
    :returns: the :class:`Cosmology` object after being assigned as current.
    """
    global __current_cosmology
    
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}
    
    if isinstance(nameorobj,basestring):
        c = __cosmo_registry[nameorobj.lower()](*args,**kwargs)
    elif isinstance(nameorobj,Cosmology):
        c = nameorobj
        if c.__class__ not in __cosmo_registry.values():
            register_cosmology(c.__class__)
    
    __current_cosmology._removeParams()
    __current_cosmology._autoupdate = False
    try:
        c._exportParams()
        __current_cosmology =  c
    except:
        __current_cosmology._exportParams()
        
    c._autoupdate = bool(autoupdate)
        
    return c
    
    
def get_cosmology(name=None):
    """
    If name is None, will retreive the currently in use Cosmology instance.
    Otherwise, returns the subclass of :class:`Cosmology` with the provided
    name.
    """
    if name is None:
        return __current_cosmology
    else:
        return __cosmo_registry[name]
    
def update_cosmology():
    """
    updates the package-level variables for changes in the current Cosmology 
    object
    """
    __current_cosmology._exportParams()

def get_registry_names():
    """
    Returns the names of all cosmology types in the registry
    """
    return __cosmo_registry.keys()
