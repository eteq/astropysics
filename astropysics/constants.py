from __future__ import division
from math import pi
import numpy as np

#<--------------------------------Constants------------------------------------>
#all in cgs including gaussian esu for charge
prop=property(lambda:5)
unit_system='cgs'

G=6.673e-8
mp=1.67262171e-24
me=9.1093897e-28
e=4.8032068e-10
Ms=1.9891e33
Mj=1.8986e30
Me=5.9742e27
Rs=6.96e10
Rj=7.1492e9
Re=6.371e8
De=1.49597887e13
kb=1.3807e-16 
#sigma=5.6704e-5
c=2.99792458e10
h=6.626068E-27


#<-------------------------------Conversions----------------------------------->
ergperev=1.60217646e-12
secperyr=3.15576926e7
secpergyr=secperyr*1e9
pcpercm=1/3.08568025e18
pcperly=1/3.26

#TODO:rethink flux units and adapt to BlackbodyModel
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
    
    while not required
    
    Cosmologies should also include a sequence called "_params_" with a list of
    strings that specify the names of the cosmological parameters associated to
    be exported to the constants module
    """
    #TODO:use ABCs in py 2.6
    _params_=('H0','omega')
    
    H0=0
    omega=0
    
    def __init__(self,*args,**kwargs):
        ps=self._params_
    h = property(lambda self:self.H0/100.0)
    h7 = h70 = property(lambda self:self.H0/70.0)
    __params_cache=None
    def _getParams(self):
        if self.__params_cache is None:
            import inspect
            pars= [cls._params_ for cls in inspect.getmro(self.__class__) if 
                    hasattr(cls,'_params_')]
            s=set()
            for p in pars:
                s.update(p)
            self.__params_cache = tuple(s)
        return self.__params_cache
    params=property(_getParams)
    
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
        
    def rhoC(self):
        """
        critical density
        """
        return 3*self.H0*self.H0/(8*pi*G)
    
    def rho(self):
        """
        mean density in this cosmology
        """
        return self.omega*self.rhoC()
    
class MLRCosmology(Cosmology):
    """
    A cosmology with a lambda, a matter density, and a radiation density
    
    default values are approximately LambdaCDM
    """
    _params_=('omegaR','omegaM','omegaL')
    
    H0=72
    omega = property(lambda self:self.omegaR+self.omegaM+self.omegaL)
    omegaR=0 #radiation density
    omegaM=0.3 #matter density
    omegaL=0.7 #dark energy density
    
    def H(self,z):
        z=np.array(z)
        M,L,R=self.omegaM,self.omegaL,self.omegaR
        K=1-M-L-R
        a=1/(1+z)
        H=self.H0*(R*a**-4 + M*a**-3 + L + K*a**-2)**0.5
        
    
class WMAP5Cosmology(MLRCosmology):
    _params_=('t0','sigma8')
    t0=13.69 #Gyr
    sigma8=.796
    H0=70.1
    omegaB=0.046
    omegaC=0.233
    omegaL=0.721
    omegaM=property(lambda self:self.omegaB+self.omegaC)

class WMAP3Cosmology(MLRCosmology):
    _params_=('t0','sigma8')
    t0=13.69 #Gyr
    sigma8=.776
    H0=70.4
    omegaB=0.044
    omegaC=0.224
    omegaL=0.732
    omegaM=property(lambda self:self.omegaB+self.omegaC)


__current_cosmology=WMAP5Cosmology() #default value
__current_cosmology._exportParams()
__cosmo_registry={}

def register_cosmology(cosmocls,name=None):
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

def choose_cosmology(name,*args,**kwargs):
    """
    generate a new cosmology with the args and kwargs going into the initializer
    """
    global __current_cosmology
    c = __cosmo_registry[name.lower()](*args,**kwargs)
    __current_cosmology._removeParams()
    c._exportParams()
    
    __current_cosmology =  c
    
def get_cosmology(name=None):
    """
    if name is None, will retreive the currently in use Cosmology instance.
    Otherwise, returns the Class object for the requested cosmology
    """
    if name is None:
        return __current_cosmology
    else:
        return __cosmo_registry[name]
        
