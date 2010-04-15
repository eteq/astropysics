#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains builtin function models that use the framework of the
:mod:`core` module.

Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.models.builtins
   :parts: 1
   
Module API
^^^^^^^^^^

"""

from __future__ import division,with_statement

from core import *
from ..spec import HasSpecUnits as _HasSpecUnits
from math import e


class ConstantModel(FunctionModel1DAuto):
    """
    the simplest model imaginable - just a constant value
    """
    def f(self,x,C=0):
        return C*np.ones_like(x)
    
    def derivative(self,x,dx=1):
        return np.zeros_like(x)
    
    def integrate(self,lower,upper,**kwargs):
        if 'jac' in kwargs and kwargs['jac'] is not None:
            return FunctionModel1D.integrate(self,lower,upper,**kwargs)
        else:
            return self.C*(upper-lower)
    
class LinearModel(FunctionModel1DAuto):
    """
    y=mx+b linear fit
    
    A number of different fitting options are available - the `fittype`
    attribute determines which type should be used on calls to
    fitData.  It can be:
    
    * 'basic': analytically compute the parameters from simple least-squares
       form.  If `weights` are provided they will be interpreted as
       1/yerr or sqrt((x/xerr)^2+(y/yerr)^2) if a 2-tuple
    * 'yerr': same as weighted version of basic, but `weights`
      are interpreted as y-error instead of inverse
    * 'fiterrxy': allows errors in both x and y using the chi^2
      from the fitexy algorithm from numerical recipes.  `weights`
      must be a 2-tuple (xerr,yerr)
    """
    
    def f(self,x,m=1,b=0):
        return m*x+b
    
    def _linearFit(self,x,y,fixedpars=(),weights=None,**kwargs):
        """
        does least-squares fit on the x,y data
        
        fixint and fixslope can be used to specify the intercept or slope of the 
        fit or leave them free by having fixint or fixslope be False or None
        
        lastfit stores ((m,b),dm,db)
        """  
        if self.fittype == 'basic':
            if weights is None:
                (m,b),merr,berr,dy = self.fitBasic(x,y,'m' in fixedpars,'b' in fixedpars)
            else:
                if len(fixedpars)>0:
                    raise ValueError('cannot fix parameters and do weighted fit')
                weights = np.array(weights,copy=False)
                if len(weights.shape) == 2:
                    weights = ((x/weights[0])**2+(y/weights[1])**2)**0.5
                m,b,merr,berr = self.fitWeighted(x,y,1/weights,**kwargs)
        elif self.fittype == 'yerr':
            weights = np.array(weights,copy=False) if weights is not None else None
            if weights is not None and len(weights.shape) == 2:
                weights = weights[1]
            m,b,merr,berr = self.fitWeighted(x,y,weights,**kwargs)
        elif self.fittype == 'fiterrxy':
            if len(fixedpars)>0:
                raise ValueError('cannot fix parameters and do fiterrxy fit')
            if weights is None:
                xerr = yerr = None
            else:
                weights = np.array(weights,copy=False)
                if len(weights.shape)==1:
                    xerr = yerr = 1/(weights*2**-0.5)
                else:
                    xerr,yerr = 1/weights
            m,b = self.fitErrxy(x,y,xerr,yerr,**kwargs)
            merr = berr = None
        else:
            raise ValueError('invalid fittype %s'%self.fittype)
            
        return ((m,b),merr,berr)
    
    _fittypes = {'basic':_linearFit,'yerr':_linearFit,'fiterrxy':_linearFit}
    fittype = 'basic'
    
    def derivative(self,x,dx=1):
        return np.ones_like(x)*self.m
    
    def integrate(self,lower,upper):
        m,b = self.m,self.b
        return m*upper*upper/2+b*upper-m*lower*lower/2+b*lower
    
    def fitBasic(self,x,y,fixslope=False,fixint=False):
        """
        Does the traditional fit based on simple least-squares
        """
        N=len(x) 
        if not fixint and not fixslope:
            if len(y)!=N:
                raise ValueError('data arrays are not same length!')
            sxsq=np.sum(x*x)
            sx,sy=np.sum(x),np.sum(y)
            sxy=np.sum(x*y)
            delta=N*sxsq-sx**2
            m=(N*sxy-sx*sy)/delta
            b=(sxsq*sy-sx*sxy)/delta
            dy=(y-m*x-b).std(ddof=2)
            dm=dy*(sxsq/delta)**0.5 
            db=dy*(N/delta)**0.5 
            
        elif not fixint:
            
            m,dm=self.m,0
            
            b=np.sum(y-m*x)/N 
            
            dy=(y-m*x-b).std(ddof=1)
            #db= sum(dy*dy)**0.5/N
            db = dy
            
        elif not fixslope:
            b,db=self.b,0
            
            sx=np.sum(x)
            sxy=np.sum(x*y)
            sxsq=np.sum(x*x)
            m=(sxy-b*sx)/sxsq
            
            dy=(y-m*x-b).std(ddof=1) 
            #dm=(np.sum(x*dy*x*dy))**0.5/sxsq
            dm = dy*sxsq**-0.5
        else:
            raise ValueError("can't fix both slope and intercept")
        
        return np.array((m,b)),dm,db,dy
#        if weights is not None:
#            if fixedpars or len(kwargs)>0:
#                from warnings import warn
#                warn("customized exact linear fit not yet available for fixed pars")
#                kwargs=kwargs.copy()
#                kwargs['x']=x
#                kwargs['y']=y
#                kwargs['method']='leastsq'
#                kwargs['fixedpars']=fixedpars
#                kwargs['weights']=weights
#                return FunctionModel1D.fitData(self,**kwargs)
#            m,b,dm,db = self.weightedFit(x,y,1/weights,False)
#            dy = (y-m*x-b).std(ddof=1)
#            return (np.array((m,b)),dm,db,dy)
        
        
    
    def fitWeighted(self,x,y,sigmay=None,doplot=False):
        """
        does a linear weighted least squares fit and computes the coefficients 
        and errors
        
        fit is y=B*x+A
        
        if sigma is None, the weights are all equal - otherwise, it's the stddev 
        of the y values
        
        returns B,A,sigmaB,sigmaA
        """
#        raise NotImplementedError('needs to be adapted to astro.models')
        from numpy import array,ones,sum
        if sigmay is None:
            sigmay=ones(len(x))
        else:
            sigmay=array(sigmay)
        if len(x)!=len(y)!=len(sigmay):
            raise ValueError('arrays not matching lengths')
        N=len(x)
        x,y=array(x),array(y)
        
        w=1.0/sigmay/sigmay
        delta=sum(w)*sum(w*x*x)-(sum(w*x))**2
        A=(sum(w*x*x)*sum(w*y)-sum(w*x)*sum(w*x*y))/delta
        B=(sum(w)*sum(w*x*y)-sum(w*x)*sum(w*y))/delta
        diff=y-A-B*x
        sigmaysq=sum(w*diff*diff)/(sum(w)*(N-2)/N) #TODO:check
        sigmaA=(sigmaysq*sum(w*x*x)/delta)**0.5
        sigmaB=(sigmaysq*sum(w)/delta)**0.5
        
        if doplot:
            from matplotlib.pyplot import plot,errorbar,legend
            errorbar(x,y,sigmay,fmt='o',label='Data')
            plot(x,B*x+A,label='Best Fit')
            plot(x,(B+sigmaB)*x+A-sigmaA,label='$1\sigma$ Up')
            plot(x,(B-sigmaB)*x+A+sigmaA,label='$1\sigma$ Down')
            legend(loc=0)
        
        return B,A,sigmaB,sigmaA
    
    def fitErrxy(self,x,y,xerr,yerr,**kwargs):
        """
        Uses the chi^2 statistic (ydata-ymodel)^2/(yerr^2+m^2*xerr^2)
        to fit the data for errors in both x and y
        
        kwargs are passed into scipy.optimize.leastsq
        
        fitting results are saved to `lastfit`
        """
        from scipy.optimize import leastsq
        if xerr is None and yerr is None:
            def chi(v,x,y,xerr,yerr):
                return y-self.f(x,*v)
        elif xerr is None:
            def chi(v,x,y,xerr,yerr):
                wsqrt = xerr
                return (y-self.f(x,*v))/wsqrt
        elif yerr is None:
            def chi(v,x,y,xerr,yerr):
                wsqrt = yerr
                return (y-self.f(x,*v))/wsqrt
        else:
            def chi(v,x,y,xerr,yerr):
                wsqrt = (yerr**2+(v[0]*xerr)**2)**0.5
                return (y-self.f(x,*v))/wsqrt
        
        kwargs.setdefault('full_output',1)
        self.lastfit = leastsq(chi,(self.m,self.b),args=(x,y,xerr,yerr),**kwargs)
        self.data = (x,y,(xerr,yerr))
        
        self._fitchi2 = np.sum(chi(self.lastfit[0],x,y,xerr,yerr)**2)
        
        return self.lastfit[0]
    
    def pointSlope(self,m,x0,y0):
        """
        sets model parameters for the given slope that passes through the point
        """
        self.m = m
        self.b = y0-m*x0
        
    def twoPoint(self,x0,y0,x1,y1):
        """
        sets model parameters to pass through two lines (identical behavior
        in fitData)
        """
        self.pointSlope((y0-y1)/(x0-x1),x0,y0)
        
    @staticmethod
    def fromPowerLaw(plmod,base=10):
        """
        Takes a PowerLawModel and converts it to a linear model assuming 
        ylinear = log_base(ypowerlaw) and xlinear = log_base(xpowerlaw) 
        
        returns the new LinearModel instance
        """
        if base == 'e' or base == 'ln':
            base = np.exp(1)
        
        logfactor=1/np.log(base)
        
        m = plmod.p
        b = logfactor*np.log(plmod.A)
        
        if plmod.data is not None:
            data = list(plmod.data)
            data[0] = logfactor*np.log(data[0])
            data[1] = logfactor*np.log(data[1])
            
        lmod = LinearModel(m=m,b=b)
        lmod.data = tuple(data)
        return lmod
        
    
class QuadraticModel(FunctionModel1DAuto):
    """
    2-degree polynomial
    """
    def f(self,x,c2=1,c1=0,c0=0):
        return c2*x*x+c1*x+c0

class PolynomialModel(FunctionModel1DAuto):
    """
    arbitrary-degree polynomial
    """
    
    paramnames = 'c'
    
    #TODO: use polynomial objects that are only updated when necessary
    def f(self,x,*args): 
        return np.polyval(np.array(args)[::-1],x)
    
    def derivative(self,x):
        return np.polyder(np.array(self.parvals)[::-1])(x)

    def integrate(self,lower,upper):
        p = np.polyint(np.array(self.parvals)[::-1])
        return p(upper)-p(lower)
    
class FourierModel(FunctionModel1DAuto):
    paramnames = ('A','B')
    #note that B0 has no effect
    
    def f(self,x,*args):
        xr = x.ravel()
        n = len(args)/2
        As = np.array(args[::2]).reshape((n,1))
        Bs = np.array(args[1::2]).reshape((n,1))
        ns = np.arange(len(As)).reshape((n,1))
        res = np.sum(As*np.sin(ns*xr),axis=0)+np.sum(Bs*np.cos(ns*xr),axis=0)
        return res.reshape(x.shape)
#        val = np.empty_like(x)
#        for n,(A,B) in enumerate(zip(args[::2],args[:1:2])):
#            val += A*sin(n*x)+B*cos(n*x)
#        return val

class GaussianModel(FunctionModel1DAuto):
    """
    Normalized 1D gaussian function
    """
    def f(self,x,A=1,sig=1,mu=0):
        tsq=(x-mu)*2**-0.5/sig
        return A*np.exp(-tsq*tsq)*(2*pi)**-0.5/sig
    
    def _getPeak(self):
        return self(self.mu)
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
        
    peak=property(_getPeak,_setPeak)
    
    __fwhmfactor = 2*(2*np.log(2))**0.5
    def _getFWHM(self):
        return self.sig*self.__fwhmfactor
    def _setFWHM(self,val):
        self.sig = val/self.__fwhmfactor
    FWHM = property(_getFWHM,_setFWHM,doc='Full Width at Half Maximum')
    
        
    def derivative(self,x,dx=1):
        sig=self.sig
        return self(x)*-x/sig/sig
    
    def integrate(self,lower,upper,**kwargs):
        if len(kwargs)>0:
            return super(GaussianModel,self).integrate(lower,upper,**kwargs)
        else:
            from scipy.special import erf
            l = (lower-self.mu)/self.sig
            u = (upper-self.mu)/self.sig
            uval = erf(u*2**-0.5)
            lval = erf(-l*2**-0.5)
            return self.A*(uval+lval)/2
        
    @property
    def rangehint(self):
        return (self.mu-self.sig*4,self.mu+self.sig*4)
    
class DoubleGaussianModel(FunctionModel1DAuto):
    """
    Two Normalized 1D gaussian functions, fixed to be of opposite sign
    A is the positive gaussian, while B is negative
    note that fitting often requires the initial condition to have the 
    upper and lower approximately correct
    """
    def f(self,x,A=1,B=1,sig1=1,sig2=1,mu1=-0.5,mu2=0.5):
        A,B=abs(A),-abs(B) #TOdO:see if we should also force self.A and self.B
        tsq1=(x-mu1)*2**-0.5/sig1
        tsq2=(x-mu2)*2**-0.5/sig2
        return (A*np.exp(-tsq1*tsq1)/sig1+B*np.exp(-tsq2*tsq2)/sig2)*(2*pi)**-0.5
    
    @property
    def rangehint(self):
        lower1 = self.mu1-self.sig1*4
        lower2 = self.mu2-self.sig2*4
        upper1 = self.mu1+self.sig1*4
        upper2 = self.mu2+self.sig2*4
        return(min(lower1,lower2),max(upper1,upper2))
    
    @staticmethod
    def autoDualModel(x,y,taller='A',wider='B',**kwargs):
        """
        generates and fits a double-gaussian model where one of the peaks is
        on top of the other and much stronger.
        the taller and wider argument must be either 'A' or 'B' for the positive
        and negative components, respectively
        kwargs go into the fitData calls
        """
        gm=GaussianModel()
        gm.fitData(x,y,**kwargs)
        dgm=DoubleGaussianModel()
        dgm.mu1=dgm.mu2=gm.mu
        if taller == 'A':
            dgm.A=gm.A
            dgm.B=gm.A/2
            dgm.sig1=gm.sig
            if wider =='A':
                dgm.sig2=gm.sig/2
            elif wider =='B':
                dgm.sig2=gm.sig*2
            else:
                raise ValueError('unrecognized wider component')
            
            dgm.fitData(x,y,fixedpars=('mu1','A','sig1'),**kwargs)
        elif taller == 'B':
            dgm.B=gm.A
            dgm.A=gm.A/2
            dgm.sig2=gm.sig
            if wider =='B':
                dgm.sig1=gm.sig/2
            elif wider =='A':
                dgm.sig1=gm.sig*2
            else:
                raise ValueError('unrecognized wider component')
            
            dgm.fitData(x,y,fixedpars=('mu2','B','sig2'),**kwargs)
        else:
            raise ValueError('unrecognized main component')
        
        dgm.fitData(x,y,fixedpars=(),**kwargs)
        return dgm
    
class LognormalModel(FunctionModel1DAuto):
    """
    A Lognormal model, e.g. f = A e^(-(log10(x)-mulog)^2/2 siglog^2)
    """
    def f(self,x,A=1,mulog=0,siglog=1):
        return A*np.exp(-0.5*((np.log(x)-mulog)/siglog)**2)
    
    @property
    def rangehint(self):
        return (np.exp(self.mulog-self.siglog*4),np.exp(self.mulog+self.siglog*4))
    
class LorentzianModel(FunctionModel1DAuto):
    def f(self,x,A=1,gamma=1,mu=0):
        return A*gamma/pi/(x*x-2*x*mu+mu*mu+gamma*gamma)
    
    def _getPeak(self):
        return self(self.mu)
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
        
    peak=property(_getPeak,_setPeak)
    
    @property
    def rangehint(self):
        return(self.mu-self.gamma*6,self.mu+self.gamma*6)
    
class VoigtModel(GaussianModel,LorentzianModel):
    """
    Convolution of a Gaussian and Lorentzian profile
    """
    def f(self,x,A=1,sig=0.5,gamma=0.5,mu=0):
        from scipy.special import wofz
        if sig == 0:
            return LorentzianModel.f(self,x,A,sig,mu)
        else:
            w=wofz(((x-mu)+1j*gamma)*2**-0.5/sig)
            return A*w.real*(2*pi)**-0.5/sig
    
    @property
    def rangehint(self):
        halfwidth = 3*(self.gamma+self.sig)
        return(self.mu-halfwidth,self.mu+halfwidth)
        
class MoffatModel(FunctionModel1DAuto):
    """
    Moffat function given by:
    A*(beta-1)/(pi alpha^2) [1+(r/alpha)^2]^-beta
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
    
class ExponentialModel(FunctionModel1DAuto):
    """
    exponential function Ae^(kx)
    """
    def f(self,x,A=1,k=1):
        return A*np.exp(k*x)
    
    @property
    def rangehint(self):
        ak = np.abs(self.k)
        return(-1.5/ak,1.5/ak)
    
class PowerLawModel(FunctionModel1DAuto):
    """
    A single power law model Ax^p 
    """
    def f(self,x,A=1,p=1):
        return A*x**p
    
    _fittypes=['linearized']
    fittype = 'leastsq'
    
    def fitLinearized(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        lm = LinearModel(m=self.p,b=np.log10(self.A))
        lm.fittype = 'basic'
        logx = np.log10(x)
        logy = np.log10(y)
        
        fixedps = []
        if 'A' in fixedpars:
            fixedps.append('b')
        if 'p' in fixedpars:
            fixedps.append('m')
        
        lm.fitData(logx,logy,tuple(fixedps),**kwargs)
        
        return np.array([10**lm.b,lm.m])
    
    @staticmethod
    def fromLinear(lmod,base=10):
        """
        Takes a LinearModel and converts it to a power law model assuming 
        ylinear = log_base(ypowerlaw) and xlinear = log_base(xpowerlaw) 
        
        returns the new PowerLawModel instance
        """
        if base == 'e' or base == 'ln':
            base = np.exp(1)
        
        p = lmod.m
        A = base**lmod.b
        
        if lmod.data is not None:
            data = list(lmod.data)
            data[0] = base**data[0]
            data[1] = base**data[1]
        else:
            data = None
            
        plmod = PowerLawModel(p=p,A=A)
        if data is not None:
            plmod.data = tuple(data)
        return plmod
    
class SinModel(FunctionModel1DAuto):
    """
    A trigonometric model A*sin(k*x+p)
    """
    def f(self,x,A=1,k=2*pi,p=0):
        return A*np.sin(k*x+p)
    
    def derivative(self,x,dx=1):
        A,k,p=self.A,self.k,self.p
        return A*k*np.cos(k*x+p)
    
    def integrate(self,lower,upper):
        A,k,p=self.A,self.k,self.p
        return A*(np.cos(k*lower+p)-np.cos(k*upper+p))/k
    
class TwoPowerModel(FunctionModel1DAuto):
    """
    A model that smoothly transitions between two power laws at the turnover 
    point xs.  a is the inner slope, b is the outer slope
    A and fxs are the same parameter - A is the absolute normalization, and fxs
    is the function's value at xs
    """
    def f(self,x,A=1,xs=1,a=1,b=2):
        return A*((x+xs)**(b-a))*(x**a)
    
    def _getFxs(self):
        A,xs,a,b=self.A,self.xs,self.a,self.b
        return A*xs**b*2**(b-a)
    
    def _setFxs(self,fxs):
        xs,a,b=self.xs,self.a,self.b
        self.A=fxs*xs**-b*2**(a-b)
    
    fxs=property(fget=_getFxs,fset=_setFxs)
    
    
class TwoSlopeModel(FunctionModel1DAuto):
    """
    This model smoothly transitions from linear with one slope to linear with
    a different slope. It is the linearized equivalent of TwoPowerModel.
    
    specifically, a*(x-xs)+(b-a)*log_base(1+base^(x-xs))+C
    """
    
    fixedpars = ('base',)
    
    def f(self,x,a=1,b=2,C=0,xs=0,base=e):
        z = x-xs
        return a*z+(b-a)*np.log(1+base**z)/np.log(base)+C
    
    @property
    def rangehint(self):
        return(self.xs-3,self.xs+3)
    
class TwoSlopeDiscontinuousModel(FunctionModel1DAuto):
    """
    This model sharply transitions between two slopes
    
    a is the slope for small x and b for large x, with xs the transition
    point.  The intercept C is for the a-slope
    """
    def f(self,x,a=1,b=2,C=0,xs=0):
        xl = x.copy()
        xl[x>xs]  = 0
        xu = x.copy()
        xu[x<=xs]  = 0
        return a*xl+b*xu+(a*xs-b*xs)*(xu!=0)+C
    
    @property
    def rangehint(self):
        return(self.xs-3,self.xs+3)
        
class TwoSlopeATanModel(FunctionModel1DAuto):
    """
    This model transitions between two asymptotic slopes with
    an additional parameter that allows for a variable transition region 
    size.  The functional form is:
    
    y = ax arctan(-(x-x0)/w) + bx arctan((x-x0)/w) + c
    
    a is the slope for small x, b for large x, c is the value at x=x0, 
    x0 is the location of the transition, and w is the width of the transition
    
    """
    def f(self,x,a=1,b=2,c=0,x0=0,w=1):
        xoff = x-x0
        tana = np.arctan(-xoff/w)/pi+0.5
        tanb = np.arctan(xoff/w)/pi+0.5
        return a*xoff*tana+b*xoff*tanb+c
    
    def criticalPoint(self):
        """
        computes the critical point where this model inflects or has a max/min
        
        if a = b, returns 0 (undefined critical point)
        """
        if self.a == self.b:
            return 0
        else:
            from scipy.optimize import newton
            bma = self.b - self.a
            bpa = self.a + self.b
            toroot = lambda u: bma*(np.arctan(u)+u/(1+u**2))+bpa*pi/2
            ur = newton(toroot,0)
            return ur*self.w+self.x0
    
    @property
    def rangehint(self):
        return self.x0-3*self.w,self.x0+3*self.w
    
class BlackbodyModel(FunctionModel1DAuto,_HasSpecUnits):
    """
    a Planck blackbody radiation model.  

    y-axis is assumed to be specific intensity
    """
    from ..constants import h,c,kb
    
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
        sets A so that the output is specific intensity/surface brightness
        """
        self.A = 1
    
    def setFlux(self,radius,distance):
        """
        sets A so that the output is the flux at the specified distance from
        a spherical blackbody with the specified radius
        """
        from .phot import intensity_to_flux
        self.A = intensity_to_flux(radius,distance)
        
    def getFlux(self,x,radius=None,distance=None):
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
        determines the radius of a spherical blackbody at the specified distance
        assuming the flux is given by the model at the given temperature
        """
        return (self.A*distance*distance/pi)**0.5
     
    def getFluxDistance(self,radius):
        """
        determines the distance to a spherical blackbody of the specified radius
        assuming the flux is given by the model at the given temperature
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
        uses the Wien Displacement Law to calculate the temperature given a
        peak wavelength or peak temperature
        
        if input is None, the peak location corresponding to the current 
        temperature is returned
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


    
class _InterpolatedModel(DatacentricModel1DAuto):
    
    _fittypes=['interp']
    fittype = 'interp'
    
    def __init__(self):
        super(_InterpolatedModel,self).__init__()
        self.i1d = lambda x:x #default should never be externally seen
    
    def f(self,x):
        if self.data is not None:
            res = self.i1d(x)
            xd,yd = self.data[0],self.data[1]
            mi,mx = np.min(xd),np.max(xd)
            res[x<mi] = yd[mi==xd][0]
            res[x>mx] = yd[mx==xd][0]
            return res
        else:
            return x
    
    def fitData(self,x,y,**kwargs):
        kwargs['savedata'] = True
        return super(_InterpolatedModel,self).fitData(x,y,**kwargs)
    
    def fitInterp(self,x,y,fixedpars=(),**kwargs):
        from scipy.interpolate import interp1d
        xi = np.argsort(x)
        self.i1d = interp1d(x[xi],y[xi],kind=self.kind,bounds_error=False)
        
        return []
        
class LinearInterpolatedModel(_InterpolatedModel):
    """
    A model that is the linear interpolation of the data, or if out of bounds, 
    fixed to the edge value.
    """
    kind = 'linear'
    
class NearestInterpolatedModel(_InterpolatedModel):
    """
    A model that is the interpolation of the data by taking the value of the 
    nearest point
    """
    kind = 'nearest'
   
class SmoothSplineModel(DatacentricModel1DAuto):
    """
    this uses a B-spline as a model for the function.  Note that by
    default the parameters are not tuned - the input smoothing and 
    degree are left alone when fitting
    
    the scipy.interpolate.UnivariateSpline class is used to
    do the calculation (in the "spline" attribute) 
    """
    def __init__(self):
        super(SmoothSplineModel,self).__init__()
        
        self._oldd = self._olds = self._ws = self._inits = None
        self.data = (np.arange(self.degree+1),np.arange(self.degree+1),self._ws)
        self.fitData(self.data[0],self.data[1])
        self._inits = self.data[:2]
    
    _fittypes=['spline']
    fittype = 'spline'
    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        from scipy.interpolate import UnivariateSpline
        
        #if the spline has never been fit to non-init data, set s appropriately
        if self._inits is not None and not (np.all(self._inits[0] == x) and np.all(self._inits[1] == y)):
            self.s = len(x)
            self._inits = None
        
        self.spline = UnivariateSpline(x,y,s=self.s,k=self.degree,w=kwargs['weights'] if 'weights' in kwargs else None)
        
        self._olds = self.s
        self._oldd = self.degree
        
        return np.array([self.s,self.degree])
        
    def fitData(self,x,y,**kwargs):
        self._oldd = self._olds = None
        
        if 'savedata' in kwargs and not kwargs['savedata']:
            raise ValueError('data must be saved for spline models')
        else:
            kwargs['savedata']=True
            
        if 'weights' in kwargs:
            self._ws = kwargs['weights']
        else:
            self._ws = None
            
        sorti = np.argsort(x)    
        return super(SmoothSplineModel,self).fitData(x[sorti],y[sorti],**kwargs)
    
    def f(self,x,s=2,degree=3):        
        if self._olds != s or self._oldd != degree:
            xd,yd,weights = self.data
            self.fitSpline(xd,yd,weights=weights)
        
        return self.spline(x)
    
    @property
    def rangehint(self):
        xd = self.data[0]
        return np.min(xd),np.max(xd)
    
    
class InterpolatedSplineModel(DatacentricModel1DAuto):
    """
    this uses a B-spline as a model for the function.  Note that by
    default the degree is left alone when fitting, as this model
    always fits the points perfectly.
    
    the scipy.interpolate.InterpolatedUnivariateSpline class is used to
    do the calculation (in the "spline" attribute) 
    """
    def __init__(self):
        super(InterpolatedSplineModel,self).__init__()
        
        self._oldd=self._olds=self._ws=None
        self.data = (np.arange(self.degree+1),np.arange(self.degree+1),self._ws)
        self.fitData(self.data[0],self.data[1])
            
            
    _fittypes = ['spline']
    fittype = 'spline'
    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        from scipy.interpolate import InterpolatedUnivariateSpline
        
        self.spline = InterpolatedUnivariateSpline(x,y,w=kwargs['weights'] if 'weights' in kwargs else None,k=self.degree)
        
        self._oldd = self.degree
        
        return np.array([self.degree])
        
    def fitData(self,x,y,**kwargs):
        self._oldd=None
        if 'savedata' in kwargs and not kwargs['savedata']:
            raise ValueError('data must be saved for spline models')
        else:
            kwargs['savedata']=True
            
        if 'weights' in kwargs:
            self._ws = kwargs['weights']
        else:
            self._ws = None
            
        sorti = np.argsort(x)    
        return super(InterpolatedSplineModel,self).fitData(x[sorti],y[sorti],**kwargs)
    
    def f(self,x,degree=3):        
        if self._oldd != degree:
            xd,yd,weights = self.data
            self.fitSpline(xd,yd,weights=weights)
        
        return self.spline(x)
    
    @property
    def rangehint(self):
        xd = self.data[0]
        return np.min(xd),np.max(xd)
    
class _KnotSplineModel(DatacentricModel1DAuto):
    """
    this uses a B-spline as a model for the function.  The knots
    parameter specifies the number of INTERIOR knots to use for the
    fit 
    
    locmethod can be:
    'cdf':the locations of the knots will be determined 
    by evenly sampling the cdf of the x-points
    'even':the knots are evenly spaced in x
    
    the scipy.interpolate.UnivariateSpline class is used to
    do the calculation (in the "spline" attribute) 
    """
    def __init__(self):
        super(_KnotSplineModel,self).__init__()
        
        self._ws = None
        
        self.data = (np.arange(self.degree+self.nknots+1),np.arange(self.degree+self.nknots+1),self._ws)
    
    @abstractmethod        
    def f(self,x):
        raise NotImplemetedError
    
    _fittypes = ['spline']
    fittype = 'spline'
    
    @abstractmethod    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):        
        from scipy.interpolate import LSQUnivariateSpline
        
        self.spline = LSQUnivariateSpline(x,y,t=self.iknots,k=int(self.degree),w=kwargs['weights'] if 'weights' in kwargs else None)
        
    def fitData(self,x,y,**kwargs):
        self._oldd=self._olds=None
        if 'savedata' in kwargs and not kwargs['savedata']:
            raise ValueError('data must be saved for spline models')
        else:
            kwargs['savedata']=True
            
        if 'weights' in kwargs:
            self._ws = kwargs['weights']
        else:
            self._ws = None
            
        sorti = np.argsort(x)    
        return super(_KnotSplineModel,self).fitData(x[sorti],y[sorti],**kwargs)
    
    @property
    def rangehint(self):
        xd = self.data[0]
        return np.min(xd),np.max(xd)

class UniformKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self._oldk = self._oldd = None
        super(UniformKnotSplineModel,self).__init__()
        self.fitData(self.data[0],self.data[1])
    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):
        self.iknots = np.linspace(x[0],x[-1],self.nknots+2)[1:-1]
        
        super(UniformKnotSplineModel,self).fitSpline(x,y,fixedpars,**kwargs)
        
        self._oldk = self.nknots
        self._oldd = self.degree
        
        return np.array([self.nknots,self.degree])
    
    def f(self,x,nknots=3,degree=3):
        if self._oldk != nknots or self._oldd != degree:
            xd,yd,weights = self.data
            self.fitSpline(xd,yd,weights=weights)
        
        return self.spline(x)
    
    

class UniformCDFKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self._oldk = self._oldd = None
        super(UniformCDFKnotSplineModel,self).__init__()
        self.fitData(self.data[0],self.data[1])
    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):
        cdf,xcdf = np.histogram(x,bins=max(10,max(2*self.nknots,int(len(x)/10))))
        mask = cdf!=0
        cdf,xcdf = cdf[mask],xcdf[np.hstack((True,mask))]
        cdf = np.hstack((0,np.cumsum(cdf)/np.sum(cdf)))
        self.iknots = np.interp(np.linspace(0,1,self.nknots+2)[1:-1],cdf,xcdf)
        
        super(UniformCDFKnotSplineModel,self).fitSpline(x,y,fixedpars,**kwargs)
        
        self._oldk = self.nknots
        self._oldd = self.degree
        
        return np.array([self.nknots,self.degree])
    
    def f(self,x,nknots=3,degree=3):
        if self._oldk != nknots or self._oldd != degree:
            xd,yd,weights = self.data
            self.fitSpline(xd,yd,weights=weights)
        
        return self.spline(x)

class SpecifiedKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self.nknots = self.__class__._currnparams
        self._oldd = None #this will force a fit at first call
        super(SpecifiedKnotSplineModel,self).__init__()
        
        self.setKnots(np.linspace(-1,1,self.nknots))
    
    def fitSpline(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        self.iknots = np.array([v for k,v in self.pardict.iteritems() if k.startswith('k')])
        self.iknots.sort()
        
        super(SpecifiedKnotSplineModel,self).fitSpline(x,y,fixedpars,**kwargs)
        
        self._oldd = self.degree
        
        retlist = list(self.iknots)
        retlist.insert(0,self.degree)
        return np.array(retlist)
    
    def setKnots(self,knots):
        if len(knots) != self.nknots:
            raise ValueError('provided knot sequence does not match the number of parameters')
        for i,k in enumerate(knots):
            setattr(self,'k'+str(i),k)
            
    def getKnots(self):
        ks = []
        for i in range(self.nknots):
            pn = 'k' + str(i)
            ks.append(getattr(self,pn))
        return np.array(ks)
    
    paramnames = 'k'
    
    degree=3 #default cubic
    def f(self,x,degree,*args):
        #TODO:faster way to do the arg check?
        if self._oldd != degree or np.any(self.iknots != np.array(args)):
            xd,yd,weights = self.data
            self.fitSpline(xd,yd,weights=weights)
        
        return self.spline(x)
    
    
    
class NFWModel(FunctionModel1DAuto):
    """
    A Navarro, Frenk, and White 1996 profile. Equivalent to
    :class:`AlphaBetaGammaModel` with alpha,beta,gamma = (1,3,1)
    
    where relevant, units have r in kpc and rho in Msun kpc^-3
    """
    xaxisname = 'r'
    yaxisname = 'rho'
        
    def f(self,r,rho0=1,rc=1):
        return rho0*rc*rc*rc*((r+rc)**(-2))*(r**-1)
    
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
    
    def integrateSpherical(self,lower,upper,*args,**kwargs):
        """
        NFW Has an analytic form for the spherical integral - if the lower 
        is not 0 or or if the keyword 'numerical' is True, this function will
        fall back to FunctionModel1D.integrateSpherical 
        """        
        if kwargs.pop('numerical',False):
            return FunctionModel1D.integrateSpherical(self,*args,**kwargs)
        else:
            x=upper/self.rc
            return 4*pi*self.rho0*self.rc**3*(np.log(1+x)-x/(1+x))
        
    def setC(self,c,Rvir=None,Mvir=None):
        """
        sets the model parameters to match a given concentration given
        virial radius and mass.
        
        If Rvir or Mvir is None, the other will be set using the z=0  
        Mvir_to_Rvir or Rvir_to_Mvir functions.  If both are None, a ValueError
        will be raised
        """
        if Rvir is None and Mvir is None:
            raise ValueError('need to specify Rvir or Mvir')
        elif Rvir is None:
            Rvir = self.Mvir_to_Rvir(Mvir,z=0)
        elif Mvir is None:
            Mvir = self.Rvir_to_Mvir(Rvir,z=0)
        
        self._c = c
        self.rc = Rvir/c
        
        self.rho0 = 1
        a0 = self.integrateSpherical(0,Rvir)
        self.rho0 = Mvir/a0
        
    def getC(self,z=0,rvir='delta',usesavedc=True):
        """
        returns the concentration (rvir/rc) if set, otherwise, calculates
        it from the virial radius assuming present cosmology at the specified 
        redshift
        
        `rvir` is the format that goes into  `getRv`
        
        `usesavedc` determines if a previously set c should override the 
        cosmologically-derived value
        """
        if hasattr(self,'_c') and usesavedc:
            return self._c
        else:
            return self.getRv(z,rvir)/self.rc
        
    def getRhoMean(self,r):
        """
        computes the mean density within the requested radius, i.e.
        M(<r)/(4/3 pi r^3)
        """
        vol = 4*pi*r**3/3
        return self.integrateSpherical(0,r)/vol
            
    def getRv(self,z=0,rvir='delta'):
        """
        get the virial radius at a given redshift, with rvir choosing the 
        type of virial radius - if 'delta' it will be computed from the 
        cosmology, otherwise, it should be a multiple of the critical 
        density
        
        units in kpc for mass in Msun
        """
        from scipy.optimize import newton
        from ..constants import get_cosmology,Ms
        cosmo = get_cosmology()
        
        if rvir == 'delta':
            overden = cosmo.deltavir(z)
        else:
            overden = float(rvir)
        
        try:
            rhoC = cosmo.rhoC(z,'cosmological')*1e-9 #Mpc^-3->kpc^-3
            rhov = overden*rhoC
        except:
            raise ValueError('current cosmology does not support critical density')
        oldcall = self.getCall()
        try:
            self.setCall('getRhoMean')
            return self.inv(rhov,self.rc)
        finally:
            self.setCall(*oldcall)
    
    def getMv(self,z=0):
        """
        gets the mass within the virial radius (for a particular redshift)
        """
        rv = self.getRv(z)
        return self.integrateSpherical(0,rv)
    
    @staticmethod
    def Rvir_to_Mvir(Rvir,z=0):
        """
        from Bullock '01
        M_sun,kpc
        """
        from ..constants import get_cosmology
        c = get_cosmology()
        
        return 1e12/c.h*(c.omega*c.deltavir(z)/97.2)*(Rvir*(1+z)/203.4/c.h)**3
    
    @staticmethod
    def Mvir_to_Rvir(Mvir,z=0):
        """
        from Bullock '01
        M_sun,kpc
        """
        from ..constants import get_cosmology
        c = get_cosmology()
        
        return 203.4/c.h*(c.omega*c.deltavir(z)/97.2)**(-1/3)*(Mvir/1e12/c.h)**(1/3)/(1+z)
    
    @staticmethod
    def Mvir_to_Vvir(Mvir,z=0):
        """
        from Bullock '01
        km/s,M_sun
        """
        from ..constants import get_cosmology
        c = get_cosmology()
        
        return 143.8*(c.omega*c.deltavir(z)/97.2)**(1/6)*(Mvir/1e12/c.h)**(1/3)*(1+z)**0.5
    
    @staticmethod
    def Vvir_to_Mvir(Vvir,z=0):
        """
        from Bullock '01
        km/s,M_sun
        """
        from ..constants import get_cosmology
        c = get_cosmology()
        
        return (c.omega*c.deltavir(z)/97.2)**-0.5*(1+z)**-1.5*c.h*1e12*(Vvir/143.8)**3
    
    @staticmethod
    def Vvir_to_Vmax(Vvir):
        """
        from Bullock '01
        good from 80-1200 km/s in vvir
        """
        return (Vvir/0.468)**(1/1.1)
    
    @staticmethod
    def Vmax_to_Vvir(Vmax):
        """
        from Maller & Bullock '04
        good from 80-1200 km/s in vvir
        """
        return (0.468)*Vmax**1.1
    
    @staticmethod
    def Mvir_to_Cvir(Mvir,z=0):
        """
        from Maller & Bullock '04
        M_sun
        """
        cvir = 9.6*(Mvir/1e13)**-.13*(1+z)**-1
        return cvir if cvir>4 else 4.0
    
    @staticmethod
    def Cvir_to_Mvir(Cvir,z=0):
        """
        from Maller & Bullock '04
        M_sun
        """
        return 1e13*((Cvir/9.6)*(1+z))**(-1.0/.13)
    
    @staticmethod
    def Mvir_to_Vmax(Mvir,z=0):
        return NFWModel.Vvir_to_Vmax(NFWModel.Mvir_to_Vvir(Mvir,z))
    
    @staticmethod
    def Vmax_to_Mvir(Vmax,z=0):
        return NFWModel.Vvir_to_Mvir(NFWModel.Vmax_to_Vvir(Mvir),z)
    
    @staticmethod
    def create_Mvir(Mvir,z=0):
        """
        generate a new NFWModel with the given Mvir using the Mvir_to_Cvir
        function to get the scaling
        """
        m = NFWModel()
        Cvir = m.Mvir_to_Cvir(Mvir,z)
        Rvir = m.Mvir_to_Rvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
    @staticmethod
    def create_Rvir(Rvir,z=0):
        """
        generate a new NFWModel with the given Mvir using the Mvir_to_Cvir and
        Rvir_to_Mvir functions to get the scaling
        """
        m = NFWModel()
        Mvir = m.Rvir_to_Mvir(Rvir,z)
        Cvir = m.Mvir_to_Cvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
    @staticmethod
    def create_Cvir(Cvir,z=0):
        """
        generate a new NFWModel with the given Mvir using the Cvir_to_Mvir
        function to get the scaling
        """
        m = NFWModel()
        Mvir = m.Cvir_to_Mvir(Cvir,z)
        Rvir = m.Mvir_to_Rvir(Mvir,z)
        m.setC(Cvir,Rvir=Rvir,Mvir=Mvir)
        return m
    
class AlphaBetaGammaModel(FunctionModel1DAuto):
    """
    This is a model of the form
    A*(r/rs)**-gamma*(1+(r/rs)**alpha)**((gamma-beta)/alpha) where gamma is the
    inner log slope, beta is the outer log slope, and alpha controls the
    transition region. 
    
    If a logarithmic version of the profile is desired, use::
    
    >>> m = AlphaBetaGammaModel()
    >>> m.setCall(xtrans='pow',ytrans='log')
    
    In this case the offset is given by log10(A) and the logaritmhic scale
    radius is log10(rs).
    """
    
    def f(self,r,rs=1,A=1,alpha=1,beta=2,gamma=1):
        ro = r/rs
        return A*(ro)**-gamma*(1+(ro)**alpha)**((gamma-beta)/alpha)
    
    @property
    def rangehint(self):
        return self.rs/1000,1000*self.rs
    
    

class PlummerModel(FunctionModel1DAuto):
    xaxisname = 'r'
    
    def f(self,r,rp=1.,M=1.):
        return 3*M/(4.*pi*rp**3)*(1+(r/rp)**2)**-2.5
    
    @property
    def rangehint(self):
        return 0,self.rp*2

class King2DrModel(FunctionModel1DAuto):    
    xaxisname = 'r'
    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        return A*rcsq*((r*r+rcsq)**-0.5 - (rt*rt+rcsq)**-0.5)**2
    
    @property
    def rangehint(self):
        return 0,self.rt
    
class King3DrModel(FunctionModel1DAuto):
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

class SchecterMagModel(FunctionModel1DAuto):
    xaxisname = 'M'
    yaxisname = 'phi'
    
    def f(self,M,Mstar=-20.2,alpha=-1,phistar=1.0857362047581294):
        from numpy import log,exp
        x=10**(0.4*(Mstar-M))
        return 0.4*log(10)*phistar*(x**(1+alpha))*exp(-x)
    
    @property
    def rangehint(self):
        return self.Mstar-3,self.Mstar+3
    
class SchecterLumModel(FunctionModel1DAuto):
    xaxisname = 'L'
    yaxisname = 'phi'
    
    def f(self,L,Lstar=1e10,alpha=-1.0,phistar=1.0):
        #from .phot import lum_to_mag as l2m
        #M,Mstar=l2m(L),l2m(Lstar)
        #return SchecterModel.f(self,M,Mstar,alpha,phistar)
        x = L/Lstar
        return phistar*(x**alpha)*np.exp(-x)/Lstar
    #TODO:check to make sure this is actually the right way
    
    @property
    def rangehint(self):
        return self.Lstar/10,self.Lstar*10
    
    @property
    def rangehint(self):
        return self.Lstar/3,self.Lstar*3
        
class EinastoModel(FunctionModel1DAuto):
    """
    Einasto profile ln(rho/A) = (-2/alpha)((r/rs)^alpha - 1) . :attr:`A` is the
    density where the log-slope is -2, and :attr:`rs` is the corresponding
    radius. The `alpha` parameter defaults to .17 as suggested in Navarro et al.
    2010
    """
    xaxisname = 'r'
    yaxisname = 'rho'
    
    def f(self,r,A=1,rs=1,alpha=.17):
        return A*np.exp((-2/alpha)*((r/rs)**alpha - 1))

class SersicModel(FunctionModel1DAuto):
    """
    Sersic surface brightness profile:
    Ae*exp(-b_n[(R/Re)^(1/n)-1])
    
    Ae is the value at the effective radius re
    """
    xaxisname = 'r'
    
    def f(self,r,Ae=1,re=1,n=4):
        #return EinastoModel.f(self,r,A,rs,1/n)
        #return A*np.exp(-(r/rs)**(1/n))
        return Ae*np.exp(-self.bn(n)*((r/re)**(1.0/n)-1))
    
    @property
    def rangehint(self):
        return 0,2*self.re
    
    def _getA0(self):
        return self.f(0,self.Ae,self.re,self.n)
    def _setA0(self,val):
        self.Ae *= val/self.f(0,self.Ae,self.re,self.n)
    A0 = property(_getA0,_setA0,doc='value at r=0')
    
    def getBn(self,usecache=True):
        """
        Calls `bn` function on the current sersic index.
        """
        return self.bn(self.n,usecache)
    
    _bncache={}
    _bnpoly1=np.poly1d([-2194697/30690717750,131/1148175,46/25515,4/405,-1/3])
    _bnpoly2=np.poly1d([13.43,-19.67,10.95,-0.8902,0.01945])
    @staticmethod
    def bn(n,usecache=True):
        """
        bn is used to get the half-light radius.  If n is 
        None, the current n parameter will be used
        
        the form is a fit from MacArthur, Courteau, and Holtzman 2003 
        and is claimed to be good to ~O(10^-5)
        
        if usecache is True, the cache will be searched, if False it will
        be saved but not used, if None, ignored
        """
        n = float(n) #sometimes 0d array gets passed in
        
        if n  in SersicModel._bncache and usecache:
            val = SersicModel._bncache[n]
        else:
            val = (2*n+SersicModel._bnpoly1(1/n)) if n>0.36 else SersicModel._bnpoly2(n)
            if usecache is not None:
                SersicModel._bncache[n] = val
        return val
        
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
    def f(self,r,Ae=1,re=1):
        return SersicModel.f(self,r,Ae,re,4)
    
class JaffeModel(FunctionModel1DAuto):
    """
    Jaffe (1983) profile as defined in Binney & Tremain as
    (A/4pirj^3)(rj^4/(r^2(r+rj)^2))
    where A is the total mass enclosed 
    
    """
    def f(self,r,A=1,rj=1):
        return (A/4/pi*rj**-3)*rj**4*r**-2*(r+rj)**-2
    
    @property
    def rangehint(self):
        return 0,self.rj*3
    
class HernquistModel(FunctionModel1DAuto):
    """
    Hernquist (1990) profile described as
    (A/2pi)(r0/r)(1/((r+r0)^3))
    where A is the total mass enclosed
    """
    def f(self,r,A=1,r0=1):
        return (A/2/pi)*(r0/r)*(r+r0)**-3
    
    @property
    def rangehint(self):
        return 0,self.r0*3



class MaxwellBoltzmannModel(FunctionModel1DAuto):
    xaxisname = 'v'
    
    from ..constants import me #electron
    def f(self,v,T=273,m=me):
        from ..constants import kb,pi
        return (m/(2*pi*kb*T))**0.5*np.exp(-m*v*v/2/kb/T)
    
    @property
    def rangehint(self):
        from ..constants import kb,c
        return 0,min(3*(2*kb*self.T/self.m)**0.5,c)
    
class MaxwellBoltzmannSpeedModel(MaxwellBoltzmannModel):
    from ..constants import me #electron
    xaxisname = 'v'
    
    def f(self,v,T=273,m=me):
        from ..constants import kb,pi
        return 4*pi*v*v*(m/(2*pi*kb*T))**1.5*np.exp(-m*v*v/2/kb/T)
    
    @property
    def rangehint(self):
        from ..constants import kb,c
        return 0,min(3*(2*kb*self.T/self.m)**0.5,c)
    
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

#<-------------------------------------- 2D models ---------------------------->
class Gaussian2DModel(FunctionModel2DScalarAuto):
    _fcoordsys='cartesian'
    def f(self,inarr,A=1,sigx=1,sigy=1,mux=0,muy=0):
        x,y = inarr
        xo = x-mux
        xdenom = 2*sigx*sigx
        yo = y-muy
        ydenom = 2*sigy*sigy
        return A*np.exp(-xo*xo/xdenom-yo*yo/ydenom)
    
    @property
    def rangehint(self):
        mux,muy = self.mux,self.muy
        sigx,sigy = self.sigx,self.sigy
        return (mux-4*sigx,mux+4*sigx,muy-4*sigy,muy+4*sigy)
        
class EdgeOnDiskModel(FunctionModel2DScalarAuto):
    _fcoordsys='cartesian'
    def f(self,inarr,A=1,l=2,h=1,pa=0):
        x,y = inarr
        
        sinpa,cospa = np.sin(pa),np.cos(pa)
        xr = cospa*x+sinpa*y
        yr = -sinpa*x+cospa*y
        return A*np.exp(-np.abs(xr/l)-np.abs(yr/h))
    
    @property
    def rangehint(self):
        bigscale = max(self.h,self.l)
        return (-2*bigscale,2*bigscale,-2*bigscale,2*bigscale)
    
class InclinedDiskModel(FunctionModel2DScalarDeformedRadial):
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
    A bulge modeled as a radially symmetric sersic profile
    (by default the sersic index is kept fixed - to remove
    this behavior, set the fixedpars attribute to None)
    """
    
    fixedpars = ('n',)    
    def __init__(self,Ae=1,re=1,n=4):
        super(RoundBulgeModel,self).__init__('sersic')
        self.Ae = Ae
        self.re = re
        self.n = n
    
#<-------------------------------Other Models---------------------------------->
    
class Plane(FunctionModel):
    """
    Generates a plane that follows the form
    d = a*x+b*y+c*z (e.g. (a,b,c) is the normal vector and 
    d/a, b ,or c are the intercepts
    """    
    def __init__(self,varorder='xyz',vn=(1,0,0),wn=(0,1,0),origin=(0,0,0)):
        self.varorder = varorder
        self.vn=vn
        self.wn=wn
        self.origin = origin
    
    def _getvaro(self):
        return self._varo
    def _setvaro(self,val):
        if val == 'xyz':
            self._f = self._fxyz
        elif val == 'yxz':
            self._f = self._fyxz
        elif val == 'xzy':
            self._f = self._fxzy
        elif val == 'zxy':
            self._f = self._fzxy
        elif val == 'yzx':
            self._f = self._fyzx
        elif val == 'zyx':
            self._f = self._fzyx
        else:
            raise ValueError('unrecognized variable order')
        self._varo = val
    varorder = property(_getvaro,_setvaro)
    
    def _getvn(self):
        return self._vn
    def _setvn(self,val):
        vn = np.array(val)
        if vn.shape != (3,):
            raise ValueError('vn must be a length-3 vector')
        self._vn = vn
    vn = property(_getvn,_setvn,doc='3D vector to project on to plane to get 2D basis vector 1')
    
    def _getwn(self):
        return self._wn
    def _setwn(self,val):
        wn = np.array(val)
        if wn.shape != (3,):
            raise ValueError('wn must be a length-3 vector')
        self._wn = wn
    wn = property(_getwn,_setwn,doc='3D vector to project on to plane to get 2D basis vector 2')

    def _getorigin(self):
        n = self.n
        scale = (self.d - np.dot(self._origin,n))/np.dot(n,n)
        return self._origin + scale*n
    def _setorigin(self,val):
        val = np.array(val,copy=False)
        if val.shape == (2,):
            self._origin = self.unproj(*val)[:,0]
        elif val.shape == (3,):
            self._origin = val
        else:
            raise ValueError('invalid shape for orign - must be 2-vector or 3-vector')
    origin = property(_getorigin,_setorigin)
    
    @property
    def n(self):
        """
        non-normalized unit vector
        """
        return np.array((self.a,self.b,self.c))
    
    @property
    def nhat(self):
        """
        normalized unit vector
        """
        n = np.array((self.a,self.b,self.c))
        return n/np.linalg.norm(n)
    
    def f(self,x,a=0,b=0,c=1,d=0):
        x = np.array(x,copy=False)
        shp = x.shape
        if len(shp) > 2: 
            x = x.reshape(2,np.prod(shp[1:]))
            y = self._f(x,a,b,c,d)
            return y.reshape(shp[1:])
        else:
            return self._f(x,a,b,c,d)
    
    def _fxyz(self,v,a,b,c,d):
        M = np.matrix([(a/c,b/c)])
        return d/c-(M*v).A
    def _fyxz(self,v,a,b,c,d):
        M = np.matrix((b/c,a/c))
        return d/c-(M*v).A
    def _fxzy(self,v,a,b,c,d):
        M = np.matrix((a/b,c/b))
        return d/b-(M*v).A
    def _fzxy(self,v,a,b,c,d):
        M = np.matrix((c/b,a/b))
        return d/b-(M*v).A
    def _fyzx(self,v,a,b,c,d):
        M = np.matrix((b/a,c/a))
        return d/a-(M*v).A
    def _fzyx(self,v,a,b,c,d):
        M = np.matrix((c/a,b/a))
        return d/a-(M*v).A
    
    def fitData(self,x,y,z,w=None):
        """
        least squares fit using the output variable as the dependent
        """
        from scipy.optimize import leastsq
        #reorder vars to get the right fitter
        x,y,z = eval(','.join(self._varo))
        
        xy = np.array([x,y],copy=False)
        if w is None:
            f = lambda v:(self.f(xy,*v)-z).ravel()
        else:
            f = lambda v:(self.f(xy,*v)-z).ravel()*w**0.5
        
        res = leastsq(f,self.parvals,full_output=1)
        self.lastfit = res
        
        self.parvals = res[0]
        return res[0]
    
    def distance(self,x,y,z):
        """
        compute the distance of a set of points in the 3D space from 
        the plane
        """
        shp = list(x.shape)
        x = np.array(x,copy=False).ravel()
        y = np.array(y,copy=False).ravel()
        z = np.array(z,copy=False).ravel()
        p = np.c_[x,y,z]
        
        return (np.dot(p,self.n)+self.d).reshape(shp)
    
    def proj(self,x,y,z):
        """
        project points onto the plane from the 3D space
        
        returns a 2 x N aray 
        """
        n = self.nhat
        
        vn = np.cross(np.cross(n,self.vn),n)
        wn = np.cross(np.cross(n,self.vn),n)
        
        shp = list(x.shape)
        x = np.array(x,copy=False).ravel()
        y = np.array(y,copy=False).ravel()
        z = np.array(z,copy=False).ravel()
        p = np.c_[x,y,z] - self.origin
        
        shp.insert(0,2)
        return (np.c_[np.dot(p,vn),np.dot(p,wn)].T).reshape(shp)
    
    def unproj(self,v,w):
        """
        extract points from the plane back into the 3D space
        
        returns a 3 x N array
        """
        n = self.nhat
        
        vn = np.cross(np.cross(n,self.vn),n)
        wn = np.cross(np.cross(n,self.vn),n)
        
        shp = list(v.shape)
        v = np.array(v,copy=False).ravel()
        w = np.array(w,copy=False).ravel()
        
        shp.insert(0,3)
        return (v*vn+w*wn + self.origin).reshape(shp)
    
    def plot3d(self,data=np.array([(-1,1),(-1,1),(-1,1)]),n=10,
               showdata=False,clf=True,**kwargs):
        """
        data should be 3 x N
        """
        import enthought.mayavi.mlab as M
        data = np.array(data,copy=False)
        
        xi,xx = data[0].min(),data[0].max()
        yi,yx = data[1].min(),data[1].max()
        x,y = np.meshgrid(np.linspace(xi,xx,n),np.linspace(yi,yx,n))
        
        if clf:
            M.clf()
        
        if 'color' not in kwargs:
            kwargs['color']=(1,0,0)
        if 'opacity' not in kwargs:
            kwargs['opacity'] = 0.5
            
        M.mesh(x,y,self([x,y]),**kwargs)
        if showdata:
            from operator import isMappingType
            if isMappingType(showdata):
                M.points3d(*data,**showdata)
            else:
                M.points3d(*data)
    
#register everything in this module
from inspect import isclass
for o in locals().values():
    if isclass(o) and not o.__name__.startswith('_') and issubclass(o,ParametricModel):
        if 'FunctionModel' not in o.__name__ and 'CompositeModel' not in o.__name__:
            register_model(o)
            
#cleanup namespace
del isclass,o
