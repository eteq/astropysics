from __future__ import division
from math import pi
import numpy as np


class _Parameter(object):
    
    def __init__(self,name,defaultvalue):
        self.name = name
        self.defaultvalue = defaultvalue
    
    def __get__(self,obj,owner):
        try:
            return getattr(obj,'_Parameter__'+self.name)
        except AttributeError: #if data doesn't exist, initialize to default and return it
             setattr(obj,'_Parameter__'+self.name,self.defaultvalue)
             return self.defaultvalue
    
    def __set__(self,obj,value):
        if not np.isscalar(value):
            raise ValueError('Parameters must be scalar values')
        setattr(obj,'_Parameter__'+self.name,value)
    
    def __delete__(self,obj):
        raise AttributeError("can't delete a parameter")

class _FuncMeta(type):
    def __init__(cls,name,bases,dct):
        #called on import of astro.models
        from inspect import getargspec
        
        type.__init__(name,bases,dct)
        
        args,varargs,varkw,defaults=getargspec(dct['f'])
        if varkw is not None:
            raise SyntaxError("can't have kwargs in model function")
        
        #remove self if present, and the first argument, which is the data vector
        try:
            args.remove('self')
        except ValueError: # means the function is static or class method
            try:
                args.remove('cls')
            except ValueError: # means the function is static, which is still ok
                pass
        del args[0] #remove data vector
        
        if defaults is None:
            defaults=[]
        if len(args)!=len(defaults):
            defaults=list(defaults)
            for i in range(len(args)-len(defaults)): defaults.insert(0,1)
            
        for arg,d in zip(args,defaults):
            setattr(cls,arg,_Parameter(arg,d))
            
        if varargs is not None:
            cls.__args = None
            cls.__statargs = list(args)
        else:
            cls.__args = tuple(args)
        
            
        cls.params = property(lambda cls:cls.__args)
        def _set_parvals(cls,newpars):
            if len(newpars)>len(cls.__args):
                raise ValueError('too many new parameters')
            for a,v in zip(cls.__args,newpars):
                setattr(cls,a,v)        
        def _set_pardict(cls,newpardict):
            for k in newpardict.keys():
                if k not in cls.__args:
                    raise KeyError('No parameter %s found'%k)
            for k,v in newpardict.iteritems():
                setattr(cls,k,v)
        cls.parvals = property(lambda cls:[getattr(cls,a) for a in cls.__args],fset=_set_parvals)
        cls.pardict = property(lambda cls:dict([(a,getattr(cls,a)) for a in cls.__args]),fset=_set_pardict)
        
    def __call__(cls,*args,**kwargs):
        if cls.__args is None: #this is the case for *args in the function
            args=list(args)
            try:
                nparams=args.pop(0) if 'nparams' not in kwargs else kwargs.pop('nparams')
            except IndexError:
                IndexError('No # of parameters found for variable-size function')
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls.__args])    
            obj = type.__call__(cls,**objkwargs) #object __init__ is called here
            pars = cls.__statargs
            del cls.__statargs
            for i in range(nparams):
                p='p%i'%i
                pars.append(p)
                setattr(obj,p,1) #default varargs to 1
            cls.__args = tuple(pars)
        else: #this is the case for fixed functions
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls.__args])
            obj = type.__call__(cls,**objkwargs) #object __init__ is called here
            
        obj.f(np.array([0]),*obj.parvals) #try once to check for a real function
            
        #set initial values
        if len(args)+len(kwargs) > len(obj.params):
            raise TypeError('Too many args and kwargs for the parameters')
        
        pars=list(obj.params)
        for k,v in kwargs.iteritems():
            del pars[pars.index(k)]
            setattr(obj,k,v)
        for p,v in zip(pars,args):
            setattr(obj,p,v)
        return obj

class FunctionModel(object):
    """
    This class is the base for Models with parameters based on python functions.
    
    Subclassing:
    The following method MUST be overridden in a subclass:
    *f(self,x,...)
    The following methods may be overridden for speed of some operations - pars
    should be accessed with self.pardict, self.parvals, or by properties/name :
    *integrate(self,lower,upper)
    *derivative(self,x,dx)
    *_fitCustom(x,y,fixedpars=(),**kwargs)
    
    The metaclass generates the following properties:
    Parameters for each of the inputs of the function (other than self/cls and the first)
    self.params: a tuple of the parameter names
    self.parvals: a list of the values in the parameters
    self.pardict: a dictionary of the parameter names and values
    """
    __metaclass__ = _FuncMeta
    
    defaultIntMethod = 'quad'
    
    def __call__(self,x):
        """
        call the function on the input x with the current parameters and return
        the result
        """
        return self.f(np.atleast_1d(x),*self.parvals)
        
        
    def fitData(self,x,y,method=None,fixedpars=(),contraction='sumsq',
                 updatepars=True,timer=None,**kwargs):
        """
        This will use the data to adjust the parameters to fit using any of a
        variety of methods from scipy.optimize
        
        method can be any of the optimize functions (except constrained or
        scalar optimizers) or 'custom'
        default is custom, or if not present, leastsq
        
        fixedpars is a sequence of parameter names to leave fixed
        
        contraction applies for all except the 'leastsq' method  (for which only
        'frac' is meaningful) and is the technique used to convert vectors to
        figures of merit.  this is composed of multiple string segments -
        one of:
        'sq','abs', or '' applied to each element of the vector
        one of:
        'sum','median','mean','': applied after the above
        optinally,
        'frac' to use the fractional version of the differnce vector
        if it is None/False, do nothing (function already outputs a scalar, and
        data is a scalar)
        
        kwargs go into the fitting function (don't alter 'full_output')
        
        returns the vector of the best fit parameters, and sets the parameters if
        updatepars is True
        
        the full output is available is self.lastFit
        
        see also:fitMCMC
        """
        from scipy import optimize as opt
        
        if method is None:
            #TODO: figure out a less hackish way than matching the docs
            method = 'leastsq' if self._fitCustom.__doc__ is FunctionModel._fitCustom.__doc__ else 'custom'
            
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        if timer:
            raise NotImplementedError
        
        
        ps=list(self.params)
        v=list(self.parvals) #initial guess
        if fixedpars:
            for p in fixedpars:
                i=ps.index(p)
                del ps[i]
                del v[i]
            
            #make a function of signature f(x,v) where v are the parameters to be fit
            pdict=dict([(p,getattr(self,p)) for p in fixedpars])
            def f(x,v):
                pdict.update(dict(zip(ps,v)))
                return self.f(x,**pdict)
        else:
            f=lambda x,v:self.f(x,*v)
        
        if method == 'custom':
            res = self._fitCustom(x,y,fixedpars=fixedpars,**kwargs) 
        elif method == 'leastsq':
            if 'frac' in contraction:
                g=lambda v,x,y:1-f(x,v)/y
            else:
                g=lambda v,x,y:y-f(x,v)
            res=opt.leastsq(g,v,(x,y),full_output=1,**kwargs)
        else:
            if not contraction:
                f=lambda v,x,y:y-f(x,v)
            else:
                if 'frac' in contraction:
                    if 'sq' in contraction:
                        def g1(v,x,y):
                            diff=1-f(x,v)/y
                            return diff*diff
                    elif 'abs' in contraction:
                        def g1(v,x,y):
                            diff=1-f(x,v)/y
                            return np.abs(diff)
                    else:
                        def g1(v,x,y):
                            diff=1-f(x,v)/y
                            return diff
                else:
                    if 'sq' in contraction:
                        def g1(v,x,y):
                            diff=y-f(x,v)
                            return diff*diff
                    elif 'abs' in contraction:
                        def g1(v,x,y):
                            diff=y-f(x,v)
                            return np.abs(diff)
                    else:
                        def g1(v,x,y):
                            diff=y-f(x,v)
                            return diff
                if 'sum' in contraction:
                    g=lambda v,x,y:np.sum(g1(v,x,y))
                elif 'mean' in contraction:
                    g=lambda v,x,y:np.mean(g1(v,x,y))
                elif 'median' in contraction:
                    g=lambda v,x,y:np.median(g1(v,x,y))
                
            if method == 'fmin':
                res=opt.fmin(g,v,(x,y),full_output=1)
            elif method == 'fmin_powell':
                res=opt.fmin_powell(g,v,(x,y),full_output=1)
            elif method == 'fmin_cg':
                #TODO:smartly include derivative
                res=opt.fmin_cg(g,v,args=(x,y),full_output=1)
            elif method == 'fmin_bfgs':
                #TODO:smartly include derivative
                res=opt.fmin_bfgs(g,v,args=(x,y),full_output=1)
            elif method == 'fmin_ncg':
                raise NotImplementedError
                #TODO:needs gradient and hessian
                opt.fmin_ncg
            elif method == 'anneal' or method == 'global':
                res=opt.anneal(g,v,args=(x,y),full_output=1,**kwargs)
            elif method == 'brute':
                raise NotImplementedError
                #TODO: set parrange smartly
                res=opt.brute(g,parrange,(x,y),full_output=1,**kwargs)
            else:
                raise ValueError('Unrecognzied method %s'%method)
            
        self.lastFit = res
        v=res[0] #assumes output is at least a tuple - needs "full_output=1 !"
        
        try:
            v[0]
        except IndexError: #only one parameter
            v=np.array([v])
            
        if updatepars:
            for par,newv in zip(ps,v):
                setattr(self,par,newv)
        
        return v
    
    def _fitCustom(x,y,fixedpars=(),**kwargs):
        """
        Must be overridden to use 'custom' fit technique
        
        should return a tuple where the first element gives the best fit
        parameter vector 
        """
        raise NotImplementedError('No custom fit function provided for this model')
    
    def stdData(self,x,y):
        """
        determines the standard deviation of the model from the supplied data
        """
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        return np.std(self(x)-y,ddof=len(self.params))
    
    def chi2Data(self,x,y,dy=None):
        """
        determines the chi-squared statistic for data with errors dy
        
        if dy is None, the bootstrap technique will be used to derive errors
        (possibly not correctly TODO:check)
        
        returns chi2,p-value
        """
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        n=len(y)
        m=len(self.params)
        dof=n-m
        
        if dy is None:
            ps=self.parvals
            try:
                d=self.bootstrapFits(x,y)
                dupper=dict([(k,(np.median(v)+np.std(v))) for k,v in d.iteritems()])
                dlower=dict([(k,(np.median(v)-np.std(v))) for k,v in d.iteritems()])
                self.pardict=dupper
                yu=self(x)
                self.pardict=dlower
                yl=self(x)
                dy=(yu+yl)/2.0-y
            finally:
                self.parvals=ps
        chi2=(y-self(x))/dy
        chi2=np.sum(chi2*chi2)
        
        from scipy.stats import chisqprob
        p=chisqprob(chi2,dof)
        return chi2,p
    
    def bootstrapFits(self,x,y,n=250,prefit=True,plothist=False,**kwargs):
        """
        uses the fitData function to fit the function and then uses the
        "bootstrap" technique to estimate errors - that is, use shuffling
        w/replacement
        
        x and y are the data to fit, n is the number of iterations to perform,
        prefit determines if the function is to be fit once before the
        bootstrapping begins, and kwargs go into self.fitData
        """
        
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        xyn=len(x)
        inds=(xyn*np.random.rand(n,xyn)).astype(int)
        xi,yi=x[inds],y[inds]
        vs=[]
        
        if prefit:
            self.fitData(x,y,**kwargs)
            
        for xd,yd in zip(xi,yi):
            vs.append(self.fitData(xd,yd,updatepars=False,**kwargs))
            
        d=dict([(p,v) for p,v in zip(self.params,np.array(vs).T)])
        
        if plothist:
            import matplotlib.pyplot as plt
            #want to generate as square as possible
            nr=nc=int(len(self.params)**0.5)
            if nr*nc < len(self.params):
                nr+=1
            if nr*nc < len(self.params):
                nc+=1
                
            for i,p in enumerate(self.params):
                print nr,nc
                plt.subplot(nr,nc,i+1)
                plt.hist(d[p],bins=n//20)
                plt.xlabel(p)
                plt.ylabel('N')
                plt.title('Median$=%3.3f$ $\\sigma=%3.3f$ Fit$=%3.3f$'%(np.median(d[p]),np.std(d[p]),self.pardict[p]))
        
        return d
        
    def fitMCMC(self,x,y):
        raise NotImplementedError
        
    def plot(self,lower=None,upper=None,n=100,integrate=None,clf=True,logplot='',powerx=False,powery=False,deriv=None,data=None,*args,**kwargs):
        """
        (Assumes 1D)
        plot the model function from lower to upper with n samples
        
        integrate controls whether or not to plot the integral of the function -
        if True, the base of the integral is taken to be lower, if False, upper,
        no integration if None, and otherwise, the base is the argument value
        
        extra args and kwargs go into the matplotlib plot function
        
        data is either an x,y pair of data points or a dictionary of kwargs into scatter
        
        logplot determines whether or not to plot on a log scale, and powerx and
        poweru determine if the model points should be used as powers (base 10)
        before plotting
        """
        from matplotlib import pyplot as plt
        from operator import isMappingType
        
        if (lower is None or upper is None) and data is None:
            raise ValueError("Can't choose limits for plotting without specifying or providing data")
        if data is not None:
            if lower is None:
                lower=np.min(data[0])
            if upper is None:
                upper=np.max(data[0])
        
        if 'x' in logplot:
            x = np.logspace(np.log10(lower),np.log10(upper),n)
        else:
            x = np.linspace(lower,upper,n)
        
        if powerx:
            x=10**x
            
        if integrate is None and not deriv:
            y = self(x)
        elif deriv:
            if integrate is not None:
                raise ValueError("can't do derivative and integral simultaneously")
            y=np.array([self.derivative(v,np.abs(upper-lower)/(10*n)) for v in x])
        else: #integrate is not None
            if integrate is True:
                base = lower
            elif integrate is False:
                base = upper
            else:
                base = float(integrate)
            y=np.array([self.integrate(base,v) for v in x])
            
        if powery:
            y= 10**y
        
        if clf:
            plt.clf()
        if 'x' in logplot and 'y' in logplot:
            plt.loglog(x,y,*args,**kwargs)
        elif 'x' in logplot:
            plt.semilogx(x,y,*args,**kwargs)
        elif 'y' in logplot:
            plt.semilogy(x,y,*args,**kwargs)
        else:
            plt.plot(x,y,*args,**kwargs)
        
        if data:
            if isMappingType(data):
                if clf and 'c' not in data:
                    data['c']='g'
                plt.scatter(**data)
            else:
                if clf and len(data)<4:
                    data=list(data)
                    if len(data) < 3:
                        data.append(20)
                    data.append('r')
                plt.scatter(*data)
                
        
    #Can Override:
    def integrate(self,lower,upper,method=None,n=None,jac=None,**kwargs):
        """
        (Assumes 1D)
        This will numerically integrate the model function from lower to upper
        using scipy.integrate techniques
        
        method can be specified (a function from scipy.integrate) or if None,
        self.defaultIntMethod will be used
        
        n is either the number of data values or the order for ordered methods,
        or if sequence,will be used as the integration x-values (and lower/upper
        are ignored)
        
        jac is the jacobian factor as a function f(x,*params)
        
        returns value
        some methods will also store attributes "lastIntError" and "lastIntInfo"
        
        if subclassing, this should be overwritten with a function of the form
        integrate(lower,upper,[args])
        """
        import scipy.integrate as itg
        if method is None:
            method=self.defaultIntMethod
        
        e,d=None,None
        
        ps=tuple(self.parvals)
        if jac:
            def f(x,*pars):
                return jac(x,*pars)*self.f(x,*pars)
        else:
            f=self.f
        
        if method=='quad':
            t=itg.quad(f,lower,upper,args=ps,full_output=1,**kwargs)
            if len(t) == 2:
                v,e = t
            elif len(t) == 3:
                v,e,d = t
            else:
                raise ValueError("problem with return value of scipy.integrate.quad")
        #use these for 2d and 3d
        #elif method=='dblquad':
        #    raise NotImplementedError
        #elif method=='tplquad':
        #    raise NotImplementedError
        elif method=='fixed_quad':
            v=itg.fixed_quad(f,lower,upper,ps,5 if n is None else n,**kwargs)
        elif method=='quadrature':
            v,e=itg.quadrature(f,lower,upper,ps,**kwargs)
        elif method=='romberg':
            v=itg.romberg(f,lower,upper,ps,**kwargs)
        else: #sampled techniques
            if n is None:
                n=100
            if np.isscalar(n):
                x=np.linspace(lower,upper,n)
            else:
                x=np.array(n)
            y=f(x,*ps)
            if method=='trapz':
                v=itg.trapz(y,x,**kwargs)
            elif method=='cumtrapz':
                v=itg.cumtrapz(y,x,**kwargs)
            elif method=='simps':
                v=itg.simps(y,x,**kwargs)
            elif method=='romb':
                v=itg.simps(y,np.convolve(x,[1,-1],mode='same').mean(),**kwargs)
            else:
                raise ValueError('unrecognized integration method')
        
        if e is not None:
            self.lastIntError = e
        if d is not None:
            self.lastIntInfo = d
        return v
    
    def derivative(self,x,dx=1):
        """
        the derivative at x
        
        if overridden in a subclass, the signature should be
        derivative(self,x,dx=1) , but dx may be ignored
        """
        return (self(x+dx)-self(x))/dx
    
    
    #Must Override:
    def f(self,x):
        """
        this function MUST be overriden in subclasses - default arguments will
        be taken as initial parameter values, and unspecified defaults will
        default to 1
        
        the first parameter (other than 'self' or 'cls', if provided) must be the data
        vector/matrix (and have no default value), and the rest are parameters
        """
        raise NotImplementedError
        
        
class LinearModel(FunctionModel):
    """
    y=mx+b linear fit
    """
    def f(self,x,m=1,b=0):
        return m*x+b
    
    def _customFit(x,y,fixedpars=(),**kwargs):
        """
        does least-squares fit on the x,y data
        
        fixint and fixslope can be used to specify the intercept or slope of the fit
        or leave them free by having fixint or fixslope be False or None
        
        lastFit stores (m,b,dm,db,dy)
        """
        raise NotImplementedError('Linear not done yet')
        
        fixslope = 'm' in fixedpars
        fixint = 'b' in fixedpars
        
        N=len(x) 
        if fixint is False and fixslope is False:
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
            
        elif fixslope is not False:
            m,dm=fixslope,0
            sxsq=np.sum(x*x)
            b=np.sum(y-m*x)/N
            
            dy=(y-m*x-b).std(ddof=1)
            db=dy #TODO:check
            
        elif fixint is not False:
            b,db=fixslope,0
            sxy=np.sum(x*y)
            sxsq=np.sum(x*x)
            m=(sxy-N*b)/sxsq
            
            dy=(y-m*x-b).std(ddof=1)
            dm=dy*(sxsq/N)**0.5
        else:
            raise ValueError("can't fix both slope and intercept")
        
        self.lastFit=(m,b,dm,db,dy)
        
        return np.array((m,b))
    
    def derivative(self,x,dx=1):
        return np.ones_like(x)*m
    
    def integrate(self,lower,upper):
        return m*upper*upper/2+b*upper-m*lower*lower/2+b*lower
    
    def weightedFit(x,y,sigmay=None,doplot=False):
        """
        does a linear weighted least squares fit and computes the coefficients and errors
        fit is y=B*x+A
        if sigma is None, the weights are all equal - otherwise, it's the stddev of the y values
        returns B,A,sigmaB,sigmaA
        """
        raise NotImplementedError('needs to be adapted to astro.models')
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
    
class QuadraticModel(FunctionModel):
    """
    2-degree polynomial
    """
    def f(self,x,c2=1,c1=0,c0=0):
        return c2*x*x+c1*x+c0

class PolynomialModel(FunctionModel):
    """
    arbitrary-degree polynomial
    """
    #TODO: use polynomial objects that are only updated when necessary
    def f(self,x,*args): 
        return np.polyval(np.array(args)[::-1],x)
    
    def derivative(self,x):
        return np.polyder(np.array(self.parvals)[::-1])(x)

    def integrate(self,lower,upper):
        p = np.polyint(np.array(self.parvals)[::-1])
        return p(upper)-p(lower)

class GaussianModel(FunctionModel):
    """
    Normalized 1D gaussian function
    """
    def f(self,x,A=1,sig=1,mu=0):
        tsq=(x-mu)*2**-0.5/sig
        return A*np.exp(-tsq*tsq)*(2*pi)**-0.5/sig
    
    def derivative(self,x,dx=1):
        sig=self.sig
        return self(x)*-x/sig/sig
    
class LorentzianModel(FunctionModel):
    def f(self,x,A=1,gamma=1,mu=0):
        return gamma/pi/(x*x-2*x*mu+mu*mu+gamma*gamma)
    
class VoigtModel(GaussianModel,LorentzianModel):
    """
    Convolution of a Gaussian and Lorentzian profile
    """
    def f(self,x,A=1,sig=1,gamma=1,mu=0):
        from scipy.special import wofz
        w=wofz((x+1j*gamma)*2**-0.5/sig)
        return A*w.real*(2*pi)**-0.5/sig
    
class ExponentialModel(FunctionModel):
    """
    exponential function Ae^(kx)
    """
    def f(self,x,A=1,k=1):
        return A*np.exp(k*x)
    
class PowerLawModel(FunctionModel):
    """
    A single power law model Ax^p+B 
    """
    def f(self,x,A=1,p=1,B=0):
        return A*x**p+B
    
class TwoPowerModel(FunctionModel):
    """
    A model that smoothly transitions between two power laws at the turnover 
    point xs.  a is the inner slope, b is the outer slope
    """
    def f(self,x,A=1,xs=1,a=1,b=2):
        return A*((x+xs)**(b-a))*(x**a)
    
class SplineModel(FunctionModel):
    """
    this uses the scipy.interpolate.UnivariateSpline class as a model for the function
    """
    from scipy.interpolate import UnivariateSpline
    def __init__(self,xdata=None,ydata=None,degree=3):
        from warnings import warn
        if xdata is None and ydata is None:
            warn('Spline models should be provided with initial xdata and ydata (initializing to (0,0 and 1,1))')
            self._x,self._y = np.array((0,1)),np.array((0,1))
        elif xdata is None or ydata is None:
            raise ValueError('need both xdata and ydata')
        
        self._x,self._y = np.array(xdata),np.array(ydata)
        if xdata.shape != ydata.shape:
            raise ValueError("xdata and ydata don't match")
        self.s = xdata.size
        self._olds=None
        self.degree = degree
            
    def _fitCustom(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        self.spline = UnivariateSpline(self._x,self._y,s=self.s,k=self.degree)
        self._olds = self.s
        return np.array([s])
        
    def fitData(x,y,*args,**kwargs):
        self._x=x
        self._y=y
        self._olds=None
        return FunctionModel.fitData(x,y,*args,**kwargs)
    
    def f(self,x,s=2):
        if self._olds != s:
            self.spline = UnivariateSpline(self._x,self._y,s=s,k=self.degree)
            self._olds = s
        
        return self.spline(x)
    
class NFWModel(TwoPowerModel):
    """
    A Navarro, Frenk, and White 1996 profile
    """
    def f(self,x,rho0=1,rc=1): #TODO: see if its better to handcode the A->rho0
        return TwoPowerModel.f(self,x,rho0*rc*rc,rc,-1,-3)
    
    
    @staticmethod
    def Delta(z):
        """
        Virial overdensity - value is from Maller&Bullock 2004
        """
        return 360.0/(1.0+z)
    
    @staticmethod
    def Rvir_to_Mvir(Rvir,z=0,h=.72,Omega0=1):
        """
        M_sun,kpc
        """
        return 1e12/h (Omega0*NFWModel.Delta(z)/97.2)*(Rvir*(1+z)/203.4/h)**3
    
    @staticmethod
    def Mvir_to_Rvir(Mvir,z=0,h=.72,Omega0=1):
        """
        M_sun,kpc
        """
        return 203.4/h*(Omega0*NFWModel.Delta(z)/97.2)**(-1/3)*(Mvir/1e12/h)**(1/3)/(1+z)
    
    @staticmethod
    def Mvir_to_Vvir(Mvir,z=0,h=.72,Omega0=1):
        """
        km/s,M_sun
        """
        return 143.8*(Omega0*NFWModel.Delta(z)/97.2)**(1/6)*(Mvir/1e12/h)**(1/3)*(1+z)**0.5
    
    @staticmethod
    def Vvir_to_Mvir(Vvir,z=0,h=.72,Omega0=1):
        """
        km/s,M_sun
        """
        return (Omega0*NFWModel.Delta(z)/97.2)**-0.5*(1+z)**-1.5*h*1e12*(Vvir/143.8)**3

class PlummerModel(FunctionModel):
    def f(self,r,rp=1.,M=1.):
        return 3*M/(4.*pi*rp**3)*(1+(r/rp)**2)**-2.5

class King2DModel(FunctionModel):    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        return A*rcsq*((r*r+rcsq)**-0.5 - (rt*rt+rcsq)**-0.5)**2
    
class King3DModel(FunctionModel):
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        z=((r*r+rcsq)**0.5) * ((rt*rt+rcsq)**-0.5)
        return (A/z/z/pi/rc)*((1+rt*rt/rcsq)**-1.5)*(np.arccos(z)/z-(1-z*z)**0.5)

class SchecterModel(FunctionModel):
    def f(self,M,Mstar=-20.2,alpha=-1,phistar=1.0857362047581294):
        from numpy import log,exp
        x=10**(0.4*(Mstar-M))
        return 0.4*log10(10)*phistar*(x**(1+alpha))*exp(-x)
    
class SchecterModelLum(SchecterModel):
    def f(self,L,Lstar=1e10,alpha=-1.0,phistar=1.0):
        from astro.phot import lum_to_mag as l2m
        M,Mstar=l2m(L),l2m(Lstar)
        return SchecterModel.f(self,M,Mstar,alpha,phistar)
    #TODO:check to make sure this is actually the right way
        
        
class SersicModel(FunctionModel):
    #TODO:better normalizations and rs
    def f(self,r,A=1,rs=1,n=4):
        return A*np.exp(-(r/rs)**(1/n))
        
    #TODO:versions of exponential disk and deVauculeurs
