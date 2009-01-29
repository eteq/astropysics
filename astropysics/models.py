#Erik Tollerud (etolleru@uci.edu) 2008

"""
This module contains objects and functions for fitting data to models as well as
calculations and estimates from these models.

The aim of these classes are mostly for easily and quickly generating a range of
models - subclasses with just a function "f" will do all the expected things 
right out of the box.

Currently, the main fitting algorithms are those from scipy.optimize and the 
PyMC package (http://code.google.com/p/pymc/)
"""

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

class _FuncMeta1D(type):
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
            cls._args = None
            cls.__statargs = list(args)
        else:
            cls._args = tuple(args)
        
            
        cls.params = property(lambda cls:cls._args)
        def _set_parvals(cls,newpars):
            if len(newpars)>len(cls._args):
                raise ValueError('too many new parameters')
            for a,v in zip(cls._args,newpars):
                setattr(cls,a,v)        
        def _set_pardict(cls,newpardict):
            for k in newpardict.keys():
                if k not in cls._args:
                    raise KeyError('No parameter %s found'%k)
            for k,v in newpardict.iteritems():
                setattr(cls,k,v)
        cls.parvals = property(lambda cls:[getattr(cls,a) for a in cls._args],fset=_set_parvals)
        cls.pardict = property(lambda cls:dict([(a,getattr(cls,a)) for a in cls._args]),fset=_set_pardict)
        
    def __call__(cls,*args,**kwargs):
        if cls._args is None: #this is the case for *args in the function
            args=list(args)
            try:
                nparams=args.pop(0) if 'nparams' not in kwargs else kwargs.pop('nparams')
            except IndexError:
                IndexError('No # of parameters found for variable-size function')
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._args])    
            obj = type.__call__(cls,**objkwargs) #object __init__ is called here
            pars = cls.__statargs
            del cls.__statargs
            for i in range(nparams):
                p='p%i'%i
                pars.append(p)
                setattr(obj,p,1) #default varargs to 1
            cls._args = tuple(pars)
        else: #this is the case for fixed functions
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._args])
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

class FunctionModel1D(object):
    """
    This class is the base for 1-dimensional models with parameters that are 
    implemented as python functions.
    
    Subclassing:
    The following method MUST be overridden in a subclass:
    *f(self,x,...)
    The following methods may be overridden for speed of some operations - pars
    should be accessed with self.pardict, self.parvals, or by properties/name :
    *integrate(self,lower,upper)
    *derivative(self,x,dx)
    *_customFit(x,y,fixedpars=(),**kwargs)
    
    The metaclass generates the following properties:
    Parameters for each of the inputs of the function (except self and x)
    self.params: a tuple of the parameter names
    self.parvals: a list of the values in the parameters
    self.pardict: a dictionary of the parameter names and values
    """
    __metaclass__ = _FuncMeta1D
    
    defaultIntMethod = 'quad'
    
    def __call__(self,x):
        """
        call the function on the input x with the current parameters and return
        the result
        """
        return self.f(np.atleast_1d(x),*self.parvals)
    
    def optimize(self,x0,type='min',method='fmin',**kwargs):
        """
        find an optimal value for the model - x0 is where to start the search
        type can be 'min','max','root',or 'saddle'
        method can be 'fmin' or 'fmin_powell'
        """
        from scipy.optimize import fmin,fmin_powell
        
        if type == 'min':
            g=lambda *args,**kwargs:self.f(*args,**kwargs)
        elif type == 'max':
            g=lambda *args,**kwargs:-1*self.f(*args,**kwargs)
        elif type == 'root':
            g=lambda *args,**kwargs:np.abs(self.f(*args,**kwargs))
        elif type == 'saddle':
            raise NotImplementedError
        else:
            raise ValueError('Unrecognized type')
        
        if method == 'fmin':
            res=fmin(g,x0,tuple(self.parvals),**kwargs)
        elif method == 'fmin_powell':
            res=fmin_powell(g,x0,tuple(self.parvals),**kwargs)
        else:
            raise ValueError('Unrecognized method')
        
        self.lastOpt = res
        return res[0]
        
        
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
        
        returns the vector of the best fit parameters, and sets the parameters 
        if updatepars is True
        
        the full output is available is self.lastFit
        
        see also:getMCMC
        """
        from scipy import optimize as opt
        
        if method is None:
            #TODO: figure out a less hackish way than matching the docs
            method = 'leastsq' if self._customFit.__doc__ is FunctionModel1D._customFit.__doc__ else 'custom'
            
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        if timer:
            raise NotImplementedError
        
        ps=list(self.params)
        v=list(self.parvals) #initial guess
        
        
        if method == 'custom':
            res = self._customFit(x,y,fixedpars=fixedpars,**kwargs) 
            #ensure that res is at least a tuple with parameters in elem 0
            from operator import isSequenceType
            if not isSequenceType(res[0]):
                res = (res,)
        else:
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
                
            if method == 'leastsq':
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
    
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        Must be overridden to use 'custom' fit technique
        
        fixedpars will be a sequence of parameters that should be kept fixed
        for the fitting (possibly/typically empty)
        
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
        "bootstrap" technique to estimate errors - that is, use sampling
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
        
    def getMCMC(self,x,y,priors={},datamodel=None):
        """
        Note that this function requires PyMC (http://code.google.com/p/pymc/)
        
        To specify the priors, either supply pymc.Stochastric objects, a 2-tuple 
        (uniform lower and upper), a scalar > 0 (gaussian w/ the center given by 
        the current value and sigma provided by the scalar), or 0 (poisson with 
        k set by the current value)
        
        Any missing priors will raise a ValueError
        
        datamodel can be:
        *None: a normal distribution with sigma given by the data's standard 
        deviation is used as the data model
        *a tuple: (dist,dataname,kwargs)the first element is the 
        pymc.distribution to be used as the distribution representing the data
        and the second is the name of the argument to be associated with the 
        FunctionModel1D's output, and the third is kwargs for the distribution 
        ("observed" and "data" will be ignored, as will the data argument)
        *a scalar or sequence of length == data: a normal distribution is used 
        with sigma specified by the scalar/sequence
        
        returns a pymc.MCMMC object for Markov Chain Monte Carlo sampling
        """
        import pymc
        from operator import isSequenceType
        
        if set(priors.keys()) != set(self.params):
            raise ValueError("input priors don't match function params")
        
        d={}
        for k,v in priors.iteritems():
            if isinstance(v,pymc.StochasticBase):
                d[k]=v
            elif isSequenceType(v):
                if len(v) == 2:
                    d[k]=pymc.distributions.Uniform(k,v[0].v[1])
                else:
                    raise ValueError("couldn't understand sequence "+str(v) )
            elif v > 0:
                d[k] = pymc.distributions.Normal(k,self.pardict[k],1.0/v/v)
            elif v == 0:
                d[k] = pymc.distributions.Poisson(k,self.pardict[k])
            else:
                raise ValueError("couldn't interpret prior "+str(v))
        
        funcdict=dict(d)    
        funcdict['x']=x
        funcdet=pymc.Deterministic(name='f',eval=self.f,parents=funcdict,doc="FunctionModel1D function")
        d['f'] = funcdet
        
        if type(datamodel) is tuple:
            distobj,dataarg,kwargs=datamodel
            if 'name' not in kwargs:
                kwargs['name'] = 'data'
                
            kwargs[dataarg] = funcdet
            kwargs['observed'] = True
            kwargs['data'] = y
            
            datamodel = distobj(**kwargs)
        else:
            if datamodel is None:
                sig = np.std(y)
            else:
                sig = datamodel
            datamodel = pymc.distributions.Normal('data',mu=funcdet,tau=1/sig/sig,oberved=True,data=y)
        d[datamodel.__name__ ]=datamodel
        
        return pymc.MCMC(d)
        
    def plot(self,lower=None,upper=None,n=100,integrate=None,clf=True,logplot='',
              powerx=False,powery=False,deriv=None,data=None,fit = False,*args,**kwargs):
        """
        plot the model function from lower to upper with n samples
        
        integrate controls whether or not to plot the integral of the function -
        if True, the base of the integral is taken to be lower, if False, upper,
        no integration if None, and otherwise, the base is the argument value
        
        extra args and kwargs go into the matplotlib plot function
        
        data is either an x,y pair of data points or a dictionary of kwargs into 
        scatter fit determines if the model should be fit to the data before
        plotting.
        
        logplot determines whether or not to plot on a log scale, and powerx and
        poweru determine if the model points should be used as powers (base 10)
        before plotting
        """
        from matplotlib import pyplot as plt
        from operator import isMappingType
        
        if (lower is None or upper is None) and data is None:
            raise ValueError("Can't choose limits for plotting without specifying or providing data")
        if fit:
            from operator import isMappingType
            if not isMappingType(fit):
                fit = {}
            if data is None:
                raise ValueError("Need to supply data to perform fit before plot")
            self.fitData(*data,**fit)
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
        plt.xlim(lower,upper)
        
    #Can Override:
    def integrate(self,lower,upper,method=None,n=None,jac=None,**kwargs):
        """
        This will numerically integrate the model function from lower to upper
        using scipy.integrate techniques
        
        method can be specified (a function from scipy.integrate) or if None,
        self.defaultIntMethod will be used
        
        n is either the number of data values or the order for ordered methods,
        or if sequence,will be used as the integration x-values (and lower/upper
        are ignored)
        
        jac is the jacobian factor as a function f(x,*params)
        
        returns value
        most methods will also store their result to "lastIntegrate"
        
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
            res=itg.quad(f,lower,upper,args=ps,full_output=1,**kwargs)
            if len(res) == 4:
                v,e,d,m = res
                #from warnings import warn
                #warn('Integration message: %s'%m)
                print 'Integration message:',m
        #use these for 2d and 3d
        #elif method=='dblquad':
        #    raise NotImplementedError
        #elif method=='tplquad':
        #    raise NotImplementedError
        elif method=='fixed_quad':
            res=itg.fixed_quad(f,lower,upper,ps,5 if n is None else n,**kwargs)
        elif method=='quadrature':
            res=itg.quadrature(f,lower,upper,ps,**kwargs)
        elif method=='romberg':
            res=itg.romberg(f,lower,upper,ps,**kwargs)
        else: #sampled techniques
            if n is None:
                n=100
            if np.isscalar(n):
                x=np.linspace(lower,upper,n)
            else:
                x=np.array(n)
            y=f(x,*ps)
            if method=='trapz':
                res=itg.trapz(y,x,**kwargs)
            elif method=='cumtrapz':
                res=itg.cumtrapz(y,x,**kwargs)
            elif method=='simps':
                res=itg.simps(y,x,**kwargs)
            elif method=='romb':
                res=itg.simps(y,np.convolve(x,[1,-1],mode='same').mean(),**kwargs)
            else:
                raise ValueError('unrecognized integration method')
        
        
        self.lastIntegrate = res
        return res if np.isscalar(res) else res[0]
    
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
        
        the first parameter (other than 'self' or 'cls', if provided) must be 
        the data vector/matrix (and have no default value), and the rest are 
        parameters
        """
        raise NotImplementedError
        
        
class LinearModel(FunctionModel1D):
    """
    y=mx+b linear fit
    """
    def f(self,x,m=1,b=0):
        return m*x+b
    
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        does least-squares fit on the x,y data
        
        fixint and fixslope can be used to specify the intercept or slope of the 
        fit or leave them free by having fixint or fixslope be False or None
        
        lastFit stores ((m,b),dm,db,dy)
        """  
        
        fixslope = 'm' in fixedpars
        fixint = 'b' in fixedpars
        
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
        
        return (np.array((m,b)),dm,db,dy)
    
    def derivative(self,x,dx=1):
        return np.ones_like(x)*m
    
    def integrate(self,lower,upper):
        return m*upper*upper/2+b*upper-m*lower*lower/2+b*lower
    
    def weightedFit(self,x,y,sigmay=None,doplot=False):
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
  
class CompositeModel(FunctionModel1D):
    #TODO:initial vals
    def __init__(self,models=[],operation='+'):
        raise NotImplementedError
        self.op = operation
        ms=[]
        for m in models:
            if type(m) == str:
                m=get_model(m)
            if not issubclass(m,FunctionModel1D):
                raise ValueError('Non FunctionModel1D provided')
            ms.append(m)
        self._modelcounts=dict([(m,ms.count(m)) for m in set(ms)])
        self._models=[m() for m in ms]
        
        self.f=self.__f
        
        #TODO: set up attributes
        
    def _getParams(self):
        ps=[]
        from collections import defaultdict
        d=defaultdict(lambda:1)
        for m in self._models:
            i=d[m]
            d[m]=i+1
            mname=m.__name__.replace('Model','').replace('model','')
            for p in m.params:
                ps.append('%s_%i_%s'%(mname,i+1,p))
        return ps
    def _getParvals(self):
        ps=[]
        for m in self._models:
            for p in m.params:
                ps.append(getattr(m,p))
        return ps
    def _setParvals(self,val):
        i=0
        for m in self._models:
            for p in m.params:
                v=val[i]
                setattr(m,p,v)
                i+=1
    def _getPardict(self):
        from collections import defaultdict
        d=defaultdict(lambda:1)
        for m in self._models:
            i=d[m]
            d[m]=i+1
            mname=m.__name__.replace('Model','').replace('model','')
            ps=[]
            for p in m.params:
                ps.append(('%s_%i_%s'%(mname,i+1,p),getattr(m,p)))
        return dict(ps)
    def _setPardict(self,val):
        raise NotImplementedError
    
    params = property(_getParams)
    parvals = property(_getParvals,_setParvals)
    pardict = property(_getPardict,_setPardict)
        
    def f(self,x):
        raise RuntimeError('Placeholder function - this should never be reachable')
    def __f(self,x,*args):
        #TODO:optimize
        vals=[]
        i=0
        vals = [m(x) for m in self._models]
        if self.op == '+':
            return np.sum(vals)
        elif self.op == '*':
            return np.prod(vals)
        else:
            return ValueError('unrecognized operation')
            
        
        
#<----------------------Module functions -------------->  
def generate_composite_model(models,operation):
    return CompsiteModel(models=models,operation=operation)

__model_registry={}
def register_model(model,name=None,overwrite=False,stripmodel=True):
    if not issubclass(model,FunctionModel1D):
        raise TypeError('Supplied model is not a FunctionModel1D')
    if name is None:
        name = model.__name__.lower()
    else:
        name = name.lower()
    if stripmodel:
        name = name.replace('model','')
    if overwrite and name in __model_registry:
        raise KeyError('Model %s already exists'%name)
    
    __model_registry[name]=model
    
def get_model(modelname):
    return __model_registry[modelname]

def list_models():
    return __model_registry.keys()



#<----------------------Builtin models----------------->
class ConstantModel(FunctionModel1D):
    """
    the simplest model imaginable - just a constant value
    """
    def f(self,x,C=0):
        return C*np.ones_like(x)
    
    def derivative(self,x,dx=1):
        return np.zeros_like(x)
    
    def integrate(self,lower,upper):
        return self.C*(upper-lower)
    
class QuadraticModel(FunctionModel1D):
    """
    2-degree polynomial
    """
    def f(self,x,c2=1,c1=0,c0=0):
        return c2*x*x+c1*x+c0

class PolynomialModel(FunctionModel1D):
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

class GaussianModel(FunctionModel1D):
    """
    Normalized 1D gaussian function
    """
    def f(self,x,A=1,sig=1,mu=0):
        tsq=(x-mu)*2**-0.5/sig
        return A*np.exp(-tsq*tsq)*(2*pi)**-0.5/sig
    
    def _getPeak(self):
        return self(self.mu)[0]
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
        
    peak=property(_getPeak,_setPeak)
        
    def derivative(self,x,dx=1):
        sig=self.sig
        return self(x)*-x/sig/sig
    
class DoubleGaussianModel(FunctionModel1D):
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
            print dgm.pardict
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
            print dgm.pardict
            dgm.fitData(x,y,fixedpars=('mu2','B','sig2'),**kwargs)
        else:
            raise ValueError('unrecognized main component')
        print dgm.pardict
        dgm.fitData(x,y,fixedpars=(),**kwargs)
        print dgm.pardict
        return dgm
    
class LorentzianModel(FunctionModel1D):
    def f(self,x,A=1,gamma=1,mu=0):
        return gamma/pi/(x*x-2*x*mu+mu*mu+gamma*gamma)
    
    def _getPeak(self):
        return self(self.mu)[0]
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
        
    peak=property(_getPeak,_setPeak)
    
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
    
class ExponentialModel(FunctionModel1D):
    """
    exponential function Ae^(kx)
    """
    def f(self,x,A=1,k=1):
        return A*np.exp(k*x)
    
class PowerLawModel(FunctionModel1D):
    """
    A single power law model Ax^p+B 
    """
    def f(self,x,A=1,p=1,B=0):
        return A*x**p+B
    
class TwoPowerModel(FunctionModel1D):
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
    
class BlackbodyModel(FunctionModel1D):
    """
    a Planck blackbody radiation model
    
    inunit can be set to a variety of energy, wavelength, or frequency options
    """
    from .constants import h,c,kb
    
    def __init__(self,inunit='wl',outunit='erg'):
        self.inunit = inunit
        self.outunit = outunit
        self.stephanBoltzmannLaw = self._instanceSBLaw
    
    def f(self,x,A=1,T=5800):
        raise NotImplementedError
    
    def _flambda(self,x,A=1,T=5800):
        h,c,k=self.h,self.c,self.kb
        x=self._scaling*x
        return A*self._enscale*self._scaling*2*h*c*c*x**-5/(np.exp((h*c/(k*T*x)))-1)
    
    def _fnu(self,x,A=1,T=5800):
        h,c,k=self.h,self.c,self.kb
        x=self._scaling*x
        return A*self._enscale*self._scaling*2*h/c/c*x**3/(np.exp(h*x/(k*T))-1)
    
    def _fen(self,x,A=1,T=5800):
        x=self._scaling*x
        raise NotImplementedError
    
    def _getInunit(self):
        return self._inunit
    def _setInunit(self,inunit):
        u = inunit.lower()
        if inunit == 'wl' or inunit == 'lambda' or u == 'ang' or u == 'angstroms' or u == 'wavelength-angstrom':
            self._inunit = 'wavelength-angstrom'
            self.f = self._flambda
            self._scaling=1e-8
        elif u == 'nm' or inunit == 'wavelength-nm':
            self._inunit = 'wavelength-nm'
            self.f = self._flambda
            self._scaling=1e-7
        elif u == 'microns' or u == 'um' or u == 'wavelength-micron':
            self._inunit = 'wavelength-micron'
            self.f = self._flambda
            self._scaling=1e-4
        elif inunit == 'm' or u == 'wavelength-m':
            self._inunit = 'wavelength-m'
            self.f = self._flambda
            self._scaling=1e2
        elif inunit == 'cm' or u == 'wavelength-cm':
            self._inunit = 'wavelength-cm'
            self.f = self._flambda
            self._scaling=1
        elif inunit == 'f' or inunit == 'nu' or u == 'hz' or u == 'frequency-hz':
            self._inunit = 'frequency-Hz'
            self.f = self._fnu
            self._scaling=1
        elif u == 'thz' or inunit == 'frequency-THz':
            self._inunit = 'frequency-THz'
            self.f = self._fnu
            self._scaling=1e12
        elif u == 'e' or inunit == 'energy-eV':
            from .constants import ergperev
            self._inunit = 'energy-eV'
            self.f = self._fen
            self._scaling=ergperev
        elif u == 'erg' or inunit == 'energy-erg':
            self._inunit = 'energy-erg'
            self.f = self._fen    
            self._scaling=1
        elif inunit == 'J' or inunit == 'energy-J':
            self._inunit = 'energy-J'
            self.f = self._fen    
            self._scaling=1e-7
        else:
            raise ValueError('unrecognized inunit')
    inunit = property(_getInunit,_setInunit)
    
    def _getOutunit(self):
        return self._outunit
    def _setOutunit(self,outunit):
        outunit=' '.join(outunit.split()[:2])
        u=outunit.lower()
        if 'J' in outunit or 'joule' in u:
            self._outunit='J'
            self._enscale=1e-7
        elif 'ev' in u or 'electronvolt' in u:
            from .constants import ergperev
            self._outunit='eV'
            self._enscale=1/ergperev
        elif 'surface brightness' in u or 'mag' in u or 'sb' in u:
            raise NotImplementedError
            self._outunit='mag'
        else: # assume ergs as default for cgs
            self._outunit='erg'
            self._enscale=1
            
        if 'ang' in u:
            self._outunit+=' angstroms^-2'
            self._enscale*=1e-16
        elif 'm' in u.replace('mag','').replace('nm','').replace('angstrom','').replace('cm',''):
            self._outunit+=' m^-2'
            self._enscale*=1e4
        elif 'nm' in u:
            self._outunit+=' nm^-2'
            self._enscale*=1e-14
        else: #assume cm for cgs as default
            self._outunit+=' cm^-2'
            self._enscale*=1
            
        self._outunit+=' %s^-1 s^-1'%self._inunit.split('-')[1]
        
    outunit = property(_getOutunit,_setOutunit)
    
    def _getArea(self):
        return self.A/4/pi
    def _setArea(self,area):
        self.A=area*4*pi
    area = property(_getArea,_setArea)
    
    def setFlux(self,radius,distance):
        """
        sets A so that the output is the flux at the specified distance from
        a blackbody with the surface given by "area"
        """
        ratio=radius/distance
        self.A=pi*ratio*ratio
    def getFluxRadius(self,distance):
        """
        determines the radius of a spherical blackbody at the specified distance
        assuming the flux is given by the model at the given temperature
        """
        return (self.A*distance*distance/pi)**0.5 
    def getFluxdistance(self,radius):
        """
        determines the distance to a spherical blackbody of the specified radius
        assuming the flux is given by the model at the given temperature
        """
        return (pi*radius*radius/self.A)**0.5
    
    def _getPeak(self):
        h,k = self.h,self.kb
        if self.f == self._flambda:
            b = .28977685 #cm * K
            peakval = b/self.T/self._scaling
        elif self.f == self._fnu:
            a=2.821439 #constant from optimizing BB function
            peakval=a/h*k*self.T/self._scaling
        elif self.f == self._fen:
            raise NotImplementedError
        else:
            raise RuntimeError('Should never see this - bug in BB code')
        return self(peakval)[0]
    
    def _setPeak(self,peakval):
        self.A = 1
        self.A = peakval/self._getPeak()
    
    peak=property(_getPeak,_setPeak)
    
    def wienDisplacementLaw(self,peakval,updatepars=True):
        """
        uses the Wien Displacement Law to calculate the temperature
        """
        h,k = self.h,self.kb
        if self.f == self._flambda:
            b = .28977685 #cm * K
            T = b/peakval/self._scaling
        elif self.f == self._fnu:
            a=2.821439 #constant from optimizing BB function
            peakval=a/h*k*self.T/self._scaling
            T=peakval*self._scaling*h/a/k
        elif self.f == self._fen:
            raise NotImplementedError
        else:
            raise RuntimeError('Should never see this - bug in BB code')
        if updatepars:
            self.T=T
        return T
    
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
    
class BlackbodyOffsetModel(BlackbodyModel):
    """
    This is a Planck spectrum with y and x offsets
    """
    
    def f(self,x,A=1,T=5800,xoff=0,yoff=0):
        raise NotImplementedError
    
    def _flambda(self,x,A=1,T=5800,xoff=0,yoff=0):
        return BlackbodyModel._flambda(self,x+xoff,A,T)+yoff
    
    def _fnu(self,x,A=1,T=5800,xoff=0,yoff=0):
        return BlackbodyModel._fnu(self,x+xoff,A,T)+yoff
    
    def _fen(self,x,A=1,T=5800,xoff=0,yoff=0):
        return BlackbodyModel._fen(self,x+xoff,A,T)+yoff
    
    
class SplineModel(FunctionModel1D):
    """
    this uses the scipy.interpolate.UnivariateSpline class as a model for the function
    """
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
            
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        from scipy.interpolate import UnivariateSpline
        
        self.spline = UnivariateSpline(self._x,self._y,s=self.s,k=self.degree)
        self._olds = self.s
        return np.array([self.s])
        
    def fitData(self,x,y,*args,**kwargs):
        self._x=x
        self._y=y
        self._olds=None
        return FunctionModel1D.fitData(x,y,*args,**kwargs)
    
    def f(self,x,s=2):
        from scipy.interpolate import UnivariateSpline
        
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

class PlummerModel(FunctionModel1D):
    def f(self,r,rp=1.,M=1.):
        return 3*M/(4.*pi*rp**3)*(1+(r/rp)**2)**-2.5

class King2DModel(FunctionModel1D):    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        return A*rcsq*((r*r+rcsq)**-0.5 - (rt*rt+rcsq)**-0.5)**2
    
class King3DModel(FunctionModel1D):
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        z=((r*r+rcsq)**0.5) * ((rt*rt+rcsq)**-0.5)
        return (A/z/z/pi/rc)*((1+rt*rt/rcsq)**-1.5)*(np.arccos(z)/z-(1-z*z)**0.5)

class SchecterModel(FunctionModel1D):
    def f(self,M,Mstar=-20.2,alpha=-1,phistar=1.0857362047581294):
        from numpy import log,exp
        x=10**(0.4*(Mstar-M))
        return 0.4*np.log10(10)*phistar*(x**(1+alpha))*exp(-x)
    
class SchecterModelLum(SchecterModel):
    def f(self,L,Lstar=1e10,alpha=-1.0,phistar=1.0):
        from .phot import lum_to_mag as l2m
        M,Mstar=l2m(L),l2m(Lstar)
        return SchecterModel.f(self,M,Mstar,alpha,phistar)
    #TODO:check to make sure this is actually the right way
        
class EinastoModel(FunctionModel1D):
    def f(self,r,A=1,rs=1,alpha=.2):
        return A*np.exp(-(r/rs)**alpha)

class SersicModel(EinastoModel):
    #TODO:better normalizations and rs
    def f(self,r,A=1,rs=1,n=2):
        #return EinastoModel.f(self,r,A,rs,1/n)
        return A*np.exp(-(r/rs)**(1/n))
    
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
            raise ValuError('need data for limits or lower/upper')
        if data is not None:
            if upper is None:
                upper = np.max(data[0])
            if lower is None:
                lower = np.min(data[0])
        
        if clf:
            plt.clf()
        
        print lower,upper
        x = np.linspace(lower,upper,n)
        plt.plot(x,zpt-2.5*np.log10(self(x)))
        if data:
            skwargs={'c':'r'}
            plt.scatter(*data,**skwargs)
        
        plt.ylim(*reversed(plt.ylim()))
        plt.xlim(lower,upper)

class ExponentialDiskModel(SersicModel):
    def f(self,rz,A=1,rs=1,zs=.3):
        r=rz[0]
        z=rz[1]
        return SersicModel.f(self,r,A,rs,1)*np.exp(np.abs(z)/zs)
    
class DeVauculeursModel(SersicModel):
    def f(self,r,A=1,rs=1):
        return SersicModel.f(self,r,A,rs,4)

class MaxwellBoltzmannModel(FunctionModel1D):
    from .constants import me #electron
    def f(self,v,T=273,m=me):
        from .constants import kb,pi
        return (m/(2*pi*kb*T))**0.5*np.exp(-m*v*v/2/kb/T)
    
class MaxwellBoltzmanSpeedModel(MaxwellBoltzmannModel):
    from .constants import me #electron
    def f(self,v,T=273,m=me):
        from .constants import kb,pi
        return 4*pi*v*v*(m/(2*pi*kb*T))**1.5*np.exp(-m*v*v/2/kb/T)
    
#register all Models in this module
for o in locals().values():
    if type(o) == type(FunctionModel1D) and issubclass(o,FunctionModel1D) and o != FunctionModel1D and o!= CompositeModel:
        register_model(o)
