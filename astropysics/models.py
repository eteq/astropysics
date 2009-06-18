#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

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
try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
from .spec import HasSpecUnits as _HasSpecUnits


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
            raise ValueError('Parameters must be scalar values - type '+str(type(value))+' provided')
        setattr(obj,'_Parameter__'+self.name,value)
    
    def __delete__(self,obj):
        raise AttributeError("can't delete a parameter")

class _FuncMeta1D(ABCMeta):
    def __init__(cls,name,bases,dct):
        #called on import of astro.models
        from inspect import getargspec
        
        super(_FuncMeta1D,cls).__init__(name,bases,dct)
        
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
            
    def __call__(cls,*args,**kwargs):
        if cls._args is None: #this is the case for *args in the function
            args=list(args)
            try:
                nparams=args.pop(0) if 'nparams' not in kwargs else kwargs.pop('nparams')
            except IndexError:
                IndexError('No # of parameters found for variable-size function')
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._args])    
            obj = super(_FuncMeta1D,_FuncMeta1D).__call__(cls,**objkwargs) #object __init__ is called here
            pars = cls.__statargs
            del cls.__statargs
            for i in range(nparams):
                p='p%i'%i
                pars.append(p)
                setattr(obj,p,1) #default varargs to 1
            cls._args = tuple(pars)
        else: #this is the case for fixed functions
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._args])
            obj = super(_FuncMeta1D,_FuncMeta1D).__call__(cls,**objkwargs) #object __init__ is called here
            
        obj.f(np.array([0]),*obj.parvals) #try once to check for a working function
            
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
    *inv(yval,*args,**kwargs)
    *_customFit(x,y,fixedpars=(),weights=None,**kwargs)
    The following attributes may be set for additional information:
    xAxisName
    yAxisName
    
    The metaclass generates Parameters for each of the 
    inputs of the function (except self and x)
    
    The initializer's arguments specify initial values for the parameters,
    and any non-parameter kwargs will be passed into the __init__ method
    """
    __metaclass__ = _FuncMeta1D
    
    defaultIntMethod = 'quad'
    defaultInvMethod = 'brentq'
    xAxisName = None
    yAxisName = None
    
    def __call__(self,x):
        """
        call the function on the input x with the current parameters and return
        the result
        """
        arr = np.array(x,copy=False,dtype=float)
        res = self.f(np.atleast_1d(arr),*self.parvals)
        if len(arr.shape) > 0:
            return res
        else:
            return res[0]
        
    @property
    def params(self):
        """
        a tuple of the parameter names
        """
        return self._args
    
    def _getParvals(self):
        return [getattr(self,a) for a in self._args]
    def _setParvals(self,newpars):
        if len(newpars)>len(self._args):
            raise ValueError('too many new parameters')
        for a,v in zip(self._args,newpars):
            setattr(self,a,v)        
    parvals = property(_getParvals,_setParvals,doc='a list of the values in the parameters') 
    
    def _getPardict(self):
        return dict([(a,getattr(self,a)) for a in self._args])
    def _setPardict(self,newpardict):
        for k in newpardict.keys():
            if k not in self._args:
                raise KeyError('No parameter %s found'%k)
        for k,v in newpardict.iteritems():
            setattr(self,k,v)
    pardict = property(_getPardict,_setPardict,doc='a dictionary of the parameter names and values')
    
    def inv(self,yval,*args,**kwargs):
        """
        Find the x value matching the requested y-value.  
        (note that yval must be a scalar)
        
        typical usage is inv(yval,a,b,method='brentq')
        
        the kwarg 'method' can be any of the root finders from scipy.optimize 
        scalar solvers are recommended (default is brentq or newton).  'method' 
        can also be a function that should take f(g(x),*args,**kwargs) and 
        return where g(x) is 0
        
        *args and **kwargs are passed into the minimizer
        """    
        import scipy.optimize
        
        if len(args) == 0:
            args=[0]
            method = kwargs.pop('method','newton')
        elif len(args) == 1:
            method = kwargs.pop('method','newton')
        else:
            method = kwargs.pop('method','brentq')
            
        if isinstance(method,basestring):
            f = getattr(scipy.optimize,method)
        else:
            f = method
        
        if kwargs.pop('abs',False):
            g = lambda(x):np.abs(self(x)-yval)
        else:
            g = lambda(x):self(x)-yval
        
        return f(g,*args,**kwargs)
        
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
        
        
    def fitData(self,x,y,method=None,fixedpars=None,weights=None,contraction='sumsq',
                 updatepars=True,savedata=True,timer=None,**kwargs):
        """
        This will use the data to adjust the parameters to fit using any of a
        variety of methods from scipy.optimize
        
        method can be any of the optimize functions (except constrained or
        scalar optimizers) or 'custom'
        default is custom, or if not present, leastsq
        
        fixedpars is a sequence of parameter names to leave fixed.  If it is
        None, the fixed parameters are inferred from self.fixedpars
        
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
        
        if savedata is true, the input fit data will be available in 
        self.fitteddata as an (x,y) tuple
        
        see also:getMCMC
        """
        from scipy import optimize as opt
        
        if method is None:
            #TODO: figure out a less hackish way than matching the docs
            method = 'leastsq' if self._customFit.__doc__ is FunctionModel1D._customFit.__doc__ else 'custom'
            
        if fixedpars is None:
            fixedpars = self.fixedpars if hasattr(self,'fixedpars') else ()
            
        if x.shape != y.shape:
            raise ValueError("x and y don't match")
        
        if timer:
            raise NotImplementedError
        
        ps=list(self.params)
        v=list(self.parvals) #initial guess
        
        w = 1 if weights is None else weights
        
        if method == 'custom':
            res = self._customFit(x,y,fixedpars=fixedpars,weights=weights,**kwargs) 
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
                    #return self.f(x,**pdict)
                    params = [pdict[a] for a in self._args]
                    return self.f(x,*params)
            else:
                f=lambda x,v:self.f(x,*v)
                
            if method == 'leastsq':
                if 'frac' in contraction:
                    g=lambda v,x,y:w*(1-f(x,v)/y)
                else:
                    g=lambda v,x,y:w*(y-f(x,v))
                res=opt.leastsq(g,v,(x,y),full_output=1,**kwargs)
            else:
                if not contraction:
                    g=lambda v,x,y:w*(y-f(x,v))
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
                        g=lambda v,x,y:np.sum(w*g1(v,x,y))
                    elif 'mean' in contraction:
                        g=lambda v,x,y:np.mean(w*g1(v,x,y))
                    elif 'median' in contraction:
                        g=lambda v,x,y:np.median(w*g1(v,x,y))
                    
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
                
        if savedata: 
            self.fitteddata = (x,y)
        
        return v
    
    def _customFit(self,x,y,fixedpars=(),weights=None,**kwargs):
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
              powerx=False,powery=False,deriv=None,data='auto',fit = False,*args,**kwargs):
        """
        plot the model function from lower to upper with n samples
        
        integrate controls whether or not to plot the integral of the function -
        if True, the base of the integral is taken to be lower, if False, upper,
        no integration if None, and otherwise, the base is the argument value.
        It can also be 's', 'spherical','c',or 'circular' for the 
        appropriate types of integration
        
        extra args and kwargs go into the matplotlib plot function
        
        data is either an x,y pair of data points or a dictionary of kwargs into 
        scatter fit determines if the model should be fit to the data before
        plotting. If it is 'auto', the last data input to fitData will be used
        (if savedata was true).  If it evaluates to False, no data will be 
        plotted and lower and upper must be set
        
        logplot determines whether or not to plot on a log scale, and powerx and
        poweru determine if the model points should be used as powers (base 10)
        before plotting
        """
        from matplotlib import pyplot as plt
        from operator import isMappingType
        
        if fit:
            from operator import isMappingType
            if not isMappingType(fit):
                fit = {}
            if data is None:
                raise ValueError("Need to supply data to perform fit before plot")
            self.fitData(*data,**fit)
            
        if data is 'auto':
            if hasattr(self,'fitteddata') and self.fitteddata:
                data = self.fitteddata
            else:
                data = None
        if data is not None:
            if lower is None:
                lower=np.min(data[0])
            if upper is None:
                upper=np.max(data[0])
                
        if (lower is None or upper is None) and data is None:
            raise ValueError("Can't choose limits for plotting without specifying or providing data")
        
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
            intfunc = self.integrate
            if integrate is True:
                base = lower
            elif integrate is False:
                base = upper
            elif isinstance(integrate,basestring):
                if integrate.startswith('s'): #spherical
                    intfunc = self.integrateSpherical
                    base = 0
                elif integrate.startswith('c'): #circular
                    intfunc = self.integrateCircular
                    base = 0
                else:
                    raise ValueError('unrecognized integrate string')
            else:
                base = float(integrate)
            y=np.array([intfunc(base,v) for v in x])
            
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
        
        if np.any(data):
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
        
        if self.xAxisName:
            plt.xlabel(self.xAxisName)
        if self.yAxisName:
            plt.ylabel(self.yAxisName)
    
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
    def integrateCircular(self,*args,**kwargs):
        """
        This is an alias for self.integrate with jacobian set for a azimuthally
        symmetric 2D radial profile
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*2.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*2.0*pi
        return self.integrate(*args,**kwargs)
    
    def integrateSpherical(self,*args,**kwargs):
        """
        This is an alias for self.integrate with jacobian set for a spherically
        symmetric 3D radial profile
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*x*4.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*x*4.0*pi
        return self.integrate(*args,**kwargs)
        
    def derivative(self,x,dx=1):
        """
        the derivative at x
        
        if overridden in a subclass, the signature should be
        derivative(self,x,dx=1) , but dx may be ignored
        """
        return (self(x+dx)-self(x))/dx
    
    @abstractmethod
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
  
  
class _CompMeta1D(_FuncMeta1D):
#    def __init__(cls,name,bases,dct):
#        super(_CompMeta1D,cls).__init__(name,bases,dct)
    def __call__(cls,*args,**kwargs):
        obj = super(_FuncMeta1D,_CompMeta1D).__call__(cls,*args,**kwargs) #object __init__ is called here
        return obj
        
class CompositeModel1D(FunctionModel1D):
    """
    This model contains a group of FunctionModel1D objects and evaluates them
    as a single model.
    
    The models can either be FunctionModel1D objects, FunctionModel1D classes,
    or a string (in the later two cases, new instances will be generated)
    
    parameter names are of the form 'A-0' and 'A-1' where A is the parameter
    name and the number is the sequential number of the model with that
    parameter
    
    the operations can be any valid python operator, or a sequence of operators
    to apply (e.g. ['+','*','+'] will do mod1+mod2*mod3+mod4
    
    any extra kwargs will be used to specify default values of the parameters
    """
    __metaclass__ = _CompMeta1D
    
    #TODO:initial vals
    def __init__(self,models=[],operation='+',**kwargs):
        from inspect import isclass
        
        self.__dict__['_args'] = tuple() #necessary for later when getattr is invoked
        
        mods = []

        for m in models:
            if isinstance(m,FunctionModel1D):
                pass
            elif isinstance(m,basestring):
                m = get_model(m)()
            elif isclass(m) and issubclass(m,FunctionModel1D):
                m = m()
            else:
                raise ValueError('Supplied object is not a function model')
            mods.append(m)
            
        self._models = tuple(mods)
        
        if isinstance(operation,basestring):
            operation = [operation for i in range(len(mods)-1)]
        elif len(operation) != len(mods)-1:
            raise ValueError('impossible number of operations')
        
        self._ops = tuple(operation)
        oplst = ['mval[%i]%s'%t for t in enumerate(self._ops)]
        oplst += 'mval[%i]'%len(oplst)
        self._opstr = ''.join(oplst)
        
        args = []
        argmodi = []
        for i,m in enumerate(self._models):
            for p in m.params:
                args.append(p+str(i))
                argmodi.append(i)
        self._args = tuple(args)
        self._argmodi = tuple(argmodi)
        
        self._filters = None
        
        for k,v in kwargs:
            setattr(self,k,v)
    
    def __getattr__(self,name):
        if name in self._args:
            i = self._args.index(name)
            j = self._argmodi[i]
            return getattr(self._models[j],name[:-len(str(j))])
        raise AttributeError("'%s' has no attribute '%s'"%(self,name))
    
    def __setattr__(self,name,val):
        if name in self._args:
            i = self._args.index(name)
            j = self._argmodi[i]
            setattr(self._models[j],name[:-len(str(j))],val)
        else:
            self.__dict__[name] = val
    
    def f(self,x,*args):
        #TODO: find out if the commented part can be sometimes skipped somehow
        for p,a in zip(self.params,args):
            setattr(self,p,a)
        mval = [m.f(x,*m.parvals) for m in self._models]
        res = eval(self._opstr)
        if self._filters is None:
            return res
        else:
            for filter in self._filters:
                res = filter(res)
            return res
        
    def addFilter(self,filter):
        """
        This adds a function to be applied after the model is evaluated
        """
        if self._filters is None:
            self._filters = []
        
        if not callable(filter):
            raise ValueError('input filter is not a function')
        
        self._filters.append(filter)
        
    def clearFilters(self):
        """
        this clears all previously added filters
        """
        self._filters = None
        
    def addLowerBoundFilter(self,bound):
        def bndfunc(x):
            x[x<bound] = bound
            return x
        self.addFilter(bndfunc)
        
    def fitDataFixed(self,*args,**kwargs):
        """
        Calls fitData with kwargs and args, but ignores fixedpars 
        argument and uses 'fixedmods' or 'freemods' kwarg to 
        determine which parameters should be fixed
        
        freemods or fixedmods should be a sequence of indecies of 
        models for which the parameters should be left free or held
        fixed
        """
        fps = []
        if 'fixedmods' in kwargs and 'freemods' in kwargs:
            raise TypeError('fitDataFixed cannot have both fixedmod and freemod arguments')
        elif 'fixedmods' in kwargs:
            for i in kwargs.pop('fixedmods'):
                stri = str(i)
                for p in self._models[i].params:
                    fps.append(p+stri)
        elif 'freemods' in kwargs:
            fps.extend(self.params)
            for i in kwargs.pop('freemods'):
                stri=str(i)
                for p in self._models[i].params:
                    fps.remove(p+stri)
                
        else:
            raise TypeError('fitDataFixed must have fixedmods or freemods as arguments')
            
        
        if len(args)>=4:
            args[3] = fps
        else:
            kwargs['fixedpars'] = fps
        print fps
        return self.fitData(*args,**kwargs)
            
    @property
    def models(self):
        return self._models
    
    @property
    def ops(self):
        return self._ops  
        
#<----------------------Module functions -------------->  
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
    
    def integrate(self,lower,upper,**kwargs):
        if 'jac' in kwargs and kwargs['jac'] is not None:
            return FunctionModel1D.integrate(self,lower,upper,**kwargs)
        else:
            return self.C*(upper-lower)
    
class LinearModel(FunctionModel1D):
    """
    y=mx+b linear fit
    """
    
    def f(self,x,m=1,b=0):
        return m*x+b
    
    def _customFit(self,x,y,fixedpars=(),weights=None,**kwargs):
        """
        does least-squares fit on the x,y data
        
        fixint and fixslope can be used to specify the intercept or slope of the 
        fit or leave them free by having fixint or fixslope be False or None
        
        lastFit stores ((m,b),dm,db,dy)
        """  
        if weights is not None:
            if fixedpars or len(kwargs)>0:
                raise NotImplementedError("can't fix pars with weighted linear fit yes")
            m,b,dm,db = self.weightedFit(x,y,1/weights,False)
            dy = (y-m*x-b).std(ddof=1)
            return (np.array((m,b)),dm,db,dy)
        
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
        return self(self.mu)
    
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
        return self(self.mu)
    
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
    
class TwoSlopeModel(FunctionModel1D):
    """
    This model smoothly transitions from linear with one slope to linear with
    a different slope. It is the linearized equivalent of TwoPowerModel.
    
    specifically, a*x+(b-a)*log(1+exp(x))+c
    """
    def f(self,x,a=1,b=2,c=0,xs=0):
        z = x-xs
        return a*z+(b-a)*np.log(1+np.exp(z))+c
    
class BlackbodyModel(FunctionModel1D,_HasSpecUnits):
    """
    a Planck blackbody radiation model.  

    y-axis is assumed to be specific intensity
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
    
    def wienDisplacementLaw(self,peakval,updatepars=True):
        """
        uses the Wien Displacement Law to calculate the temperature
        """
        h,k = self.h,self.kb
        if self.f == self._flambda:
            b = .28977685 #cm * K
            T = b/peakval/self._xscaling
        elif self.f == self._fnu:
            a=2.821439 #constant from optimizing BB function
            peakval=a/h*k*self.T/self._xscaling
            T=peakval*self._xscaling*h/a/k
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
    
   
class SmoothSplineModel(FunctionModel1D):
    """
    this uses a B-spline as a model for the function.  Note that by
    default the parameters are not tuned - the input smoothing and 
    degree are left alone when fitting
    
    the scipy.interpolate.UnivariateSpline class is used to
    do the calculation (in the "spline" attribute) 
    """
    def __init__(self):
        super(SmoothSplineModel,self).__init__()
        
        self._oldd=self._olds=self._ws=None
        self.fitteddata=(np.arange(self.degree+1),np.arange(self.degree+1))
            
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        from scipy.interpolate import UnivariateSpline
        
        self.spline = UnivariateSpline(x,y,s=self.s,k=self.degree,w=kwargs['weights'] if 'weights' in kwargs else None)
        
        self._olds = self.s
        self._oldd = self.degree
        
        return np.array([self.s,self.degree])
        
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
        return super(SmoothSplineModel,self).fitData(x[sorti],y[sorti],**kwargs)
    
    def f(self,x,s=2,degree=3):        
        if self._olds != s or self._oldd != degree:
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)
    
    
class KnotSplineModel(FunctionModel1D):
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
    def __init__(self,locmethod='even'):
        super(KnotSplineModel,self).__init__()
        
        self._oldd=self._oldk=self._ws=None
        self.locmethod = locmethod
        self.fitteddata=(np.arange(self.degree+1),np.arange(self.degree+1))
            
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        from scipy.interpolate import LSQUnivariateSpline
        if self.locmethod == 'evencdf':
            cdf,xcdf = np.histogram(x,bins=max(10,max(2*self.nknots,int(len(x)/10))))
            mask = cdf!=0
            cdf,xcdf = cdf[mask],xcdf[np.hstack((True,mask))]
            cdf = np.hstack((0,np.cumsum(cdf)/np.sum(cdf)))
            self.iknots = np.interp(np.linspace(0,1,self.nknots+2)[1:-1],cdf,xcdf)
        elif self.locmethod == 'even':
            self.iknots = np.linspace(x[0],x[-1],self.nknots+2)[1:-1]
        else:
            raise ValueError('unrecognized locmethod %s'%self.locmethod)
        self.spline = LSQUnivariateSpline(x,y,t=self.iknots,k=self.degree,w=kwargs['weights'] if 'weights' in kwargs else None)
        
        self._oldk = self.nknots
        self._oldd = self.degree
        
        return np.array([self.nknots,self.degree])
        
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
        return super(KnotSplineModel,self).fitData(x[sorti],y[sorti],**kwargs)
    
    def f(self,x,nknots=3,degree=3):        
        if self._oldk != nknots or self._oldd != degree:
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)
    
    
    
class NFWModel(TwoPowerModel):
    """
    A Navarro, Frenk, and White 1996 profile
    """
    def __init__(self):
        TwoPowerModel.__init__(self)
        self.a = -1
        self.b = -3
        
    def f(self,x,rho0=1,rc=1): #TODO: see if its better to handcode the A->rho0
        return TwoPowerModel.f(self,x,rho0*rc*rc*rc,rc,self.a,self.b)
    
    def integrateSpherical(self,*args,**kwargs):
        """
        NFW Has an analytic form for the spherical integral - if the inner 
        is not 0 or the form is anything other than integrateSpherical(0,r) 
        fall back to FunctionModel1D.integrateSpherical (or if the keyword
        'numerical' is true)
        """
        
        if kwargs.pop('numerical',False) or len(kwargs)>0 or len(args)!=2 or args[0]!=0.0:
            return FunctionModel1D.integrateSpherical(self,*args,**kwargs)
        else:
            x=args[1]/self.rc
            return 4*pi*self.rho0*self.rc**3*(np.log(1+x)-x/(1+x))
    def setRc(self,Rs,c,M):
        """
        sets the model parameters according to the given Rs/rc, concentration
        and mass enclosed
        """
        self.rc = Rs
        Router = c*Rs
        
        self.rho0 = 1
        a0 = self.integrateSpherical(0,Router)
        self.rho0 = M/a0
    
    def getRv(self,z=0):
        """
        get the virial radius at a given redshift (uses NFWModel.Delta(z))
        
        WARNING: may not be working right unit-wise
        """
        from scipy.optimize import newton
        
        try:
            from .constants import get_cosmology
            rhoC = get_cosmology().rhoC()
        except:
            raise ValueError('current cosmology does not support critical density')
        
        rhov = self.Delta(z)*rhoC
        return self.inv(rhov,1)
        
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
        return 0.4*log(10)*phistar*(x**(1+alpha))*exp(-x)
    
class SchecterModelLum(SchecterModel):
    def f(self,L,Lstar=1e10,alpha=-1.0,phistar=1.0):
        #from .phot import lum_to_mag as l2m
        #M,Mstar=l2m(L),l2m(Lstar)
        #return SchecterModel.f(self,M,Mstar,alpha,phistar)
        x = L/Lstar
        return phistar*(x**alpha)*np.exp(-x)/Lstar
    #TODO:check to make sure this is actually the right way
        
class EinastoModel(FunctionModel1D):
    def f(self,r,A=1,rs=1,alpha=.2):
        return A*np.exp(-(r/rs)**alpha)

class SersicModel(FunctionModel1D):
    """
    Sersic surface brightness profile:
    A*exp(-b_n[(R/Re)^(1/n)-1])
    """
    def f(self,r,A=1,re=1,n=2):
        #return EinastoModel.f(self,r,A,rs,1/n)
        #return A*np.exp(-(r/rs)**(1/n))
        return A*np.exp(-self.bn(n)*((r/re)**(1.0/n)-1))
    
    _bncache={}
    _bnpoly1=np.poly1d([-2194697/30690717750,131/1148175,46/25515,4/405,-1/3])
    _bnpoly2=np.poly1d([13.43,-19.67,10.95,-0.8902,0.01945])
    def bn(self,n,usecache=True):
        """
        bn is used to get the appropriate half-light radius
        
        the form is a fit from MacArthur, Courteau, and Holtzman 2003 
        and is claimed to be good to ~O(10^-5)
        
        if usecache is True, the cache will be searched, if False it will
        be saved but not used, if None, ignored
        """
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
        
        print lower,upper
        x = np.linspace(lower,upper,n)
        plt.plot(x,zpt-2.5*np.log10(self(x)))
        if data:
            skwargs={'c':'r'}
            plt.scatter(*data,**skwargs)
        
        plt.ylim(*reversed(plt.ylim()))
        plt.xlim(lower,upper)

class ExponentialDiskModel(SersicModel):
    def f(self,rz,A=1,re=1,zs=.3):
        r=rz[0]
        z=rz[1]
        return SersicModel.f(self,r,A,re,1)*np.exp(np.abs(z)/zs)
    
class DeVaucouleursModel(SersicModel):
    def f(self,r,A=1,re=1):
        return SersicModel.f(self,r,A,re,4)

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
    if type(o) == type(FunctionModel1D) and issubclass(o,FunctionModel1D) and o != FunctionModel1D and o!= CompositeModel1D:
        register_model(o)

#<---------------------------Other modelling techniques------------------------>

class Pca(object):
    """
    A basic class for Principal Component Analysis (PCA).
    
    p is the number of dimensions, while N is the number of data points
    """
    _colors=('r','g','b','c','y','m','k') #defaults
    
    def __calc(self):
        A = self.A
        M=A-np.mean(A,axis=0)
        N=M/np.std(M,axis=0)
        
        self.M = M
        self.N = N
        self._eig = None
    
    def __init__(self,data,names=None):
        """
        p X N matrix input
        """
        from warnings import warn
        A = np.array(data).T
        n,p = A.shape
        self.n,self.p = n,p
        if p > n:
            warn('p > n - intentional?')
        self.A = A
        self._origA=A.copy()
            
        self.__calc()
        
        self._colors= np.tile(self._colors,int((p-1)/len(self._colors))+1)[:p]
        if names is not None and len(names) != p:
            raise ValueError('names must match data dimension')
        self.names = None if names is None else tuple([str(n) for n in names])
        
        
    def cov(self):
        return np.cov(self.N.T)
        
    def eig(self):
        if self._eig is None:
            res = np.linalg.eig(self.cov())
            sorti=np.argsort(res[0])[::-1]
            res=(res[0][sorti],res[1][:,sorti])
            self._eig=res
        return self._eig
    def evals(self):
        return self.eig()[0]
    def evecs(self):
        return self.eig()[1]
    
    def energies(self):
        v=self.evals()
        return v/np.sum(v)
        
    def plot2d(self,ix=0,iy=1,clf=True):
        import matplotlib.pyplot as plt
        x,y=self.N[:,ix],self.N[:,iy]
        if clf:
            plt.clf()
        plt.scatter(x,y)
        vals,evs=self.eig()
        #evx,evy=evs[:,ix],evs[:,iy]
        xl,xu=plt.xlim()
        yl,yu=plt.ylim()
        dx,dy=(xu-xl),(yu-yl)
        for val,vec,c in zip(vals,evs.T,self._colors):
            plt.arrow(0,0,val*vec[ix],val*vec[iy],head_width=0.05*(dx*dy/4)**0.5,fc=c,ec=c)
        #plt.arrow(0,0,vals[ix]*evs[ix,ix],vals[ix]*evs[iy,ix],head_width=0.05*(dx*dy/4)**0.5,fc='g',ec='g')
        #plt.arrow(0,0,vals[iy]*evs[ix,iy],vals[iy]*evs[iy,iy],head_width=0.05*(dx*dy/4)**0.5,fc='r',ec='r')
        if self.names is not None:
            plt.xlabel('$'+self.names[ix]+'/\\sigma$')
            plt.ylabel('$'+self.names[iy]+'/\\sigma$')
        
    def plot3d(self,ix=0,iy=1,iz=2,clf=True):
        import enthought.mayavi.mlab as M
        if clf:
            M.clf()
        z3=np.zeros(3)
        v=(self.evecs()*self.evals())
        M.quiver3d(z3,z3,z3,v[ix],v[iy],v[iz],scale_factor=5)
        M.points3d(self.N[:,ix],self.N[:,iy],self.N[:,iz],scale_factor=0.3)
        if self.names:
            M.axes(xlabel=self.names[ix]+'/sigma',ylabel=self.names[iy]+'/sigma',zlabel=self.names[iz]+'/sigma')
        else:
            M.axes()
        
    def sigclip(self,sigs):
        if np.isscalar(sigs):
            sigs=sigs*np.ones(self.N.shape[1])
        n = self.N.shape[0]
        m = np.all(np.abs(self.N) < sigs,axis=1)
        self.A=self.A[m]
        self.__calc()
        return n-sum(m)
        
    def reset(self):
        self.A = self._origA.copy()
        self.__calc()
        
        
    def project(self,enthresh=None,nPCs=None,cumen=None,vals=None):
        """
        projects the normalized values onto the components
        
        enthresh, nPCs, and cumen determine how many PCs to use
        
        if vals is None, the normalized data vectors are the values to project
        
        returns n,p(>threshold) dimension array
        """
        nonnones = sum([e != None for e in (enthresh,nPCs,cumen)])
        if nonnones == 0:
            m = slice(None)
        elif nonnones > 1:
            raise ValueError("can't specify more than one threshold")
        else:
            if enthresh is not None:
                m = self.energies() > enthresh
            elif nPCs is not None:
                m = slice(None,nPCs)
            elif cumen is not None:
                m = np.cumsum(self.energies()) <  cumen
            else:
                raise RuntimeError('Should be unreachable')
        
        if vals is None:
            vals = self.N.T
        else:
            if self.N.T.shape[0] != vals.shape[0]:
                raise ValueError("shape for vals doesn't match")
        proj = np.matrix(self.evecs()).T*vals
        return proj[m].T
            
    def deproject(self,A,normed=True):
        """
        input is an n X q array, where q <= p
        
        output is p X n
        """
        A=np.atleast_2d(A)
        n,q = A.shape
        p = self.A.shape[1]
        if q > p :
            raise ValueError("q > p")
        
        evinv=np.linalg.inv(np.matrix(self.evecs()).T)
        
        zs = np.zeros((n,p))
        zs[:,:q]=A
        
        proj = evinv*zs.T
        
        if normed:
            return np.array(proj.T).T
        else:
            mns=np.mean(self.A,axis=0)
            sds=np.std(self.M,axis=0)
            return (np.array(proj.T)*sds+mns).T
    
    def subtractPC(self,pc,vals=None):
        """
        pc can be a scalar or any sequence of pc indecies
        
        if vals is None, the source data is self.A, else whatever is in vals
        (which must be p x m)
        """
        if vals is None:
            vals = self.A
        else:
            vals = vals.T
            if vals.shape[1]!= self.A.shape[1]:
                raise ValueError("vals don't have the correct number of components")
        
        pcs=self.project()
        zpcs=np.zeros_like(pcs)
        zpcs[:,pc]=pcs[:,pc]
        upc=self.deproject(zpcs,False)
        
        A = vals.T-upc
        B = A.T*np.std(self.M,axis=0)
        return B+np.mean(self.A,axis=0)






del ABCMeta,abstractmethod,abstractproperty #clean up namespace