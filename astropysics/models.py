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

from __future__ import division,with_statement
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

class DimensionError(Exception):
    """
    This exception indicates a problem with the dimensionality of 
    an input or model
    """
    def __init__(self,message):
        super(DimensionError,self).__init__(message)

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
    
class Model(object):
    """
    The superclass of all models
    
    subclasses should implement all abstract properties and methods:
    *__call__: determine the value at the location - note that
    dimension size-checking should not be performed here - the checkDims 
    method is for that purpose
    *params: a sequence of names for the parameters of the model
    *parvals: a sequence of values for the parameters of the model
    
    subclasses should also have the following attribute
    *ndims: a tuple of the form (indims,outdims)
    
    optional overrides:
    *pardict: a dictionary with keys as parameter names and values as 
    the value for that parameter
    *inv: compute the inverse of the function 
    """
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def __call__(self,x):
        raise NotImplementedError
        
    def checkDims(self,x,strict=False):
        """
        calls the function and checks the dimension - raises a DimensionError
        if the shapes are not appropriate
        """
        indim = len(np.array(x,copy=False).shape)
        if indim > self.indims:
            raise DimensionError('input dimension is too large')
        elif strict and indim < self.indims:
            raise DimensionError('input dimension is too small')
        
        y = self(x)
        
        outdim = len(np.array(y,copy=False).shape)
        if outdim > self.outdims:
            raise DimensionError('output dimension is too large')
        elif strict and outdim < self.outdims:
            raise DimensionError('output dimension is too small')
        
        return y
    
    
    @property
    def indims(self):
        return self.ndims[0]     
    @property
    def outdims(self):
        return self.ndims[1]
    
    params = abstractproperty(doc='a sequence of the parameter names')     
    parvals = abstractproperty(doc='a sequence of the values in the parameters')       
    @property
    def pardict(self):
        """
        a dictionary of the parameter names and values
        """
        return dict([t for t in zip(self.params,self.parvals)])      
    
    def inv(output,*args,**kwargs):
        """
        Compute the inverse of this model for the requested output
        
        It should return None if the model is not invertible, or raise a
        DimensionError if indims != outdims
        """
        if self.indims != self.outdims:
            raise DimensionError('inverse only meaningful if input and output dimensions are the same')
        return None
    
    

class _FuncMeta(ABCMeta):
    def __init__(cls,name,bases,dct):
        #called on import of astro.models
        from inspect import getargspec
        
        super(_FuncMeta,cls).__init__(name,bases,dct)
        if 'f' in dct:
            ffunc = dct['f']
        else:
            for b in bases:
                if 'f' in b.__dict__:
                    ffunc = b.__dict__['f']
                    break
            else:
                raise KeyError('missing f function')
        args,varargs,varkw,defaults=getargspec(ffunc)
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
            cls._statargs = list(args)
        else:
            cls._args = tuple(args)
            
    def __call__(cls,*args,**kwargs):
        if cls._args is None: #this is the case for *args in the function
            args=list(args)
            parsname = getattr(cls,'parsname') if hasattr(cls,'parsname') else 'p'
            
            try:
                nparams=args.pop(0) if 'nparams' not in kwargs else kwargs.pop('nparams')
            except IndexError:
                raise IndexError('No # of parameters found for variable-size function')
            objkwargs = dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._statargs if not k.startswith(parsname)])
            cls._currnparams = nparams #this is set in case the constructor needs to know for some reason
            obj = super(_FuncMeta,_FuncMeta).__call__(cls,**objkwargs) #object __init__ is called here
            del cls._currnparams
            pars = list(cls._statargs)
            
            for i in range(nparams):
                p = parsname+str(i)
                pars.append(p)
                if not hasattr(cls,p):
                    setattr(obj,p,1) #default varargs to 1
                else: #let class variables specify defaults for varargs
                    setattr(obj,p,getattr(cls,p))
            obj._args = tuple(pars)
        else: #this is the case for fixed functions
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._args])
            obj = super(_FuncMeta,_FuncMeta).__call__(cls,**objkwargs) #object __init__ is called here
            
        #TODO: smarter function check?
        #obj.f(np.zeros([1 for i in range(cls._ndims[0])]),*obj.parvals) #try once to check for a working function
            
        #set initial values
        if len(args)+len(kwargs) > len(obj.params):
            raise TypeError('Too many args and kwargs for the parameters')
        
        pars=list(obj.params)
        for k,v in kwargs.iteritems():
            del pars[pars.index(k)]
            setattr(obj,k,v)
        for p,v in zip(pars,args):
            setattr(obj,p,v)
            
        #bind _callf to the standard call if it isn't messed with by the function constructor
        if not hasattr(obj,'_callf'):
            obj._callf = obj.f 
        return obj
    
    
class FunctionModel(Model):
    __metaclass__ = _FuncMeta
    
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
    
    def __call__(self,x):
        """
        call the function on the input x with the current parameters and return
        the result
        """
        arr = np.array(x,copy=False,dtype=float)
        res = self._callf(np.atleast_1d(arr),*self.parvals)
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
        
        returns a pymc.MCMC object for Markov Chain Monte Carlo sampling
        """
        import pymc
        from operator import isSequenceType
        from inspect import getargspec
        from types import MethodType
        
        if set(priors.keys()) != set(self.params):
            raise ValueError("input priors don't match function params")
        
        d={}
        for k,v in priors.iteritems():
            if isinstance(v,pymc.StochasticBase):
                d[k]=v
            elif isSequenceType(v):
                if len(v) == 2:
                    d[k]=pymc.distributions.Uniform(k,v[0],v[1])
                else:
                    raise ValueError("couldn't understand sequence "+str(v) )
            elif v > 0:
                d[k] = pymc.distributions.Normal(k,self.pardict[k],1.0/v/v)
            elif v == 0:
                d[k] = pymc.distributions.Poisson(k,self.pardict[k])
            else:
                raise ValueError("couldn't interpret prior "+str(v))
        
        funcdict=dict(d)
        if type(self.f) is MethodType:
            xstr = getargspec(self.f)[0][1]
        else:
            assert callable(self.f),'object function is not a callable'
            xstr = getargspec(self.f)[0][0] 
        funcdict[xstr]=x
        funcdet=pymc.Deterministic(name='f',eval=self.f,parents=funcdict,doc="FunctionModel1D function")
        d['f'] = funcdet
        
        if type(datamodel) is tuple:
            distobj,dataarg,kwargs=datamodel
            if 'name' not in kwargs:
                kwargs['name'] = 'data'
                
            kwargs[dataarg] = funcdet
            kwargs['observed'] = True
            kwargs['value'] = y
            
            datamodel = distobj(**kwargs)
        else:
            if datamodel is None:
                sig = np.std(y)
            else:
                sig = datamodel
            datamodel = pymc.distributions.Normal('data',mu=funcdet,tau=1/sig/sig,observed=True,value=y)
        d[datamodel.__name__ ]=datamodel
        
        return pymc.MCMC(d)

class FunctionModel1D(FunctionModel):
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
    
    if f has a variable number of arguments, the first argument of the 
    constructor is taken to be the number of arguments that particular 
    instance should have, and the 'parsname' class attribute can be used
    to specify the prefix for the name of the parameters
    """
    ndims = (1,1)
    
    defaultIntMethod = 'quad'
    defaultInvMethod = 'brentq'
    xAxisName = None
    yAxisName = None
    fitteddata = None
    
    
    
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
            res = fmin(g,x0,tuple(self.parvals),**kwargs)
        elif method == 'fmin_powell':
            res = fmin_powell(g,x0,tuple(self.parvals),**kwargs)
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
        
        the full output is available is self.lastfit
        
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
                
            if fixedpars:
                for p in fixedpars:
                    i=ps.index(p)
                    del ps[i]
                    del v[i]
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
            
        self.lastfit = res
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
                plt.subplot(nr,nc,i+1)
                plt.hist(d[p],bins=n//20)
                plt.xlabel(p)
                plt.ylabel('N')
                plt.title('Median$=%3.3f$ $\\sigma=%3.3f$ Fit$=%3.3f$'%(np.median(d[p]),np.std(d[p]),self.pardict[p]))
        
        return d
        
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
            if self.fitteddata:
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
        
        #TODO:vectorize when needed
        
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
                from warnings import warn
                warn('Integration message: %s'%m)
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
    def integrateCircular(self,lower,upper,*args,**kwargs):
        """
        This is an alias for self.integrate with jacobian set for a azimuthally
        symmetric 2D radial profile
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*2.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*2.0*pi
        return self.integrate(lower,upper,*args,**kwargs)
    
    def integrateSpherical(self,lower,upper,*args,**kwargs):
        """
        This is an alias for self.integrate with jacobian set for a spherically
        symmetric 3D radial profile
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*x*4.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*x*4.0*pi
        return self.integrate(lower,upper,*args,**kwargs)
        
    def derivative(self,x,dx=1):
        """
        the derivative at x
        
        if overridden in a subclass, the signature should be
        derivative(self,x,dx=1) , but dx may be ignored
        """
        return (self(x+dx)-self(x))/dx
    
    def setCall(self,type=None,xtrans=None,ytrans=None,**kwargs):
        """
        sets the type of function evaluation to occur when the model is called
        
        type can be:
        *None: basic function evaluation
        *derivative: derivative at the location
        *integrate: integral - specify 'upper' or 'lower' kwarg and
                    the evaluation location will be treated as the
                    other bound.  If neither is given, lower=0 is assumed
        *integrateCircular: same as integrate, but using polar jacobian
        *integrateSpherical: same as integrate, but using spherical jacobian
        
        xtrans and ytrans are functions applied to either the input call (for x)
        or the output (for y) - they can also be strings:
        *'log':base-10 logarithm
        *'ln':base-e logarithm
        *'pow':10**
        *'exp':e**
        
        kwargs are passed into the type requested, and the call will occur on 
        
        note that if the model object is called directly in an overridden 
        method that is used in the new call, it probably won't work
        """
        from types import MethodType
        import inspect
        from functools import partial
        
        transmap={'log':np.log10,'ln':np.log,'pow':partial(np.power,10),'exp':np.exp}
        xts,yts = xtrans,ytrans
        if isinstance(xtrans,basestring):
            xtrans = transmap[xtrans]
        if isinstance(ytrans,basestring):
            ytrans = transmap[ytrans]
            
        if type is None:
            if xtrans and ytrans: 
                self._callf = MethodType(lambda self,x,*pars:ytrans(self.f(xtrans(x),*pars)),self,self.__class__)
            elif xtrans:
                self._callf = MethodType(lambda self,x,*pars:self.f(xtrans(x),*pars),self,self.__class__)
            elif ytrans:
                self._callf = MethodType(lambda self,x,*pars:ytrans(self.f(x,*pars)),self,self.__class__)
            else:
                self._callf = self.f
        else:
            try:
                newf = getattr(self,type)
                
                if not callable(newf):
                    raise AttributeError
            except AttributeError:
                raise AttributeError('function %s not found in %s'%(type,self))
                
            if 'integrate' in type:
                if 'upper' in kwargs and 'lower' in kwargs:
                    raise ValueError("can't do integral with lower and upper both specified")
                elif 'upper' in kwargs:
                    xkw = 'lower'
                elif  'lower' in kwargs:
                    xkw = 'upper'
                else: 
                    kwargs['lower']=0
                    xkw = 'upper'
                    #raise ValueError("can't do integral without lower or upper specified")
                
                #TODO:undo this once integrateSpherical is vectorized
                def vecf(x):
                    kwargs[xkw] = x
                    return newf(**kwargs)
                vecf = np.vectorize(vecf)
                
                
                if xtrans and ytrans: 
                    def callfunc(self,x,*pars):
                        self.parvals = pars
                        return ytrans(vecf(xtrans(x)))
                elif xtrans:
                    def callfunc(self,x,*pars):
                        self.parvals = pars
                        return vecf(xtrans(x))
                elif ytrans:
                    def callfunc(self,x,*pars):
                        self.parvals = pars
                        return ytrans(vecf(x))
                else:
                    def callfunc(self,x,*pars):
                        #TODO:test cost of par-setting
                        self.parvals = pars
                        return vecf(x)
                
            else:    
                
                fargs, fvarargs, fvarkws, fdefaults = inspect.getargspec(newf)
                    
                for arg in fargs:
                    if arg not in kwargs and arg !='self':
                        xkw = arg
                        break
                else:
                    raise ValueError('function with no arguments attempted in setCall')
                if xtrans and ytrans: 
                    def callfunc(self,x,*pars):
                        #TODO:test cost of par-setting
                        self.parvals = pars
                        kwargs[xkw] = xtrans(x)
                        return ytrans(newf(**kwargs))
                elif xtrans:
                    def callfunc(self,x,*pars):
                        #TODO:test cost of par-setting
                        self.parvals = pars
                        kwargs[xkw] = xtrans(x)
                        return newf(**kwargs)
                elif ytrans:
                    def callfunc(self,x,*pars):
                        #TODO:test cost of par-setting
                        self.parvals = pars
                        kwargs[xkw] = x
                        return ytrans(newf(**kwargs))
                else:
                    def callfunc(self,x,*pars):
                        #TODO:test cost of par-setting
                        self.parvals = pars
                        kwargs[xkw] = x
                        return newf(**kwargs)
                
            self._callf = MethodType(callfunc,self,self.__class__)
            
        self._callftype = (type,xts,yts)
    def getCall(self):
        """
        returns the type of evaluation to perform when this model is called - 
        a string like that of the type passed into `setCall`, or None if
        the model function itself is to be called.
        """
        if hasattr(self,'_callftype'):
            return self._callftype
        else:
            return None
  
class _CompMeta1D(_FuncMeta):
#    def __init__(cls,name,bases,dct):
#        super(_CompMeta1D,cls).__init__(name,bases,dct)
    def __call__(cls,*args,**kwargs):
        obj = super(_FuncMeta,_CompMeta1D).__call__(cls,*args,**kwargs) #object __init__ is called here
        
        #bind _callf to the standard call if it isn't messed with by the function constructor
        if not hasattr(obj,'_callf'):
            obj._callf = obj.f 
        return obj
        
class CompositeModel1D(FunctionModel1D):
    """
    This model contains a group of FunctionModel1D objects and evaluates them
    as a single model.
    
    The models can either be FunctionModel1D objects, FunctionModel1D classes,
    or a string (in the later two cases, new instances will be generated)
    
    parameter names are of the form 'A0' and 'A1' where A is the parameter
    name and the number is the sequential number of the model with that
    parameter.  If autoshorten is True, the suffix will be removed if there
    is only one of that parameter
    
    parnames is a dictionary that maps from names of the form 'A0' to a 
    different name for the parameter.  If a parameter with the name already
    exists, a ValueError will be raised.
    
    the operations can be any valid python operator, or a sequence of operators
    to apply (e.g. ['+','*','+'] will do mod1+mod2*mod3+mod4), or a string
    giving the expression to be evaluated, with  'm' in places where
    the evaluated models (in order) should go (e.g. 'm + m * m + m' will do
    mod1+mod2*mod3+mod4)
    
    
    
    any extra kwargs will be used to specify default values of the parameters
    """
    __metaclass__ = _CompMeta1D
    
    #TODO:initial vals
    def __init__(self,models=[],operation='+',parnames={},autoshorten=True,
                  **parvals):
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
            if len(operation.strip()) == 1:
                operation = [operation for i in range(len(mods)-1)]
            else:
                operation = operation.split('m')[1:-1]
        elif len(operation) != len(mods)-1:
            raise ValueError('impossible number of operations')
        
        
        self._ops = tuple(operation)
            
        oplst = ['mval[%i]%s'%t for t in enumerate(self._ops)]
        oplst += 'mval[%i]'%len(oplst)
        self._opstr = ''.join(oplst)
        
        args = []
        argmap = {}
        #add all parameters with suffixes for their model number
        for i,m in enumerate(self._models):
            for p in m.params:
                argn = p+str(i)
                args.append(argn)
                argmap[argn] = (i,p) 
                
        if parnames:
            for k,v in parnames.iteritems():
                try:
                    i = args.index(k)
                except ValueError:
                    raise KeyError('parameter %s not present'%k)
                if v in args:
                    raise ValueError('attempted to specify a replacement parameter name that already exists')
                args[i] = v
                argmap[v] = argmap[k]
                del argmap[k]
        else:
            parnames = {}
                        
        #remove suffixes if they are unique
        if autoshorten:
            cargs = list(args)
            for i,a in enumerate(cargs):
                if a not in parnames.values():
                    argname = argmap[a][1]
                    for j,a2 in enumerate(cargs):
                        if a2.startswith(argname) and i != j:
                            break
                    else:
                        argmap[argname] = argmap[a]
                        del argmap[a]
                        args[i] = argname
            
        self._args = tuple(args)
        self._argmap = argmap
        
        self._filters = None
        
        for k,v in parvals.iteritems():
            setattr(self,k,v)
    
    def __getattr__(self,name):
        if name in self._args:
            i,n = self._argmap[name]
            return getattr(self._models[i],n)
        raise AttributeError("'%s' has no attribute '%s'"%(self,name))
    
    def __setattr__(self,name,val):
        if name in self._args:
            i,n = self._argmap[name]
            setattr(self._models[i],n,val)
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
        return self.fitData(*args,**kwargs)
            
    @property
    def models(self):
        return self._models
    
    @property
    def ops(self):
        return self._ops  
    
class Grid1DModels(object):
    """
    A set of models that are matching types and represnt a grid of parameter
    values.  The main purpose is to allow the grid to be inverted to extract
    parameter values from the grid
    """
    
    def __init__(self,models,nmodels=None,extraparams={},**params):
        """
        models can be a sequence of models (must be the same type) or a model
        type.  If a type, the number of models to generate will be inferred
        from the parameters (or nmodels if the parameters are all scalar)
        """
        from operator import isSequenceType
        
        n = nmodels
        if isSequenceType(models) and not isinstance(models,basestring):
            for i,m in enumerate(models):
                if not isinstance(m,Model):
                    raise TypeError('Invalid model given as sequence index %i'%i)
            if n is not None and n != len(models):
                raise ValueError('n does not match the number of models')
            n = len(models)
            modcls= None
        else:
            modcls = get_model(models)
            d = params.copy()
            d.update(extraparams)
            for p,pv in d.iteritems():
                if isSequenceType(pv):
                    nn = len(pv)
                    if n is not None:
                        if n != nn:
                            raise ValueError('param %s does not match previous value'%p)
                    else:
                        n = len(pv)
        if n is None:
            raise ValueError('must specify n if nothing else is a sequence')
                
        if modcls is not None:
            models = [modcls() for i in range(n)]
            
        extraparams = dict([(p,pv) if isSequenceType(pv) else (p,pv*np.ones(n))for p,pv in extraparams.iteritems()])    
        params = dict([(p,pv) if isSequenceType(pv) else (p,pv*np.ones(n))for p,pv in params.iteritems()]) 
        
        for i,m in enumerate(models):
            for p,pv in params.iteritems():
                setattr(m,p,pv[i])
            
        self.models = models
        self.extraparams = extraparams
        
        self.interpmethod = 'lineary'
        
    def getParam(self,x,y,parname):
        """
        computes the value of the requested parameter that interpolates onto
        the provided point
        """        
        isscalar = np.isscalar(x)
        
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        if x.shape != y.shape:
            raise ValueError("x and y do not match")
        
        if self.interpmethod == 'lineary':
            modys = np.array([m(x) for m in self.models])
            if parname in self.extraparams:
                ps = self.extraparams[parname]
            else:
                ps = [getattr(m,parname) for m in self.models]
            
            res = np.array([np.interp(y[i],modys[:,i],ps) for i,xi in enumerate(x)])
        elif self.interpmethod == 'linearx':
            raise NotImplementedError('linearx not ready')
        else:
            raise NotImplementedError('invalid interpmethod')
        
        if isscalar:
            return res[0]
        else:
            return res
        
    def plot(self,*args,**kwargs):
        from matplotlib import pyplot as plt
        
        if kwargs.pop('clf',True):
            plt.clf()
        kwargs['clf'] = False
        for m in self.models:
            m.plot(*args,**kwargs)
            
        
#<-------------------------------Module functions ----------------------------->  
__model_registry={}
def register_model(model,name=None,overwrite=False,stripmodel=True):
    """
    register a model at the module package level for get_model and list_model
    
    model is the class object
    
    name is the name to assign (class name will be used if this is None)
    
    if overwrite is True, if a model already exists with the provided name, it
    will be silently overwritten, otherwise a KeyError will be raised
    
    if stripmodel is true, the characters 'model' will be stripped from
    the name before the model is assigned
    """
    if not issubclass(model,Model):
        raise TypeError('Supplied model is not a Model')
    if name is None:
        name = model.__name__.lower()
    else:
        name = name.lower()
    if stripmodel:
        name = name.replace('model','')
    if overwrite and name in __model_registry:
        raise KeyError('Model %s already exists'%name)
    
    __model_registry[name]=model
    
def get_model(model):
    """
    returns the class object for the requested model in the model registry
    
    the model can be a name, instance, or class object
    
    KeyError is raise if a name is provided and it is incorrect, all other 
    invalid inputs raise a TypeError
    """
    from inspect import isclass
    
    if isinstance(model,basestring):    
        return __model_registry[model]
    elif isclass(model):
        if issubclass(model,Model):
            return model
        else:
            raise TypeError('class object is not a Model')
    elif issubclass(model.__class__,Model):
        return model.__class__
    else:
        raise TypeError('attempted to get invalid model')

def list_models(ndims=None,include=None,exclude=None):
    """
    lists the registered Model objects in the package
    
    to get only models that accept a certain dimensionality, set ndims 
    to the model dimension wanted
    
    include is a sequence of model names to include in the list (e.g. the 
    function will just validate that they are valid models and return the
    strings) - if not in the registry, a ValueError will be raised
    
    exclude is a sequence of model names that should be excluded (if any 
    are not in the registry, a ValueError will be raised)
    """
    from operator import isSequenceType
    
    res = __model_registry.keys()
    
    if include is not None:
        if exclude is not None:
            raise TypeError("can't specify both included models and excluded models")
        if isinstance(include,basestring) or not isSequenceType(include):
            include = [include]
        for m in include:
            if m not in res:
                raise ValueError('modelname to include %s not a valid model name'%m)
            res = include
    elif exclude is not None:
        if isinstance(exclude,basestring) or not isSequenceType(exclude):
            exclude = [exclude]
        for m in exclude:
            if m not in res:
                raise ValueError('modelname to exclude %s not a valid model name'%m)
            res.remove(m)
    
    if ndims is not None:
        if np.isscalar(ndims):
            ndims = (ndims,ndims)
        res = [m for m in res if __model_registry[m].ndims == ndims]
    
    return res


def intersect_models(m1,m2,bounds=None,nsample=1024,full_output=False,**kwargs):
    """
    determine the points where two models intersect
    
    if bounds is None, the bounds will be determined from the model data if
    any is saved.  Otherwise, it should be (min,max) for the region to look
    for points
    
    returns a sorted array of points where the two models intersect on the 
    interval (up to a resolution in nsample), or if full_output is True,
    returns array,scipy.optimize.zeros.RootResults
    
    kwargs are passed into scipy.optimize.brentq
    """
    from scipy.optimize import brentq
    
    if bounds is None:
        data1 = m1.fitteddata if hasattr(m1,'fitteddata') else None
        data2 = m2.fitteddata if hasattr(m2,'fitteddata') else None
        
        if data1 is None and data2 is None:
            raise ValueError('must supply bounds if neither model has data')
        elif data1 is None:
            data2 = data2[0]
            bounds = (np.min(data2),np.max(data2))
        elif data2 is None:
            data1 = data1[0]
            bounds = (np.min(data1),np.max(data1))
        else: #both are valid
            data1 = data1[0]
            data2 = data2[0]
            bounds = (min(np.min(data1),np.min(data2)),max(np.max(data1),np.max(data2)))
    
    xsample = np.linspace(bounds[0],bounds[1],nsample)
    
    diff = m1(xsample) - m2(xsample)
    transis = np.convolve(diff>0,[1,-1]) #1-> crossing up between index and the one before
    transis[0] = 0
    transis[-1] = 0
    
    inters = []
    reses = []
    kwargs['full_output'] = 1
    for i in np.where(transis)[0]:
        t = brentq(lambda x:m1(x)-m2(x),xsample[i-1],xsample[i],**kwargs)
        inters.append(t[0])
        reses.append(t[1])
        
    arr = np.array(inters)
    arr.sort()
    
    if full_output:
        return arr,reses
    else:
        return arr
    
def offset_model(model,pname='C',**kwargs):
    """
    generator function that takes a Model type (class or name) and uses it
    to construct a CompositeModel with an extra parameter that is a constant
    offset (with a name given by the pname argument)
    
    kwargs are used to set the parameter values
    """
    return CompositeModel1D((model,ConstantModel),('+'),{'C1':pname},**kwargs)
    
def scale_model(model,pname='A',**kwargs):
    """
    generator function that takes a Model type (class or name) and uses it
    to construct a CompositeModel with an extra parameter that is a constant
    offset (with a name given by the pname argument)
    
    kwargs are used to set the parameter values
    """
    mod = CompositeModel1D((model,ConstantModel),('*'),{'C1':pname},**kwargs)

    if pname not in kwargs:
        setattr(mod,pname,1.0)
    return mod

def scale_and_offset_model(model,scalepname='A',offsetpname='C',**kwargs):
    """
    generator function that takes a Model type (class or name) and uses it
    to construct a CompositeModel that adds a 'C1' parameter that scales the
    model, and a 'C2' parameter than is the offset
    
    kwargs are used to set the parameters
    """
    mod = CompositeModel1D((model,ConstantModel,ConstantModel),('*','+'),{'C1':scalepname,'C2':offsetpname},**kwargs)
    if scalepname not in kwargs:
        setattr(mod,scalepname,1.0)
    return mod
    
    
def binned_weights(values,n,log=False):
    """
    Produces an array of values of the same size as values that is generated
    by subdividing the values into n bins such that each bin has an equal share
    of the total number of values.
    
    Returns an array of values on [0,1]
    
    if log is True, the values are evenly-spaced on logarithmic intervals
    """
    
    if log:
        values = np.log(values).ravel()
    else:
        values = np.array(values,copy=False).ravel()
        
    mx,mi = np.max(values),np.min(values)
    
    n,edges = np.histogram(values)
    ws = np.zeros_like(values)
    
    wsr = ws.ravel()
    for i,w in enumerate(1.0/n):
        m = (edges[i]<=values) & (values<edges[i+1])
        wsr[m] = w 
    wsr[edges[-1]==values] = w
    
    return ws

#<---------------------------------Builtin models------------------------------>
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
        
        lastfit stores ((m,b),dm,db,dy)
        """  
        if weights is not None:
            if fixedpars or len(kwargs)>0:
                from warnings import warn
                warn("customized exact linear fit not yet available for fixed pars")
                kwargs=kwargs.copy()
                kwargs['x']=x
                kwargs['y']=y
                kwargs['method']='leastsq'
                kwargs['fixedpars']=fixedpars
                kwargs['weights']=weights
                return FunctionModel1D.fitData(self,**kwargs)
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
        m,b = self.m,self.b
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
    
    parsname = 'c'
    
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
        
class MoffatModel(FunctionModel1D):
    """
    Moffat function given by:
    (beta-1)/(pi alpha^2) [1+(r/alpha)^2]^-beta
    """
    def f(self,r,alpha=1,beta=4.765):
        roa=r/alpha
        return (beta-1)/(pi*alpha**2)*(1+roa*roa)**-beta
    
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
    
class SinModel(FunctionModel1D):
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
    def f(self,x,a=1,b=2,C=0,xs=0):
        z = x-xs
        return a*z+(b-a)*np.log(1+np.exp(z))+C
    
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
        self.fitData(*self.fitteddata)
            
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
    
    
class InterpolatedSplineModel(FunctionModel1D):
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
        self.fitteddata=(np.arange(self.degree+1),np.arange(self.degree+1))
        self.fitData(*self.fitteddata)
            
    def _customFit(self,x,y,fixedpars=(),**kwargs):
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
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)
    
class _KnotSplineModel(FunctionModel1D):
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
        
        self.fitteddata=(np.arange(self.degree+self.nknots+1),np.arange(self.degree+self.nknots+1))
    
    @abstractmethod        
    def f(self,x):
        raise NotImplemetedError
    
    @abstractmethod    
    def _customFit(self,x,y,fixedpars=(),**kwargs):        
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

class UniformKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self._oldk = self._oldd = None
        super(UniformKnotSplineModel,self).__init__()
        self.fitData(*self.fitteddata)
    
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        self.iknots = np.linspace(x[0],x[-1],self.nknots+2)[1:-1]
        
        super(UniformKnotSplineModel,self)._customFit(x,y,fixedpars,**kwargs)
        
        self._oldk = self.nknots
        self._oldd = self.degree
        
        return np.array([self.nknots,self.degree])
    
    def f(self,x,nknots=3,degree=3):
        if self._oldk != nknots or self._oldd != degree:
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)
    
    

class UniformCDFKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self._oldk = self._oldd = None
        super(UniformCDFKnotSplineModel,self).__init__()
        self.fitData(*self.fitteddata)
    
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        cdf,xcdf = np.histogram(x,bins=max(10,max(2*self.nknots,int(len(x)/10))))
        mask = cdf!=0
        cdf,xcdf = cdf[mask],xcdf[np.hstack((True,mask))]
        cdf = np.hstack((0,np.cumsum(cdf)/np.sum(cdf)))
        self.iknots = np.interp(np.linspace(0,1,self.nknots+2)[1:-1],cdf,xcdf)
        
        super(UniformCDFKnotSplineModel,self)._customFit(x,y,fixedpars,**kwargs)
        
        self._oldk = self.nknots
        self._oldd = self.degree
        
        return np.array([self.nknots,self.degree])
    
    def f(self,x,nknots=3,degree=3):
        if self._oldk != nknots or self._oldd != degree:
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)

class SpecifiedKnotSplineModel(_KnotSplineModel):
    def __init__(self):
        self.nknots = self.__class__._currnparams
        self._oldd = None #this will force a fit at first call
        super(SpecifiedKnotSplineModel,self).__init__()
        
        self.setKnots(np.linspace(-1,1,self.nknots))
    
    def _customFit(self,x,y,fixedpars=(),**kwargs):
        """
        just fits the spline with the current s-value - if s is not changed,
        it will execute very quickly after
        """
        self.iknots = np.array([v for k,v in self.pardict.iteritems() if k.startswith('k')])
        self.iknots.sort()
        
        super(SpecifiedKnotSplineModel,self)._customFit(x,y,fixedpars,**kwargs)
        
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
    
    parsname = 'k'
    
    degree=3 #default cubic
    def f(self,x,degree,*args):
        #TODO:faster way to do the arg check?
        if self._oldd != degree or np.any(self.iknots != np.array(args)):
            xd,yd = self.fitteddata
            self._customFit(xd,yd,weights=self._ws)
        
        return self.spline(x)
    
    
    
class NFWModel(FunctionModel1D):
    """
    A Navarro, Frenk, and White 1996 profile
    """
        
    def f(self,x,rho0=1,rc=1):
        #return TwoPowerModel.f(self,x,rho0*rc*rc*rc,rc,-1,-3)
        return rho0*rc*rc*rc*((x+rc)**(-2))*(x**-1)
    
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
        sets the model parameters to match a given concentration 
        
        if Rvir or Mvir are None, the Rvir/Mvir relation in this model 
        (Maller&Bullock 2004) will be used to infer the relation
        """
        if Rvir is None and Mvir is None:
            raise ValueError('Must specify Rvir,Mvir, or both')
        elif Rvir is None:
            Rvir = self.Mvir_to_Rvir(Mvir)
        elif Mvir is None:
            Mvir = self.Rvir_to_Mvir(Rvir)
        
        self.rc = Rvir/c
        
        self.rho0 = 1
        a0 = self.integrateSpherical(0,Rvir)
        self.rho0 = Mvir/a0
            
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
        return 1e12/h*(Omega0*NFWModel.Delta(z)/97.2)*(Rvir*(1+z)/203.4/h)**3
    
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

class King2DrModel(FunctionModel1D):    
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        return A*rcsq*((r*r+rcsq)**-0.5 - (rt*rt+rcsq)**-0.5)**2
    
class King3DrModel(FunctionModel1D):
    def f(self,r,rc=1,rt=2,A=1):
        rcsq=rc*rc
        z=((r*r+rcsq)**0.5) * ((rt*rt+rcsq)**-0.5)
        return (A/z/z/pi/rc)*((1+rt*rt/rcsq)**-1.5)*(np.arccos(z)/z-(1-z*z)**0.5)

class SchecterMagModel(FunctionModel1D):
    def f(self,M,Mstar=-20.2,alpha=-1,phistar=1.0857362047581294):
        from numpy import log,exp
        x=10**(0.4*(Mstar-M))
        return 0.4*log(10)*phistar*(x**(1+alpha))*exp(-x)
    
class SchecterLumModel(FunctionModel1D):
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
        
        print '?',lower,upper
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
    
    
class Plane(FunctionModel):
    """
    Generates a plane that follows the form
    d = a*x+b*y+c*z (e.g. (a,b,c) is the normal vector and 
    d/a, b ,or c are the intercepts
    """
    ndims = (2,1)
    
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
    
#register all Models in this module
for o in locals().values():
    if type(o) == type(FunctionModel1D) and issubclass(o,FunctionModel1D) and o != FunctionModel1D and o!= CompositeModel1D and not o.__name__.startswith('_'):
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