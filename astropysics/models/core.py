#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module holds the (mostly abstract) classes for the data modeling framework
used in astropysics, seperated from the implementations of specific models found
in the :mod:`builtins` module.


Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.models.core
   :parts: 1
   
Module API
^^^^^^^^^^

"""

#TODO: refactor so that FunctionModels that interpolate over data have a common superclass, then add linear interpolator, polynomial interpolator
#TODO: refine ModelGrid, tying into util grid inverter

from __future__ import division,with_statement
from ..constants import pi
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

class ModelTypeError(Exception):
    """
    This exception indicates a problem with the value of the input 
    or output of a model. 
    """
    def __init__(self,message):
        super(ModelTypeError,self).__init__(message)
    
class ParametricModel(object):
    """
    The superclass of all models with parameters. Subclasses should implement
    abstract properties and methods:
    
    * :meth:`__call__`
        Takes a single argument as the model input, returns the model output
    * :meth:`params`
        a sequence of names for the parameters of the model
    * :meth:`parvals`
        a sequence of values for the parameters of the model
    
    Optional overrides:
    
    * :meth:`pardict`
        a dictionary with keys as parameter names and values as the value for
        that parameter
    * :meth:`inv`
        compute the inverse of the model
        
    """
    __metaclass__ = ABCMeta
    
    _data = None
    
    @abstractmethod
    def __call__(self,x):
        raise NotImplementedError
       
    params = abstractproperty(doc='a sequence of the parameter names')     
    parvals = abstractproperty(doc='a sequence of the values in the parameters')       
    
    def _getPardict(self):
        return dict([t for t in zip(self.params,self.parvals)]) 
    def _setPardict(self,val):
        self.parvals = [val[p] for p in self.params]
    pardict = property(_getPardict,_setPardict,doc="""
        A dictionary of the parameter names and values.
        """)
        
    def _getData(self):
        return self._data
    def _setData(self,val):
        if val is None:
            self._data = None
        else:
            try:
                if not (2 <= len(val) <= 3):
                    raise ValueError('input model data in invalid form - must be 2 or 3-sequence')
                else:
                    ind = np.array(val[0],copy=False)
                    outd = np.array(val[1],copy=False)
                    if len(val) < 3 or val[2] is None:
                        ws = None
                    else:
                        ws = np.array(val[2],copy=False)
                    self._data = (ind,outd,ws)
            except TypeError,e:
                e.args = ('invalid type for model data',)
                raise
    data = property(_getData,_setData,doc='data should be either None, or (datain,dataout,weights)')
    
    def inv(self,output,*args,**kwargs):
        """
        Compute the inverse of this model for the requested output
        
        It should raise a ModelTypeError if the model is not invertable for
        the provided data set
        """
        raise ModelTypeError('Model is not invertible')
    
    
    

class _AutoParameter(object):
    """
    Parameter representation used for objects instantiated with the
    AutoParamsMeta magic.
    """
    def __init__(self,name,defaultvalue):
        self.name = name
        self.defaultvalue = defaultvalue
    
    def __get__(self,obj,owner):
        try:
            return getattr(obj,'_AutoParameter__'+self.name)
        except AttributeError: #if data doesn't exist, initialize to default and return it
             setattr(obj,'_AutoParameter__'+self.name,self.defaultvalue)
             return self.defaultvalue
    
    def __set__(self,obj,value):
        if isinstance(value,np.ndarray) and value.shape == ():
            pass
        elif not np.isscalar(value):
            raise ValueError('Parameters must be scalar values - type '+str(type(value))+' provided')
        setattr(obj,'_AutoParameter__'+self.name,value)
    
    def __delete__(self,obj):
        raise AttributeError("can't delete a parameter")

class AutoParamsMeta(ABCMeta):
    """
    Metaclass used to auto-generate parameters from the :meth:`FunctionModel.f`
    method of :class:`FunctionModel` subclasses.
    
    This generates attributes for each of the inputs of the method
    :meth:`FunctionModel.f` (except 'self' and 'x')
    
    When a model is instantiated, the arguments specify initial values for the
    parameters, and any non-parameter kwargs will be passed into the
    :meth:`__init__` method of the model class.
    
    If f is of the form f(self,x,*args), the first argument of the constructor
    is taken to be the number of arguments that particular instance should have,
    and the :attr:`paramnames` class attribute can be used to specify the prefix
    for the name of the parameters (and can be an iterator, in which case one of
    each name will be made for each var arg, in sequence order). default values
    can be given by adding class variables of the form :attr:`_param0_default`
    """
    def __init__(cls,name,bases,dct):
        #called on import of astro.models
        from inspect import getargspec
        
        super(AutoParamsMeta,cls).__init__(name,bases,dct)
        for c in cls.mro():
            if 'f' in c.__dict__:
                ffunc = c.__dict__['f']
                if callable(ffunc):
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
            setattr(cls,arg,_AutoParameter(arg,d))
            
        if varargs is not None:
            cls._pars = None
            cls._statargs = list(args)
        else:
            cls._pars = tuple(args)
            
    def __call__(cls,*args,**kwargs):
        if cls._pars is None: #this is the case for *args in the function
            if len(args)>1:
                raise ValueError('must pass parameters with keywords at initializer if variable argument numbers are used')
            
            if hasattr(cls,'paramnames'):
                paramsnames = getattr(cls,'paramnames')
            else:
                paramsnames =  ('p',)
            if isinstance(paramsnames,basestring):
                paramsnames = (paramsnames,)
                
            if hasattr(cls,'paramoffsets'):
                paramsoffsets = getattr(cls,'paramoffsets')
            else:
                paramsoffsets = 0
            
            try:
                if 'nparams' in kwargs:
                    nparams = kwargs.pop('nparams')
                else:
                    args = list(args)
                    nparams = args.pop(0)
            except IndexError:
                raise IndexError('No # of parameters found for variable-size function')
            objkwargs = dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._statargs if not any([k.startswith(paramsname) for paramsname in paramsnames])])
            cls._currnparams = nparams #this is set in case the constructor needs to know for some reason
            obj = super(AutoParamsMeta,AutoParamsMeta).__call__(cls,**objkwargs) #object __init__ is called here
            del cls._currnparams
            pars = list(cls._statargs)
            
            for i in range(nparams):
                for paramsname in paramsnames:
                    p = paramsname+str(i+paramsoffsets)
                    pars.append(p)
                    if hasattr(cls,'_'+p+'_default'): #let class variables of the form "_param_default" specify defaults for varargs
                        setattr(obj,p,getattr(cls,'_'+p+'_default'))
                    else:
                        setattr(obj,p,1) #default varargs to 1
            obj._pars = tuple(pars)
        else: #this is the case for fixed functions
            objkwargs=dict([(k,kwargs.pop(k)) for k in kwargs.keys() if k not in cls._pars])
            try:
                obj = super(AutoParamsMeta,AutoParamsMeta).__call__(cls,**objkwargs) #object __init__ is called here
            except TypeError,e:
                if len(e.args)>0 and 'object.__new__() takes no parameters' == e.args[0]:
                    raise TypeError('invalid parameter in model constructor:'+str(objkwargs.keys()))
                else:
                    raise
            
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
    
    
class FunctionModel(ParametricModel):
    """
    A base class for models that use a single function to map a
    numpy array input to an output (the function must be specified as the 
    "f" method)
    
    The following attributes *must* be defined:
    
    * :meth:`f`
        A method that takes an array as the first argument and returns another
        array.
    * :attr:`_pars`
        A sequence of strings that define the names of attributes that are to be
        treated as the parameters for the model. If these are not defined, they
        will be initialized to the :attr:`defaultparval` value
    
    * :meth:`_filterfunc`
        A function that performs any processing of inputs and outputs. Should
        call :meth:`f` and have the same type of return value.
      
    If type-checking is to be performed on the input, it should be performed in
    :meth:`__call__` (see docstring for :meth:`__call__` for syntax), generally
    using :func:`astropysics.utils.check_type`, but any mandatory input
    *conversion* should be done in :meth:`_filterfunc`.
    """
    
    defaultparval = 1
    """
    The base default value if a parameter cannot get a default any other way.
    """
    
    fixedpars = tuple()
    """
    A sequence of the parameter names that by default should be kept fixed.
    """
    
    @abstractmethod
    def f(self,x,*params):
        """
        This abstract method *must* be overriden in subclasses. It is
        interpreted as the function that takes an array (or array-like) as the
        first argument and returns an array of output values appropriate for the
        model.
        
        model parameters are passed in as the rest of the arguments, in the 
        order they are given in the _pars sequence.
        """
        raise NotImplementedError
    
    
    def _filterfunc(self,*args,**kwargs):
        """
        This single-use function is present to avoid need for setting the filter
        function in an initializer.
        """
        self._filterfunc = self.f
        return self.f(*args,**kwargs)
    
    def __call__(self,x):
        """
        Call the model function on the input x with the current parameters and
        return the result.
    
        :except ModelTypeError:
            If the input or output are incorrect type or dimensionality for this
            model.
        
        ..note ::
            If a subclass overrides this method to do type-checking, it should 
            either call this method or call :meth:`_filterfunc` with an array
            input and raise a ModelTypeError if there is a type problem 
        """
        
        arr = np.array(x,copy=False,dtype=float)
        return self._filterfunc(arr,*self.parvals)
    
    _pars=tuple() #default empty parameter sequence
        
    @property
    def params(self):
        """
        A tuple of the parameter names. (read-only)
        """
        return self._pars
    
    def _getParvals(self):
        try:
            return tuple([getattr(self,a) for a in self._pars])
        except AttributeError:
            self._autoInitPars()
            return self._getParvals()
    def _setParvals(self,newpars):
        if len(newpars)>len(self._pars):
            raise ValueError('too many new parameters')
        for a,v in zip(self._pars,newpars):
            setattr(self,a,v)        
    parvals = property(_getParvals,_setParvals,doc="""
    The values of the parameters in the same order as :attr:`params`
    """)
    
    def _getPardict(self):
        try:
            return dict([(a,getattr(self,a)) for a in self._pars])
        except AttributeError:
            self._autoInitPars()
            return self._getPardict()
    def _setPardict(self,newpardict):
        for k in newpardict.keys():
            if k not in self._pars:
                raise KeyError('No parameter %s found'%k)
        for k,v in newpardict.iteritems():
            setattr(self,k,v)
    pardict = property(_getPardict,_setPardict,doc="""
    A dictionary mapping parameter names to the associated values.
    """)
    def _autoInitPars(self):
        """
        Called to generate parameters if they are missing.
        """
        for p in self._pars:
            setattr(self,p,self.defaultparval)
    
    fittype ='leastsq'
    """
    The currently selected fitting technique.
    """
    
    _optfittypes = ('leastsq','fmin','fmin_powell','fmin_cg','fmin_bfgs',
                 'fmin_ncg','anneal','global','brute')
    @property
    def fittypes(self):
        """
        A Sequence of the available valid values for the :attr:`fittype`
        attribute.  (Read-only)
        """
        ls = list()
        for cls in self.__class__.__mro__:
            if hasattr(cls,'_fittypes'):
                ls.extend(cls._fittypes)
        ls.extend(self._optfittypes)
        return tuple(ls)
    
    def fitData(self,x=None,y=None,fixedpars='auto',weights=None,savedata=True,
                 updatepars=True,fitf=False,contraction='sumsq',**kwargs):
        """
        Fit the provided data using algorithms from scipy.optimize, and adjust
        the model parameters to match.
        
        The fitting technique is sepcified by the :attr:`fittype` attribute of
        the object, which by default can be any of the optimization types in the
        :mod:`scipy.optimize` module (except for scalar minimizers)
        
        The full fitting output is available in :attr:`lastfit` attribute after
        this method completes.
        
        :param x: 
            The input values at which to evaluate the model. Valid shapes are
            those that this model will accept.
        :type x: array-like
        :param y: 
            The expected output values for the model at the given `x` values.
            Valid shapes are those that this model will output.
        :type y: array-like
        :param fixedpars: 
            Parameter names to leave fixed. If 'auto' the fixed parameters are
            inferred from self.fixedpars (or all will be free parameters if
            self.fixedpars is absent). If None, all parameters will be free.
        :type fixedpars: sequence of strings, 'auto' or None
        :param weights: 
            Weights to use for fitting, statistically interpreted as inverse
            errors (*not* inverse variance). May be one of the following forms:
            
            * None for equal weights
            * an array of points that must match the output
            * a 2-sequence of arrays (xierr,yierr) such that xierr matches the
              x-data and yierr matches the y-data
            * a function called as f(params) that returns an array of weights 
              that match one of the above two conditions

        :param savedata: 
            If True, `x`,`y`,and`weights` will be saved to self.data. Otherwise,
            data will be discarded after fitting.
        :type savedata: bool
        :param updatepars:  
            If True, sets the parameters on the object to the best-fit values.
        :type updatepars: bool
        :param fitf: 
            If True, the fit is performed directly against the :meth:`f` method
            instead of against the model as evaluated if called (as altered
            using :meth:`setCall`).
        :type fitf: bool
        :param contraction: 
            Only applies for optimize-based methods and is the technique used to
            convert vectors to figures of merit. this is composed of multiple
            string segments:
            
            #. these will be applied to each element of the vector first:
        
               - 'sq': square
               - 'abs':absolute value
               - '':raw value
              
            #. the contraction is performed after:
            
               - 'sum': sum of all the elements
               - 'median': median of all the elements
               - 'mean': mean of all the elements
               - 'prod': product of all the elements 
               
            * optionally,the string can include 'frac' to use the fractional 
              version of the differnce vector instead of raw values.  For the 
              leastsq method, this is the only applicable value
        
        `kwargs` are passed into the fitting function.
        
        :returns: array of the best fit parameters
        
        :except ModelTypeError: 
            If the output of the model does not match the shape of `y`.
        
        .. seealso:: :meth:`getMCMC`
        """
        from scipy import optimize as opt
        from operator import isMappingType
        from functools import partial
        
        self._fitchi2 = None #clear saved chi-squared if it exists
        
        if x is None:
            if hasattr(self,'data') and self.data is not None:
                x = self.data[0]
            else: 
                raise ValueError('No x data provided and no fitted data already present')
        else:
            x = np.array(x,copy=False)
        
        if y is None:
            if hasattr(self,'data') and self.data is not None:
                y = self.data[1]
            else: 
                raise ValueError('No y data provided and no fitted data already present')
        else:
            y = np.array(y,copy=False)
         
        if fitf:
            fitfunc = self.f
        else:
            fitfunc = self._filterfunc
        
        if fitfunc(x,*self.parvals).shape != y.shape:
            raise ModelTypeError('y array does not match output of model for input x')
        
        y = y.ravel()
        
        if self.fittype is None:
            method = self.fittypes[0]
        else:
            method = self.fittype
            
        if fixedpars is 'auto':
            fixedpars = self.fixedpars if hasattr(self,'fixedpars') else ()
        if fixedpars is None:
            fixedpars = tuple()
            
        ps=list(self.params)
        v=list(self.parvals) #initial guess
        
        if method not in self._optfittypes:
            for cls in self.__class__.__mro__:
                if hasattr(cls,'_fittypes') and isMappingType(cls._fittypes):
                    if method in cls._fittypes:
                        fitter = partial(cls._fittypes[method],self)
                        break
            else:
                fitter = 'fit'+method[0].upper()+method[1:]
                if hasattr(self,fitter):
                    fitter = getattr(self,fitter)
                else:
                    raise ValueError('could not locate fitting function for fitting method '+method)
            
            res = fitter(x,y,fixedpars=fixedpars,weights=weights,**kwargs) 
            
            #ensure that res is at least a tuple with parameters in elem 0
            from operator import isSequenceType
            if len(res)==0 or not isSequenceType(res[0]):
                res = (res,)
                
            if fixedpars:
                for p in fixedpars:
                    i=ps.index(p)
                    del ps[i]
                    del v[i]
        else:
            if weights is None:
                wf = lambda v:1
            elif callable(weights):
                wf = weights
            else:
                w = np.array(weights,copy=False)
                if w.shape == y.shape:
                    w = w.ravel()
                elif w.shape[1:] == y.shape and w.shape[0]==2:
                    w = (w[0]**2+w[1]**2)**0.5
                else:
                    raise ModelTypeError('weights shape does not match y')
                
                wf = lambda v:w
                
            kwargs['full_output'] = 1
            
            if fixedpars:
                for p in fixedpars:
                    i=ps.index(p)
                    del ps[i]
                    del v[i]
                
                #make a function of signature f(x,v) where v are the parameters to be fit
                pdict=dict([(p,getattr(self,p)) for p in fixedpars])
                def f(x,v):
                    pdict.update(dict(zip(ps,v)))
                    #return fitfunc(x,**pdict)
                    params = [pdict[a] for a in self._pars]
                    return fitfunc(x,*params).ravel()
            else:
                f=lambda x,v:fitfunc(x,*v).ravel()
                
            if method == 'leastsq':
                if 'frac' in contraction:
                    g=lambda v,x,y:wf(v)*(1-f(x,v)/y)
                else:
                    g=lambda v,x,y:wf(v)*(y-f(x,v))
                res=opt.leastsq(g,v,(x,y),**kwargs)
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
                            return np.diff
                if 'sum' in contraction:
                    g=lambda v,x,y:np.sum(wf(v)*g1(v,x,y),axis=None)
                elif 'mean' in contraction:
                    g=lambda v,x,y:np.mean(wf(v)*g1(v,x,y),axis=None)
                elif 'median' in contraction:
                    g=lambda v,x,y:np.median(wf(v)*g1(v,x,y),axis=None)
                elif 'prod' in contraction:
                    g=lambda v,x,y:np.prod(wf(v)*g1(v,x,y),axis=None)
                else:
                    raise ValueError('no valid contraction method provided')
                    
                if method == 'fmin':
                    res=opt.fmin(g,v,(x,y),**kwargs)
                elif method == 'fmin_powell':
                    res=opt.fmin_powell(g,v,(x,y),**kwargs)
                elif method == 'fmin_cg':
                    #TODO:smartly include derivative
                    res=opt.fmin_cg(g,v,args=(x,y),**kwargs)
                elif method == 'fmin_bfgs':
                    #TODO:smartly include derivative
                    res=opt.fmin_bfgs(g,v,args=(x,y),**kwargs)
                elif method == 'fmin_ncg':
                    raise NotImplementedError
                    #TODO:needs gradient and hessian
                    opt.fmin_ncg
                elif method == 'anneal' or method == 'global':
                    res=opt.anneal(g,v,args=(x,y),**kwargs)
                elif method == 'brute':
                    raise NotImplementedError
                    #TODO: set parrange smartly
                    res=opt.brute(g,parrange,(x,y),**kwargs)
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
            self.data = (x,y,weights)
        
        return v
    
    def stdData(self,x=None,y=None):
        """
        Determines the standard deviation of the model from the supplied data
        
        :param x: Input data value or None to use stored :attr:`data`
        :type x: array-like or None
        :param y: Output data value or None to use stored :attr:`data`
        :type y: array-like or None
        
        
        :returns: standard deviation of model from `y`
        
        """
        if x is None or y is None:
            if self.data is None:
                raise ValueError('must either specify data or save fitted data')
            x,y,weights = self.data
        
        if self(x).shape != y.shape:
            raise ModelTypeError('y array does not match output of model for input x')
        
        return np.std(y-self(x),ddof=len(self.params))
    
    def residuals(self,x=None,y=None,retdata=False):
        """
        Compute residuals of the provided data against the model, e.g.
        :math:`y-{\\rm model}(x)`.
        
        :param x: Input data value or None to use stored :attr:`data`
        :type x: array-like or None
        :param y: Output data value or None to use stored :attr:`data`
        :type y: array-like or None
        :param retdata: If True, returns the data along with the model.
        :type retdata: bool
        
        :returns: 
            Residuals of model from `y` or if `retdata` is True, a tuple
            (x,y,residuals).
        :rtype: array-like
        
        """
        if x is None or y is None:
            if self.data is None:
                raise ValueError('must either specify data or save fitted data')
            x,y,weights = self.data
        
        if self(x).shape != y.shape:
            raise ModelTypeError('y array does not match output of model for input x')
        if retdata:
            return x,y,y-self(x)
        else:
            return y-self(x)
    
    _fitchi2 = None #option to save so that fitData can store a chi2 if desired
    def chi2Data(self,x=None,y=None,weights=None):
        """
        Computes the chi-squared statistic for the data assuming this model.
        
        :param x: Input data value or None to use stored :attr:`data`
        :type x: array-like or None
        :param y: Output data value or None to use stored :attr:`data`
        :type y: array-like or None
        :param weights: 
            Weights to adjust chi-squared, typically for error bars.
            Statistically interpreted as the inverse error (*not* inverse
            variance). If None, any stored :attr:`data` will be used.
        :type weights: array-like or None
        
        :returns: tuple of floats (chi2,reducedchi2,p-value)
        
        If both are None, the internal data is used. In some cases the
        chi-squared statistic may be pre-computed in the fitting step rather
        than in this method.
        
        weights are taken here to be inverse-error
        
        
        returns chi2,reducedchi2,p-value
        """
        from scipy.stats import chisqprob
        
        chi2 = None
        if x is None or y is None:
            if self.data is None:
                raise ValueError('must either specify data or save fitted data')
            x,y,weights = self.data
            if self._fitchi2 is not None:
                chi2 = self._fitchi2
        
        if self(x).shape != y.shape:
            raise ModelTypeError('y array does not match output of model for input x')
        
        n=len(y)
        m=len(self.params)
        dof=n-m-1 #extra 1 because the distribution must sum to n
        
        if chi2 is None:
            if weights is None:
                mody=self(x)
                chi2 = np.sum((y-mody)**2/mody)
            else:
                while len(weights.shape)>1:
                    weights = np.sum(weights*weights,axis=0)**0.5
                mody=self(x)
                #assumes weights are inverse error, not ivar
                chi2 = np.sum((weights*(y-mody))**2)
        
        return chi2,chi2/dof,chisqprob(chi2,dof)
    
    def resampleFit(self,x=None,y=None,xerr=None,yerr=None,bootstrap=False,
                          modely=False,n=250,prefit=True, medianpars=False,
                          plothist=False,**kwargs):
        """
        Uses the fitData function to fit the function many times while 
        either using the  "bootstrap" technique (resampling w/replacement),
        monte carlo estimates for the error, or both to estimate the error
        in the fit.
        
        :param x: 
            The input data - if None, will be taken from the :attr:`data`
            attribute.
        :type x: array-like or None
        :param y: 
            The output data - if None, will be taken from the :attr:`data`
            attribute.
        :type y: array-like or None
        :param xerr: 
            Errors for the input data (assumed to be normally distributed), or
            if None, will be taken from the :attr:`data`. Alternatively, it can
            be a function that accepts the input data as the first argument and
            returns corresponding monte carlo sampled values.
        :type xerr: array-like, callable, or None
        :param yerr: 
            Errors for the output data (assumed to be normally distributed), or
            if None, will be taken from the :attr:`data`. Alternatively, it can
            be a function that accepts the output data as the first argument and
            returns corresponding monte carlo sampled values.
        :type yerr: array-like, callable, or None
        
        :param bootstrap:  
            If True, the data is also resampled (with replacement).
        :type bootstrap: bool
        
        :param modely:  
            If True, the fitting data will be generated by offseting from model
            values (evaluated at the `x`-values) instead of `y`-values.
        
        :param n: The number of times to draw samples.
        :type n: int
        
        :param prefit: 
            If True, the data will be fit without resampling once before the
            samples are recorded.
        :type prefit: bool
        
        :param medianpars:
            If True, the median from the histogram will be set as the value for
            the parameter.
        :type medianpars: bool
        
        :param plothist: 
            If True, histograms will be plotted using
            :func:`matplotlib.pyplot.hist` for each of the parameters, or if it
            is a string, only the histogram for the requested parameter will be
            shown.
        :type plothist: bool or string
        
        kwargs are passed into fitData
        
        :returns: 
            (histd,cov) where histd is a dictionary mapping parameters to their
            histograms and cov is the covariance matrix of the parameters in
            parameter order.
        
        .. note::
            If `x`, `y`, `xerr`, or `yerr` are provided, they do *not* overwrite
            the stored :attr:`data`, unlike most other methods for this class.
        """
        if x is None:
            if self.data is None:
                raise ValueError('must either specify data or use saved data')
            x = self.data[0]
        if y is None:
            if self.data is None:
                raise ValueError('must either specify data or use saved data')
            y = self.data[1]
        
        if xerr is None and self.data is not None and self.data[2] is not None:
            xerr = self.data[2][0]
        if yerr is None and self.data is not None and self.data[2] is not None:
            yerr = self.data[2][1]
        
        kwargs.setdefault('updatepars',False)
        
        if prefit:
            preupdate = kwargs['updatepars']
            kwargs['updatepars'] = True
            self.fitData(x,y,**kwargs)
            kwargs['updatepars'] = preupdate
        
        if xerr is None:
            xerr = lambda cen:cen
        elif callable(xerr):
            pass
        else:
            xsig = xerr
            xerr = lambda cen:np.random.normal(cen,xsig)
            
        if yerr is None:
            yerr = lambda cen:cen
        elif callable(yerr):
            pass
        else:
            ysig = yerr
            yerr = lambda cen:np.random.normal(cen,ysig)
        
        N=len(x)
        weights = kwargs.pop('weights',None)
        if weights is not None:
            weights = np.array(weights,copy=False)
        
        vs = []
        for i in range(n):
            xnew = xerr(x)
            ynew = yerr(self(x)) if modely else yerr(y)
            
            if bootstrap:
                sampi = (np.random.rand(N)*N).astype(int)
                xnew,ynew = xnew[sampi],ynew[sampi]
                if weights is not None:
                    #only flip last axis
                    kwargs['weights'] = weights.T[sampi].T
                
            vs.append(self.fitData(xnew,ynew,**kwargs))
            
        d=dict([(p,v) for p,v in zip(self.params,np.array(vs).T)])
        
        if plothist:
            import matplotlib.pyplot as plt
            params = (plothist,) if isinstance(plothist,basestring) else self.params
            
            #want to generate as square as possible
            nr=nc=int(len(params)**0.5)
            if nr*nc < len(params):
                nr+=1
            if nr*nc < len(params):
                nc+=1
            plt.clf()
                
            for i,p in enumerate(self.params):
                plt.subplot(nr,nc,i+1)
                plt.hist(d[p],bins=n//20,histtype='stepfilled')
                med = np.median(d[p])
                std = np.std(d[p])
                plt.axvline(self.pardict[p],c='r',ls='-')
                plt.axvline(med,c='k',ls='--')
                plt.axvline(med+std,c='k',ls=':')
                plt.axvline(med-std,c='k',ls=':')
                plt.xlabel(p)
                plt.ylabel('N')
                plt.title('Median$=%3.3f$ $\\sigma=%3.3f$ Fit$=%3.3f$'%(med,std,self.pardict[p]))
                
        #self.data = (x,y,(xerr,yerr))
        if medianpars:
            for p,v in d.iteritems():
                setattr(self,p,np.median(v))
        return d,np.cov([d[p] for p in self.params])
    
    def getMCMC(self,x,y,priors={},datamodel=None):
        """
        Generate an object to fit the data using Markov Chain Monte Carlo
        sampling. This function requires the `PyMC
        <http://code.google.com/p/pymc/>`_ package for the MCMC internals and
        sampling.
        
        :param x: Input data value
        :type x: array-like
        :param y: Output data value
        :type y: array-like
        :param priors: 
            Maps parameter names to the priors to assume for that parameter.
            There *must* be an entry for every parameter in the model. The prior
            specification can be in any of the following forms:
            
            * A :class:`pymc.Stochastric` object
            * A 2-tuple (lower,upper) for a uniform prior
            * A scalar > 0 to use a gaussian prior of the provided width 
              centered at the current value of the parameter
            * 0 for a Poisson prior with k set by the current value of the 
              parameter
            
        :type priors: dictionary
        :param datamodel:
            Specifies the model to assume for the fitting data points.  May
            be any of the following:
            
            * None
                A normal distribution with sigma given by the data's standard
                deviation.
            * A tuple (dist,dataname,kwargs)
                The first element is the pymc.distribution to be used as the
                distribution representing the data and the second is the name of
                the argument to be associated with the FunctionModel1D's output,
                and the third is kwargs for the distribution ("observed" and
                "data" will be ignored, as will the data argument)
            * A sequence
                A normal distribution is used with sigma for each data point
                specified by the sequence. The length must match the model.
            * A scalar
                A normal distribution with the given standard deviation.
        
        :except ValueError: If a prior is not provided for any parameter.
        
        :returns: A :class:`pymc.MCMC` object ready to sample for this model.
        
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
    
class CompositeModel(FunctionModel):
    """
    This model contains a group of :class:`FunctionModel` objects joined by
    standard arithmetic operations, and evaluates them as a single model. The
    models themselves are called, rather than the :meth:`f` function (and hence
    will be influenced by anything like :meth:`FunctionModel1D.setCall` calls).
    
    Generated parameter names are of the form 'A0' and 'A1' where A is the
    parameter name and the number is the sequential number of the model with
    that parameter. If `autoshorten` is True, the suffix will be removed if
    there is only one of that parameter.
    
    Note that no checking is performed here to ensure the model outputs are
    compatible - this can be done in subclasses for specific types of
    :class:`FunctionModels <FunctionModel>`.
    
    """
    
    #TODO:initial vals
    def __init__(self,models=[],operation='+',parnames={},autoshorten=True,
                  **parvals):
        """
        :param models: 
            The models objects to combine, or model types (in which case new
            models will be created)
        :type models: 
            sequence of :class:`FunctionModel` objects, :class:`FunctionModel`
            classes, or strings.
            
        :param operation: 
            The arithmetic operation(s) to join the models. (e.g. ['+','*','+']
            will do mod1+mod2*mod3+mod4). A single operator string will join all
            models with the same operator. Alternatively, a string of the form
            'm + m * m + m ...' may be used, where the 'm' will be filled in
            with the models in order.
        :type operation: 
            String or a sequence of strings
            
        :param parnames:
            Assigns new names for parameters based on their default names. The
            default parameter names are of the form 'A0' and 'A1' where A is the
            parameter name and the number is the sequential number of the model
            with that parameter.
        :type parnames: 
            Dictionary map of default parameter name string to new name string
        
        :param autoshorten: 
            If True, the numerical suffix for paraemeter names will be removed
            if there is only one of that parameter (parnames overrides this)
        :type autoshorten: bool
        
        Any additional arguments should be "parname=parval" form, setting the 
        initial values for the parameters.
        """
        from inspect import isclass
        from collections import defaultdict
        
        self.__dict__['_pars'] = tuple() #necessary for later when getattr is invoked
        
        mods = [get_model_instance(m) for m in models]
        self._models = tuple(mods)
        
        if isinstance(operation,basestring):
            if len(operation.strip()) == 1:
                operation = [operation for i in range(len(mods)-1)]
            else:
                operation = operation.split('m')[1:-1]
        elif len(operation) != len(mods)-1:
            raise ValueError('impossible number of operations')
        
        
        self._ops = ops = tuple(operation)
            
        oplst = ['mval[%i]%s'%t for t in enumerate(ops)]
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
            argnames = [argmap[a][1] for a in cargs]
            for i,a in enumerate(cargs):
                if a not in parnames.values():
                    if argnames.count(argnames[i])==1:
                        argmap[argnames[i]] = argmap[a]
                        del argmap[a]
                        args[i] = argnames[i]
            
        self._pars = tuple(args)
        self._argmap = argmap
        
        #finally apply kwargs
        for k,v in parvals.iteritems():
            setattr(self,k,v)
        
        #now generate the data structures to map CompositeModel par numbers
        #to single model par numbers -- used for quick function evaluation
        parlistpairs = []
        parlistd = defaultdict(dict)
        for i,p in enumerate(self._pars):
            modi,n = argmap[p]
            indinmodel = list(mods[modi].params).index(n)
            parlistpairs.append((modi,indinmodel))
            parlistd[modi][indinmodel] = mods[modi].parvals[indinmodel] #value not important but we give it current value just for checking
            
        parlists = []
        for i in range(len(parlistd.keys())):
            d = parlistd[i]
            innerlist = [d[j] for j in range(len(d.keys()))]
            parlists.append(innerlist)
        
        self._parlistmaps = tuple(parlistpairs)   
        self._parlists = parlists
        
        
    
    def __getattr__(self,name):
        if name in self._pars:
            i,n = self._argmap[name]
            return getattr(self._models[i],n)
        raise AttributeError("'%s' has no attribute '%s'"%(self,name))
    
    def __setattr__(self,name,val):
        if name in self._pars:
            i,n = self._argmap[name]
            setattr(self._models[i],n,val)
        else:
            super(CompositeModel,self).__setattr__(name,val)
    
    @property
    def models(self):
        return self._models
    
    @property
    def ops(self):
        return self._ops  
    
    def f(self,x,*args):
        """
        The composite model function. 
        """
        #TODO: switch back to this system if the mapping technique is too confusing
#        for p,a in zip(self.params,args):
#            setattr(self,p,a)
#        mval = [m.f(x,*m.parvals) for m in self._models]
        if len(args)!=len(self._pars):
            raise ValueError('incorrect number of parameters')
        
        parlists = self._parlists
        for a,(i,j) in zip(args,self._parlistmaps):
            parlists[i][j] = a
        
        #mval = [m.f(x,*parlists[i]) for i,m in enumerate(self._models)]
        #TODO:speed-up/cacheing?
        for m,pl in zip(self._models,parlists):
            m.parvals = pl
        mval = [m(x) for m in self._models]
        return eval(self._opstr)
        
    
    def fitDataFixed(self,*args,**kwargs):
        """
        Fits data using :meth:`FunctionModel.fitData`, but allows this
        :class:`CompositeModel` to hold all parameters of one of the sub-models
        fixed instead of fixing a list of parameters.
        
        The provided arguments will be passed into
        :meth:`FunctionModel.fitData`, except for `fixedpars`. Instead the
        `fixedmods` or `freemods` kwargs are used to determine which parameters
        should be fixed. Either can be specified, but not both.
        
        :param freemods: 
            The indecies of the models (0-based) for which the parameters should
            be treated as fitting parameters.  All other models' parameters will
            be held fixed.
        :type freemods: sequence of ints
        :param fixedmods:
            The indecies of the models (0-based) for which the parameters should
            be held constant. All other models' parameters will be treated as
            free fitting parameters.
        :type fixedmods: sequence of ints
        
        :returns: same return value as :meth:`FunctionModel.fitData`
        
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
    
    
                

class FunctionModel1D(FunctionModel):
    
    """
    This class is the base for 1-dimensional models that are implemented
    as python functions.
    
    **Subclassing**
    The following method *must* be overridden in a subclass:
    
    * f(self,x,...)
    
    The following methods may be overridden for speed of some operations - pars
    should be accessed with self.pardict, self.parvals, or by properties/name:
    
    * integrate(self,lower,upper)
        integrate the model from lower to upper
    * derivative(self,x,dx)
        derivative of the model at x, optinally using spacing dx
    * inv(yval,*args,**kwargs)
        inverse of the model at the requested yvalue
    
    The following attributes may be set for additional information:
    
    * :attr:`xaxisname`
        name of the input axis for this model
    * :attr:`yaxisname`
        name of the output axis for this model
    * :attr:`rangehint`
        a hint for the relevant range for this model as a (lower,upper) tuple
        (or None if no particular range is relevant)
        
    """    
    defaultIntMethod = 'quad'
    defaultInvMethod = 'brentq'
    xaxisname = None
    """
    Name of the input axis for this model.
    """
    yaxisname = None
    """
    Name of the output axis for this model.
    """
    rangehint = None
    """
    A hint for the relevant range for this model as a (lower,upper) tuple (or
    None if no particular range is relevant)
    """
    
    def __call__(self,x):
        """
        call the function on the input x,which will be flattened.  If not 1D,
        the output will be reshaped to the input
        """
        arr = np.array(x,copy=False,dtype=float)
        res = self._filterfunc(arr.ravel(),*self.parvals)
        return res.reshape(arr.shape)
    
    def inv(self,yval,*args,**kwargs):
        """
        Find the x value matching the requested y-value. The inverse is computed
        using root finders from the :mod:`scipy.optimize` module.
        
        :param yval: the output y-value at which to compute the inverse
        :type yval: float
        
        Other args and kwargs are those appropriate for the chosen root-finder,
        except for the keyword `method` which can be a name of any of the root
        finders from :mod:`scipy.optimize`. `method` can also be a function that
        should take f(g(x),*args,**kwargs) and return the x value at which g(x)
        is 0.
        
        The default method depends on the input arguments as follows:
        
        * inv(yval)
            Uses :func:`scipy.optimize.newton`
        * inv(yval,x0)
            Uses :func:`scipy.optimize.newton`, starting the search at x0
        * inv(yval,a,b)
            Uses :func:`scipy.optimize.brentq`, searching the bracketing
            interval [a,b] for the lower and upper edges of the search range.
            
        :returns: the x-value at which the model equals the given `yval`
        
        **Examples**
        
        These examples use Newton's, Brent's, and Ridder's method, respectively.
        
        .. testsetup::
            
            from astropysics.models.builtins import QuadraticModel
        
        .. doctest::
        
            >>> m = QuadraticModel()
            >>> '%.2f'%m.inv(2)
            '1.41'
            >>> '%.2f'%m.inv(9,2,4)
            '3.00'
            >>> '%.2f'%m.inv(16,3,5,method='ridder')
            '4.00'
        
        """    
        import scipy.optimize
        
        if not np.isscalar(yval):
            raise ModelTypeError('generic 1D inverter can only accept scalar inputs')
        
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
    
    def minimize(self,x0,method='fmin',**kwargs):
        """
        Finds a local minimum for the model.
        
        :param x0: The location to start the search
        :type x0: float
        :param method: 
            Can be 'fmin' or 'fmin_powell', to use `scipy.optimize.fmin` and
            `scipy.optimize.fmin_powell`.
        :type method: string
        
        kwargs are passed into the `method` function
        
        :returns: a x value where the model is a local minimum
        """
        return self._optimize(x0,'min',method,**kwargs)
    
    def maximize(self,x0,method='fmin',**kwargs):
        """
        Finds a local maximum for the model.
        
        :param x0: The location to start the search
        :type x0: float
        :param method: 
            Can be 'fmin' or 'fmin_powell', to use `scipy.optimize.fmin` and
            `scipy.optimize.fmin_powell`.
        :type method: string
        
        kwargs are passed into the `method` function
        
        :returns: a x value where the model is a local maximum
        """
        return self._optimize(x0,'max',method,**kwargs)
    
    def findroot(self,x0,method='fmin',**kwargs):
        """
        Finds a root for the model (i.e. location where the model is 0)
        
        :param x0: The location to start the search
        :type x0: float
        :param method: 
            Can be 'fmin' or 'fmin_powell', to use `scipy.optimize.fmin` and
            `scipy.optimize.fmin_powell`.
        :type method: string
        
        kwargs are passed into the `method` function
        
        :returns: the x value where the model is 0
        """
        return self._optimize(x0,'root',method,**kwargs)
    
    def findval(self,val,x0,method='fmin',**kwargs):
        """
        Finds where the model is equal to the value given by the val argument
        
        x0 is the location to start the search
        method can be 'fmin' or 'fmin_powell' (from scipy.optimize)
        kwargs are passed into the scipy.optimize function
        """
        kwargs['valtofind'] = val
        return self._optimize(x0,'val',method,**kwargs)
        
    def _optimize(self,x0,type,method,**kwargs):
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
        elif type == 'val':
            val = kwargs.pop('valtofind')
            g=lambda *args,**kwargs:np.abs(self.f(*args,**kwargs)-val)
        elif type == 'saddle':
            raise NotImplementedError
        else:
            raise ValueError('Unrecognized optimization type')
        
        if method == 'fmin':
            res = fmin(g,x0,tuple(self.parvals),**kwargs)
        elif method == 'fmin_powell':
            res = fmin_powell(g,x0,tuple(self.parvals),**kwargs)
        else:
            raise ValueError('Unrecognized method')
        
        self.lastOpt = res
        return res[0]
    
    def _getInvertedRangehint(self):
        """
        computes the range hint taking into account input filters
        """
        rh = self.rangehint
        
        if rh is None:
            return None
        else:
            #invert the rangehint if the input call type is a recognized string
            callin = self.getCall()[1] if self.getCall() is not None else None
            if callable(callin):
                from warnings import warn
                warn('Could not infer inverse of input filter, so assumed range may be incorrect')
            elif callin == 'log':
                rh = 10**np.array(rh)
            elif callin == 'ln':
                rh = np.exp(rh)
            elif callin == 'pow':
                rh = np.log10(rh)
            elif callin == 'exp':
                rh = np.log(rh)
            
            return rh
        
    def plot(self,lower=None,upper=None,n=100,clf=True,logplot='',data='auto',
                  *args,**kwargs):
        """
        Plot the model function and possibly data and error bars with
        :mod:`matplotlib.pyplot`. The plot will reflect any changes applied with
        :meth:`setCall`.
        
        :param lower: 
            The starting x value for the plot. If None, the bound will be
            inferred from the `data` argument to this function or the
            :attr:`rangehint` attribute of the model.
        :type lower: scalar or None
        :param upper: 
            The ending x value for the plot. If None, the bound will be inferred
            from the `data` argument to this function or the :attr:`rangehint`
            attribute of the model.
        :type upper: scalar or None
        :param n: The number of samples for the plot
        :type n: int
        :param clf: 
            If True, the figure will be cleared before the plot is drawn. 
        :type clf: boolean
        :param logplot: Sets which axes are logarithmic
        :type logplot: '','x','y', or 'xy' string
        :param data: 
            Determines what (if any) data to display. If None, no data is
            displayed. If 'auto', data from the :attr:`data` attribute of the
            :class:`FunctionModel1D` will be used if present, or nothing if the
            :attr:`data` attribute is empty or None. Otherwise, the data to be
            plotted (and possibly errors) can be provided as arrays or as a
            dictionary of inputs to the :func:`matplotlib.pyplot.scatter`
            function.
        :type data: 
            'auto', None, or array like of the form (x,y) or (x,y,(xerr,yerr))
            
        """
        from matplotlib import pyplot as plt
        from operator import isMappingType
        
        isinter = plt.isinteractive()
        try:
            plt.ioff()
            if data is 'auto':
                if self.data:
                    data = self.data
                else:
                    data = None
                    
            if data is not None:
                if lower is None:
                    lower=np.min(data[0])
                if upper is None:
                    upper=np.max(data[0])
            
            if lower is None or upper is None:
                rh = self._getInvertedRangehint()
                if rh is None:
                    raise ValueError("Can't choose limits for plotting without data or a range hint")
                
                
                    
                if lower is None:
                    lower = rh[0]
                if upper is None:
                    upper = rh[1]
            
            if 'x' in logplot:
                x = np.logspace(np.log10(lower),np.log10(upper),n)
            else:
                x = np.linspace(lower,upper,n)
            
            y = self(x)
            
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
            
            if data is not None:
                if isMappingType(data):
                    if clf and 'c' not in data:
                        data['c']='g'
                    plt.scatter(**data)
                else:
                    plt.scatter(data[0],data[1],c='r',zorder=10)
                    if data[2] is not None:
                        if len(data[2].shape)>1:
                            plt.errorbar(data[0],data[1],data[2][1],data[2][0],fmt=None,ecolor='k')
                        else:
                            plt.errorbar(data[0],data[1],data[2],fmt=None,ecolor='k')
            plt.xlim(lower,upper)
            
            if self.xaxisname:
                plt.xlabel(self.xaxisname)
            if self.yaxisname:
                plt.ylabel(self.yaxisname)
                
            if isinter:
                    plt.draw()
                    plt.show()
        finally:
            plt.interactive(isinter)
    
    #Can Override:
    def integrate(self,lower,upper,method=None,n=100,jac=None,**kwargs):
        """
        Numerically compute the definite integral of the model using
        :mod:`scipy.integrate` functions. The integral computed is:
        
        .. math::
            \\int_l^u \\! f(x) \\, j(x) \\, dx
        
        where :math:`j(x)` is the jacobian set from the `jac` argument.
        
        
        
        :param lower: the lower limit of the integral
        :type lower: float
        :param upper: the upper limit of the integral
        :type upper: float
        :param method: 
            The name of a function from :mod:`scipy.integrate`, or if None,
            the class attribute :attr:`defaultIntMethod` will be used.
        :type method: string or None
        :param n: 
            Only has an effect for integration techniques that use samples. If
            an integer, it specifies the number of evenly spaced samples. If it
            as a sequence, it is an array of samples to use for the integral
            (this renders `lower` and `upper` meaningless and their values have
            no effect)
        :type n: int
        :param jac: The jacobian factor to include in the integrand.
        :typ jac: a callable f(x,*params) or None
            
        :returns: The value of the computed definite integral.
        

        .. note::
            Integration methods will store their full output to the attribute
            :attr:`lastintegrate` upon completion.
            
        .. note::
            If overridden in a subclass, the signature should be
            integrate(self,lower,upper,*args,**kwargs), but the args and kwargs
            will typically be ignored where an analytic solution is available.
            
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
        
        
        self.lastintegrate = res
        return res if np.isscalar(res) else res[0]
    
    def integrateCircular(self,lower,upper,*args,**kwargs):
        """
        This is a convinience for self.integrate with jacobian set for a
        azimuthally symmetric 2D radial profile.  
        
        .. math::
            \\int_l^u \\! f(x) \\, 2 \\pi x \\, dx
        
        If a `jac` keyword is provided, it is taken as an additional factor to
        multiply the circular jacobian.
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*2.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*2.0*pi
        return self.integrate(lower,upper,*args,**kwargs)
    
    def integrateSpherical(self,lower,upper,*args,**kwargs):
        """
        This is a convinience for self.integrate with jacobian set for a
        spherically symmetric 3D radial profile.  
        
        .. math::
            \\int_l^u \\! f(x) \\, 4 \\pi x^2 \\, dx
        
        If a `jac` keyword is provided, it is taken as an additional factor to
        multiply the circular jacobian.
        """
        if 'jac' in kwargs:
            kwargs['jac'] = lambda x,*params:kwargs['jac'](x,*params)*x*x*4.0*pi
        else:
            kwargs['jac'] = lambda x,*params:x*x*4.0*pi
        return self.integrate(lower,upper,*args,**kwargs)
        
    def derivative(self,x,dx=None):
        """
        The numerically estimated derivative at x of the form:
        
        .. math::
            \\frac{df}{dx} \\approx \\frac{f(x+\\Delta x)-f(x)}{\\Delta x}
        
        :param x: the value at which to compute the derivative
        :type x: float or array-like
        :param dx: 
            The spacing to assume for the numerically-computed derivative
            (:math:`\\Delta x` above). If None, and an array is provided for
            `x`, spacing will be inferred from the spacing between the `x`
            values. Otherwise, it defaults to 1.
        
        .. note::
            If overridden in a subclass, the signature should be
            derivative(self,x,dx=1) , but dx will typically be ignored if an
            analytic solution is available. 
        """
        if dx is None:
            x = np.array(x,copy=False)
            if len(x.shape)>0:
                dx = np.convolve(x,[1,-1],mode='valid')
                dx = np.insert(dx,0,np.mean(dx))/2
            else:
                dx = 1
        return (self(x+dx)-self(x))/dx
    
    def pixelize(self,xorxl,xu=None,n=None,edge=False,sampling=None):
        """
        This method integrates over the model for a number of ranges
        to get a 1D "pixelized" version of the model.  
        
        :param xoroxl: 
            If array, specifies the location of each of the pixels. If float,
            specifies the lower edge of the pixelized section.
        :type xorxl: array-like or float
        
        :param xu: 
            Specifies the upper edge of the pixelized section. Ignored if
            `xorxl` is array-like.
        :type xu: float
        :param n:
            Specifies the number of pixels. Ignored if `xorxl` is array-like.
        :type n: int
        :param edge: 
            If True, pixel locations are for the edges (xorxl lower edges, xu
            upper edge), or if False, they are pixel centers.
        :type edge: bool  
        :param sampling: 
            If None, each pixel will be computed by integrating. Otherwise, 
            the number of samples to use in each pixel.
        :type sampling: int or None
        
        :returns: 
            Integrated values for the function as an array of size `n` or
            matching `xorxl`.
        
        """
        if xu is None and n is not None or xu is not None and n is None:
            raise TypeError('must specify both xu and n')
        
        if xu is None and n is None:
            xorxl = np.array(xorxl)
            
            if not edge:
                dx = (xorxl[-1]-xorxl[0])/(len(xorxl)-1)
                edges = np.concatenate((xorxl-dx/2,(xorxl[-1]+dx/2,)))
                n = xorxl.size
            else:
                edges = xorxl
                n = xorxl.size-1
        else:
            if edge:
                edges = np.linspace(xorxl,xu,n+1)
            else:
                dx = (xu-xorxl)/(n-1)
                edges = np.linspace(xorxl,xu+dx,n+1)-dx/2
        
        if sampling is None:
            dx = np.convolve(edges,[1,-1],mode='valid')
            arr = [self.integrate(edges[i],edges[i+1]) for i in range(n)]
            arr = np.array(arr)/dx
        else:
            arr = [np.mean(self(np.linspace(edges[i],edges[i+1],sampling))) for i in range(n)]
            arr = np.array(arr)
            
        return arr
                
    
    #TODO:move this up into FunctionModel
    def setCall(self,calltype=None,xtrans=None,ytrans=None,**kwargs):
        """
        Sets the type of function evaluation to occur when the model is called.
        Changes the output of a function call to be :math:`t_y(g(f(t_x(x))))`,
        where :math:`g(x)` is set by `calltype`, :math:`t_x` is given by the
        `xtrans` argument, and :math:`t_y` is given by `ytrans`.
        
        :param calltype:  
            Specifies what should be output when the model is called. Can be:
        
            * None
                basic function evaluation
            * 'derivative'
                derivative at the location (see :meth:`derivative`)
            * 'integrate'
                integral - specify `upper` or `lower` kwargs and the evaluation
                location will be treated as the other bound. If neither is
                given, lower = 0 is assumed. (see :meth:`integrate`)
            * 'integrateCircular'
                Polar integration -- see :meth:`integrateCircular`
            * 'integrateSpherical'
                Spherical integration -- see :meth:`integrateSpherical`
            * Any other string that is the name of a method on this object.  It
              will be called with the function value as its first argument.
              
        :param xtrans: 
            Transformations applied to the input value before being passed into
            the model.  See below for valid forms.
        :type xtrans: string or None
        :param ytrans:
            Transformations applied to the output value. See below for valid
            forms.
        :type ytrans: string or None
        
        Any kwargs are passed into the function specified in `calltype`.
        
        `xtrans` and `ytrans` transformation functions can accept the following
        values:
            
            * None
                :math:`t(x) = x`
            * 'log'
                :math:`t(x) = \\log_{10}(x)`
            * 'ln'
                :math:`t(x) = \\ln(x)`
            * 'log##.#'
                :math:`t(x) = \\log_{\\#\\#.\\#}(x)`
            * 'pow'
                :math:`t(x) = 10^x`
            * 'exp'
                :math:`t(x) = e^x`
            * 'pow##.#'
                :math:`t(x) = (\\#\\#.\\#)^x`
        
        .. warning:: 
            There may be unintended consequences of this method due to methods
            using the call value instead of the default function evaluation
            result. You have been warned...
        
        """
        from types import MethodType
        import inspect
        from functools import partial
        
        transmap={'log':np.log10,'ln':np.log,'pow':partial(np.power,10),'exp':np.exp}
        xts,yts = xtrans,ytrans
        if isinstance(xtrans,basestring):
            if 'log' in xtrans and xtrans!='log':
                basex = float(xtrans.replace('log',''))
                if basex==10:
                    xtrans = transmap['log']
                else:
                    sclx = np.log(basex)
                    xtrans = lambda v:np.log(v)/sclx
            elif 'pow' in xtrans and xtrans!='pow':
                basex = float(xtrans.replace('log',''))
                xtrans = lambda v:basex**np.log(v)
            else:
                xtrans = transmap[xtrans]
        if isinstance(ytrans,basestring):
            if 'log' in ytrans and ytrans!='log':
                basey = float(ytrans.replace('log',''))
                scly = np.log(basey)
                if basey==10:
                    ytrans = transmap['log']
                else:
                    scly = np.log(base)
                    ytrans = lambda v:np.log(v)/scly
            elif 'pow' in ytrans and ytrans!='pow':
                basey = float(ytrans.replace('log',''))
                ytrans = lambda v:basey**np.log(v)
            else:
                ytrans = transmap[ytrans]
            
        if calltype is None:
            if xtrans and ytrans: 
                self._filterfunc = MethodType(lambda self,x,*pars:ytrans(self.f(xtrans(x),*pars)),self,self.__class__)
            elif xtrans:
                self._filterfunc = MethodType(lambda self,x,*pars:self.f(xtrans(x),*pars),self,self.__class__)
            elif ytrans:
                self._filterfunc = MethodType(lambda self,x,*pars:ytrans(self.f(x,*pars)),self,self.__class__)
            else:
                self._filterfunc = self.f
        else:
            try:
                newf = getattr(self,calltype)
                
                if not callable(newf):
                    raise AttributeError
            except AttributeError:
                raise AttributeError('function %s not found in %s'%(calltype,self))
                
            if 'integrate' in calltype:
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
                
            self._filterfunc = MethodType(callfunc,self,self.__class__)
            
        self._filterfunctype = (calltype,xts,yts)
        
    def getCall(self):
        """
        Retreives infromation about the calling function.
        
        :returns:
            The type of evaluation to perform when this model is called - a
            string like that of the type passed into :meth:`setCall`, or None if
            the model function itself is to be called.
        """
        if hasattr(self,'_filterfunctype'):
            return self._filterfunctype
        else:
            return None
        

class FunctionModel1DAuto(FunctionModel1D):
    """
    A :class:`FunctionModel1D` that infers the parameters from the :meth:`f`
    function.
    
    Equivalent to simply setting the metaclass to :class:`AutoParamsMeta`.
    """
    __metaclass__ = AutoParamsMeta

class CompositeModel1D(FunctionModel1D,CompositeModel):
    """
    This model is a composite model of FunctionModel1D models with a few 
    extra additions for 1D models.
    """
    def __init__(self,*args,**kwargs):
        """
        see `CompositeModel.__init__` for arguments
        """
        super(CompositeModel1D,self).__init__(*args,**kwargs)
        for m in self._models:
            if not isinstance(m,FunctionModel1D):
                raise ModelTypeError('Input model %s is not a 1D model'%m)
        self._filters = None
    #__init__.__doc__ = CompositeModel.__init__.__doc__
    
    def f(self,x,*pars):
        res = CompositeModel.f(self,x,*pars)
        
        if self._filters is None:
            return res
        else:
            for filter in self._filters:
                res = filter(res)
            return res
        
    @property
    def rangehint(self):
        rhints = []
        for m in self._models:
            if not hasattr(m,'rangehint'):
                return None
            else:
                rhint = m.rangehint
                if rhint is None:
                    return None
                else:
                    rhints.append(rhint)
        rhints = np.array(rhints)
        return np.min(rhints[:,0]),np.max(rhints[:,1])
#        mx,mi = np.max(rhints,0),np.min(rhints,0)
#        return mi[0],mx[1],mi[2],mx[3]
    
    #TODO: remove
    def addFilter(self,filter):
        """
        :deprecated:
        
        This adds a function to be applied after the model is evaluated
        """
        from warnings import warn
        warn('addFilter/clearFilters is deprecated in favor of setCall',DeprecationWarning)
        
        if self._filters is None:
            self._filters = []
        
        if not callable(filter):
            raise ValueError('input filter is not a function')
        
        self._filters.append(filter)
        
    def clearFilters(self):
        """
        :deprecated:
        
        this clears all previously added filters
        """
        from warnings import warn
        warn('addFilter/clearFilters is deprecated in favor of setCall',DeprecationWarning)

        self._filters = None
        
    def addLowerBoundFilter(self,bound):
        """
        :deprecated:
        """
        def bndfunc(x):
            x[x<bound] = bound
            return x
        self.addFilter(bndfunc)
        
class DatacentricModel1D(FunctionModel1D):
    """
    A FunctionModel1D that *requires* data to compute its value.
    """
    def __call__(self,x):
        if self.data is None:
            raise ValueError('DataModel1D must have data to execute')
        arr = np.array(x,copy=False,dtype=float)
        res = self._filterfunc(arr.ravel(),*self.parvals)
        return res.reshape(arr.shape)
    

class DatacentricModel1DAuto(DatacentricModel1D):
    """
    A :class:`DatacentricModel1D` that infers the parameters from the :meth:`f`
    function.
    
    Equivalent to simply setting the metaclass to :class:`AutoParamsMeta`.
    """
    __metaclass__ = AutoParamsMeta
        

class ModelSequence(object):
    """
    A group of models with attached parameters that can be used to infer the
    parameter values at an arbitrary point in the space defined by the models.  
    """
    def __init__(self,models,extraparams=None,outputcontraction=None,
                      interpolation='linear',interpolationdirection='y',
                      offgrid=None):
        """
        `models` must be a sequence of models or a 2-tuple with the second 
        element a dictionary.  If a model sequnce, they must all have the same 
        parameter names and have compatible inputs.  Otherwise, the first 
        element specifies the type of model to use, and the second is a 
        dictionary mapping parameter names to a sequence of parameter values.

        The resultant models are expected to be a "sequence" in the sense that 
        they can be interpolated across to get a meaningful parameter value 
        between two models (e.g. double-valued model grids are not meaningful)
        
        `extraparams` are parameters not a part of the model that are still 
        available to be interpolated over.  It must be None or a dictionary with 
        values that are the same length as the model sequence.
        
        `output contraction` can be either a function to call on the output of 
        the models returning a scalar, None to do nothing to the output, 'dist' 
        to use the euclidian distance (e.g. sqrt(sum(output**2)), or 'sqdist' to 
        use the square of the distance.
        
        `interpolation` may be either 'linear' for simple linear interpolation 
        along the y-axis, or a 1D Model that will be fit with fitData on the 
        contracted outputs.
        
        `interpolationdirection` determines the direction along which 
        interpolation is to be performed, and can be:
        
        * 'y': interpolate along the output axis (by far the fastest method)
        * 'x': interpolate along the input axis - requires that inv be defined
               for the models
        * 'perp': interpolate along the perpendicular of the tangent line for 
                  the closest model to the data point.  Only meaningful for 
                  FunctionModel1D Models.
                  
        `offgrid` determines the behavior if a location is requested that is 
        beyond the limits of the model sequence.  It can be:
        
        * 'raise': raises a ValueError
        * 'warn': issue a warning
        * None: just output whatever the interpolating function provides - for 
          linear interpolation, this takes the closest value.
          
        """        
        from operator import isMappingType
        
        if len(models)==2 and isMappingType(models[1]):
            modtype = get_model(models[0])
            params = np.array(models[1].values())
            params = [dict([(models[1].keys()[i],v)] for v in t) for i,t in enumerate(params.T)]
            models = [modtype(**params[i]) for m in range(len(params))]
        
        params = None
        
        for m in models:
            if params is None:
                params = m.params
            else:
                if m.params != params:
                    raise ValueError('model %s does not match parameters for other models'%m)
        
        if extraparams is not None:
            self._extraparams = {}
            for n,ps in extraparams.iteritems():
                arr = np.array(ps)
                if extraparams[n].size != len(models):
                    raise ValueError('too many/few extra parameters for parameter %s'%n)
                self._extraparams[n] = arr
        else:
            self._extraparams = None
        
        self._params = params
        self._models = tuple(models)
        self._extraparams = extraparams
        
        self.outputcontraction = outputcontraction
        self.interpolation = interpolation
        self.interpolationdirection = interpolationdirection
        self.offgrid = offgrid
        
    
    def _getOutputcontraction(self):
        return self._outcont
    def _setOutputcontraction(self,val):
        if val is None:
            val = lambda x:x
        elif isinstance(val,basestring):
            if val == 'dist':
                val = lambda x:np.sum(x**2)**0.5
            elif val == 'sqdist':
                val = lambda x:np.sum(x**2)
            else:
                raise ValueError('invalid string for outputcontraction')
        elif not callable(val):
            raise TypeError('invalid type for outputcontraction')
        self._outcont = val
    outputcontraction = property(_getOutputcontraction,_setOutputcontraction,doc=None)
    
    def _getInterpolation(self):
        return self._interp
    def _setInterpolation(self,val):
        if isinstance(val,basestring):
            if val == 'linear' or val == 'linearinterp':
                self._interp = 'linear'
            else:
                raise ValueError('invalid interpolation string value %s'%val)
        elif isinstance(val,ParametricModel):
            self._interp = val
        else:
            raise TypeError('invalid type for interpolation')
    interpolation = property(_getInterpolation,_setInterpolation,doc=None)
    
    def _getInterpolationdirection(self):
        return self._interpdir
    def _setInterpolationdirection(self,val):
        if val not in ('x','y','perp'):
            raise ValueError('invalid interpolation direction %s requested'%val)
        self._interpdir = val
    interpolationdirection = property(_getInterpolationdirection,_setInterpolationdirection,doc=None)
    
    def _getOffgrid(self):
        return self._offgrid
    def _setOffgrid(self,val):
        if val not in ('raise','warn',None):
            raise ValueError("invalid value for offgrid - must be 'raise', 'warn', or None")
        self._offgrid = val
    offgrid = property(_getOffgrid,_setOffgrid,doc=None)
    
    
    @property
    def params(self):
        """
        a tuple of the possible parameter names for this ModelSequence
        """
        return tuple(self._params)
    
    @property
    def extraparams(self):
        """
        A tuple of the possible extra parameter names for this ModelSequence or
        None if there are none
        """
        return None if self._extraparams is None else tuple(self._extraparams) 
    
    def getParamArray(self,parnames):
        """
        Generates an array of values for the requested parameter corresponding to
        each of the models.
        
        `parnames` can be a single parameter name (output is then a scalar), a
        stirng with a comma-seperated list of parameters, a sequence of 
        parameter names, or None to get a dictionary mapping parameter names
        to their arrays
        """    
        scalarout = dictout = False
        if parnames is None:
            parnames = self._models[0].params
            if self._extraparams is not None:
                parnames.extend(self.extraparams)
            dictout = True
        elif isinstance(parnames,basestring):
            parnames = parnames.split(',')
            if len(parnames) == 1:
                scalarout = True
            
        res = []
        for p in parnames:
            if p in self.params:
                res.append(np.array([getattr(m,p) for m in self._models]))
            elif p in self._extraparams:
                res.append(np.array(self._extraparams[p]))
            else:
                raise ValueError('invalid parameter name %s requested'%p)
            
        if dictout:
            return dict([t for t in zip(parnames,res)])
        elif scalarout:
            return res[0]
        else:
            return res
        
    def getParam(self,x,y,parnames=None,contracty=True):
        """
        Computes the value of a requested parameter at the provided point in 
        the input/output space of the models.
        
        `parnames` can be a single parameter name (output is then a scalar), a
        stirng with a comma-seperated list of parameters, a sequence of 
        parameter names, or None to get a dictionary mapping parameter names
        to their value at the provided point
        
        `x` and `y` should be inputs and outputs, respectively, of the Model
        objects that define this sequence.
        
        if `contracty` is True, the output contraction will be applied to the 
        provided y-value.  Otherwise, it should be a scalar.
        """
        scalarout = dictout = False
        if parnames is None:
            parnames = self._models[0].params
            if self._extraparams is not None:
                parnames.extend(self.extraparams)
            dictout = True
        elif isinstance(parnames,basestring):
            parnames = parnames.split(',')
            if len(parnames) == 1:
                scalarout = True
            
        for p in parnames:
            if p not in self._params and p not in self._extraparams:
                raise ValueError('invalid parameter %s'%p)
            
        if contracty:
            y = self._outcont(y)
        
        if self._interpdir == 'y':
            interpin = np.array([self._outcont(m(x)) for m in self._models])
            interpat = y
        elif self._interpdir == 'x':
            interpin = np.array([m.inv(y) for m in self._models])
            interpat = x
        elif self._interpdir == 'perp':
            raise NotImplementedError('perp interp direction not working yet')
            from .builtins import LinearModel
            ys = np.array([m(x)-y for m in self._models])
            closestmodel = self._models[np.where(ys==np.min(ys))[0][0]]
            dmod = closestmodel.derivative(x)
            linmod = LinearModel(m=-1/dmod,b=0)
            linmod.b = y - linmod(x)
            interpin = (1 + linmod.m)**0.5*np.array([intersect_models(m,linmod,bounds=bnds(x-dmod,x+dmod)) for m in self._models]).ravel()
            interpat = x
        else:
            raise RuntimeError('This point should never be reachable in getParam!')
        
        #take parnames sequence and for each one add the result to the res list
        res = []
        sorti = np.argsort(interpin)
        interpin = interpin[sorti]
        if self._offgrid is not None:
            if interpat <  interpin[0] or interpat > interpin[-1]:
                if self._offgrid == 'warn':
                    from warnings import warn
                    warn('off grid in getParam input %s'%interpat)
                elif self._offgrid == 'raise':
                    raise ValueError('off grid in getParam input %s'%interpat)
                else:
                    raise RuntimeError('This point should never be reachable in getParam!')
        for p in parnames:
            if p in self._params:
                interpout = np.array([getattr(m,p) for m in self._models])[sorti]
            else: #instead in _extraparams
                interpout = self._extraparams[p][sorti]
                
            if self._interp == 'linear':
                res.append(np.interp(interpat,interpin,interpout))
            else:
                self._interp.fitData(interpin,interpout)
                res.append(self._interp(interpat))
            
        if scalarout:
            return res[0]
        elif dictout:
            return dict([t for t in zip(parnames,res)])
        else:
            return res
        
    def getParams(self,xs,ys,parnames=None,contracty=True):
        """
        Get parameters for an array of  inputs - `xs` and `ys` should be 
        matching sequences which each item matches the input and ouptut types of 
        the models.  For more info, see ModelSequence.getParam
        """
        if len(xs) != len(ys):
            raise ValueError("xs and ys don't match") 
        return np.array([self.getParam(x,y,parnames,contracty) for x,y in zip(xs,ys)])
        
    def plot1D(self,x1=None,x2=None,legend=False,clf=True,n=100,**kwargs):
        """
        Plots the models in this ModelSequence assuming the models are 1D
        
        `x1` and `x2` are the upper and lower limits to the plot.  If not 
        specified, the rangehints for the models will be used
        
        The output contraction function will be used to compute the y-axis.
        
        `legend` determines if the legend with parameter value should be shown
        on the plot
        
        kwargs are passed into the models' plot method (except for n)
        """
        from matplotlib import pyplot as plt
        
        isinter = plt.isinteractive()
        try:
            plt.ioff()
            if clf:
                plt.clf()
                
              
            if x1 is None:
                x1 = [m._getInvertedRangehint()[0] for m in self._models if m._getInvertedRangehint() is not None]
                if len(x1)>0:
                    x1 = np.min(x1)
                else:
                    raise ValueError('could not infer range from models - must provide x1')
            if x2 is None:  
                x2 = [m._getInvertedRangehint()[1] for m in self._models if m._getInvertedRangehint() is not None]
                if len(x2)>0:
                    x2 = np.max(x2)
                else:
                    raise ValueError('could not infer range from models - must provide x2')  
            kwargs['n'] = n
            kwargs['clf'] = False
            
            #allow legend to be a string 'par1,par2' specifying which parameters to show
            if legend and isinstance(legend,basestring):
                labelpars = legend.split(',')
            else:
                labelpars = list(m.pardict)
                if self._extraparams is not None:
                    labelpars.extend(self._extraparams)
            
            fixedlabel = kwargs.pop('label',None)        
            for i,m in enumerate(self._models):
                if fixedlabel is None:
                    label = ['%s=%.2g'%t for t in m.pardict.iteritems() if t[0] in labelpars]
                    for ep,pv in self._extraparams.iteritems():
                        if ep in labelpars:
                            label.append('%s=%.2g'%(ep,pv[i]))
                    label = ','.join(label)
                    kwargs['label']=label
                elif fixedlabel is not True:
                    kwargs['label'] = fixedlabel
                    fixedlabel = True
                elif 'label' in kwargs:
                    del kwargs['label']
                m.plot(x1,x2,**kwargs)
            
            if legend:
                plt.legend()
            
            if isinter:        
                plt.show()
                plt.draw()
        finally:
            plt.interactive(isinter)


class InputCoordinateTransformer(object):
    """
    This mixin (i.e. a class intended to be subclassed to provide extra
    functionality) class converts  FunctionModel input values from one 
    coordinate system to another.
    
    In subclasses, the following should be defined at class level:
    
    * '_inputtransforms': a dictionary of dictionaries where the 
      values are functions that define the transforms as
      _intransforms[incoordsys][funccoordsys] = func
    * 'incoordsys': the default input coordinate system
    * 'fcoordsys': the coordinate system the function is defined in 
    
    to implement the transform, just call tranformCoordinates(input) and 
    feed the return value into the f-function
    """
    _inputtransforms = {}
    incoordsys = None
    _fcoordsys = None
    
    def transformCoordinates(self,x,incoordsys=None,outcoordsys=None):
        """
        Transform from the input coordinate system into that defined for the
        model function.
        
        :param incoordsys: 
            The input coordinate system name to use for this model. If None,
            :attr:`incoordsys` will be used.
        :type incoordsys: string or None
        :param outcoordsys: 
            The output coordinate system name to use for this model. If None, 
            the standard for this model will be used.
        :type outcoordsys: string 
        
        if incoordsys is None, self.`incoordsys will be used
        """
        if incoordsys is None:
            incoordsys = self.incoordsys
        if outcoordsys is None:
            outcoordsys = self._fcoordsys
        if incoordsys == outcoordsys:
            return x
        else:
            return self._inputtransforms[incoordsys][outcoordsys](*x)
        
    def addTransform(self,input,output,func):
        """
        Register a function as a transform from one coordinate system to
        another.
        
        :param input: The name for the input system of this function.
        :type input: string 
        :param output: The name for the output system of this function.
        :type output: string
        :param func: 
            The function to perform the transform. It should take one argument,
            an array with at least 2 dimensions where the first dimension is defined
            by the coordinate system.
        :type func: callable
        """
        import inspect
        from collections import defaultdict
        
        try:
            args, varargs, varkw, defaults = inspect.getargspec(func)
            if len(args)-1 > len(defaults) or varkw:
                raise TypeError('input function must take one argument')
        except TypeError:
            raise TypeError('input func is not a callable')
        
        #make sure we are in the instance
        if '_inputtransforms' not in self.__dict__:
            dd = defaultdict(dict)
            dd.update(self._inputtransforms)
            self._inputtransforms = dd
            
        dd[input][output] = func
            

class FunctionModel2DScalar(FunctionModel,InputCoordinateTransformer):
    """
    This class is a FunctionModel that maps a two-dimensional input to a
    scalar output for each 2D coordinate.  The input thus is at least a 
    two-dimensional array, with the first dimension of length 2.
    """
    from ..coords import cartesian_to_polar,polar_to_cartesian
    
    _inputtransforms = {'cartesian':{'polar':cartesian_to_polar},
                        'polar':{'cartesian':polar_to_cartesian}}
    incoordsys = 'cartesian'
    _fcoordsys = 'cartesian'
    
    def _filterfunc(self,x,*params):
        return self.f(self.transformCoordinates(x),*params)
    
    def __call__(self,x):
        """
        call the function on the input x,which will be flattened.  If not 1D,
        the output will be reshaped to the input
        """
        arr = np.array(x,copy=False,dtype=float)
        
        if len(arr.shape) > 1:
            subshape = arr.shape[1:]
        elif len(arr.shape) == 1:
            subshape = tuple()
        else:
            raise ModelTypeError('Scalar cannot be fed into 2D model')
        
        try:
            arr = arr.reshape(2,np.prod(subshape))
        except ValueError:
            raise ModelTypeError('2D model input must have first dimension of length 2')
        
        return self._filterfunc(arr,*self.parvals).reshape(subshape)
    
    def integrateCircular(self,outr,inr=0,theta=(0,2*pi),**kwargs):
        """
        Integrates using a circular jacobian from `inr` to `outr`
        
        :param outr: outer limit of integration
        :type outr: float 
        :param inr: inner limit of integration
        :type inr: float
        :param theta: 
            The angular range of the integral in radians as (lower,upper).
        :type theta: a tuple of floats
        
        kwargs are passed into :func:`scipy.integrate.dblquad`
        """
        import scipy.integrate as itg
        oldcoordsys = self.incoordsys
        try:
            self.incoordsys = 'polar'
            f = lambda th,r:r*self((r,th))
            
            res = itg.dblquad(f,inr,outr,lambda y:theta[0],lambda y:theta[1],**kwargs)
        finally:
            self.incoordsys = oldcoordsys
            
        self.lastintegrate = res
        return res[0]
    
    def integrateCartesian(self,xl,xu,yl,yu,**kwargs):
        """
        Integrates the function in a rectangular area defined in rectangular
        coordinates.
        
        :param xl: left bound of integral
        :type xl: float
        :param xu: right bound of integral
        :type xu: float
        :param yl: lower bound of integral
        :type yl: float
        :param yu: upper bound of integral
        :type yu: float
        
        kwargs are passed into :func:`scipy.integrate.dblquad`
        """
        import scipy.integrate as itg
        
        oldcoordsys = self.incoordsys
        try:
            self.incoordsys = 'cartesian'
            f = lambda y,x:self((x,y))
            
            res = itg.dblquad(f,xl,xu,lambda y:yl,lambda y:yu,**kwargs)
        finally:
            self.incoordsys = oldcoordsys
            
        self.lastintegrate = res
        return res[0]
    
    def gradient(self,x):
        raise NotImplementedError
    
    def pixelize(self,xl,xu,yl,yu,nx=100,ny=100,sampling=None):
        """
        Generates a 2D array of the model smoothed/integrated over a certain
        pixel size.
        
        :param xl: left edge of pixelized region
        :type xl: float
        :param xu: right edge of pixelized region
        :type xu: float
        :param yl: lower edge of pixelized region
        :type yl: float
        :param yu: upper edge of pixelized region
        :type yu: float
        :param nx: number of horizontal pixels
        :type nx: int
        :param ny: number of vertical pixels
        :type ny: int
        :param sampling: 
            If None, each pixel will be computed by integrating the model over
            the area of the pixel. Otherwise, sampling gives the factor to
            multiply nx and ny by to get the total number of "virtual pixels"
            that will be averaged to get the actual result at the given grid
            point. If <=1, this is equivalent to directly sampling the model.
        :type sampling: None or int
            
        :returns: An nx X ny array with the pixelized model.
        
        .. note::
            The attribute :attr:`lastintegrate` stores the result of the
            integrations if sampling is None, with lastintegrate[0] storing the
            result and lasteintegrate[1] storing the error on the integral for
            each pixel.
        """
        import scipy.integrate as itg
        
        oldcoordsys = self.incoordsys
        try:
            self.incoordsys = 'cartesian'        
            
            if sampling is None:
                xs = np.linspace(xl,xu,nx+1)
                xs = np.vstack((xs,np.roll(xs,-1))).T[:-1]
                ys = np.linspace(yl,yu,ny+1)
                ys = np.vstack((ys,np.roll(ys,-1))).T[:-1]
                
                res = np.empty((2,nx,ny))
                f = lambda y,x:self((x,y))
                
                #TODO:optimize more?
                for i,(xl,xu) in enumerate(xs):
                    for j,(yl,yu) in enumerate(ys):
                        res[:,i,j] = itg.dblquad(f,xl,xu,lambda y:yl,lambda y:yu)
                        
                self.lastintegrate = res #res[0]= result, res[1]= integration error
                return res[0]
            else:
                da = ((xu-xl)/nx)*((yu-yl)/ny)
                sampling = int(sampling)
                
                if sampling<=1:
                    xg,yg = np.mgrid[xl:xu:1j*nx,yl:yu:1j*ny]
                    return self((xg,yg))*da
                else:
                    xg,yg = np.mgrid[xl:xu:1j*nx*sampling,yl:yu:1j*ny*sampling]
                    res = self((xg,yg))
                    res = res.reshape((nx*sampling,ny,sampling)).mean(axis=-1)
                    res = res.reshape((nx,sampling,ny)).mean(axis=-2)
                    return res*da
            
        finally:
            self.incoordsys = oldcoordsys
            
    def getFluxSize(self,flux=0.5,frac=True,mode='radial',cen=(0,0),v0=1,
                    minfunc='fmin',intkwargs=None,**kwargs):
        """
        Compute the radius/area enclosing a specified amount of flux.
        
        :param flux:
            Specifies the flux value at which to compute the size (in model
            units if `frac` is False).
        :type flux: float
        :param frac: 
            If True, interprets `flux` as a fraction of the total flux instead
            of an absolute model flux value.
        :type frac: bool
        :param mode: Specifies the way to compute flux.  Can be:
        
            * 'radial' : computes the radius enclosing the specified flux - 
              return value in this case is a single scalar
            * 'square' : computes the square enclosing the specified flux - 
              return value in this case is a single scalar with the box length
            * 'rectangular' : computes the box enclosing the specified flux - 
              return value in this case is a 2-tuple (xsize,ysize)
      
        :param cen: 
            Specifies the center to assume. Currently this must be (0,0) for
            'radial' profiles.
        :type cen: 2-tuple (x,y)
        :param v0: 
            The initial guess at which to start (for methods that require it).
        :type v0: float
        :param minfunc: 
            The name of a function in :mod:`scipy.optimize` that should do the
            minimizing.
        :type minfunc: string
        :param intkwargs: Keyword arguments for the integrate method.
        :type intkwargs: dictionary
        
        kwargs are passed into the minimization function
        
        :returns: 
            The location at which the flux enclosed is given by `flux` as
            specified by the `mode` argument.
        
        """
        import scipy.optimize
        
        fmin = getattr(scipy.optimize,minfunc)
        
        if intkwargs is None:
            intkwargs = {}
        
        if mode == 'radial':
            if cen != (0,0):
                raise NotImplementedError('radial profiles must be centered on (0,0)')
            if frac:
                total = self.integrateCircular(np.inf,**intkwargs)
                flux = flux * total
            def f(r):
                intres = self.integrateCircular(r,**intkwargs)-flux
                return intres*intres
            
            if np.isscalar(v0):
                v0 = (v0,)
        elif mode == 'square':
            x0,y0 = cen
            if frac:
                total = self.integrateCartesian(-np.inf,np.inf,-np.inf,np.inf,**intkwargs)
                flux = flux * total
            def f(l):
                intres = self.integrateCartesian(x0-l,x0+l,x0-l,x0+l,**intkwargs)-flux
                return intres*intres
            
            if np.isscalar(v0):
                v0 = (v0,)
        elif mode == 'rectangular':
            x0,y0 = cen
            if frac:
                total = self.integrateCartesian(-np.inf,np.inf,-np.inf,np.inf,**intkwargs)
                flux = flux * total
            def f(ls):
                lx,ly = ls
                intres = self.integrateCartesian(x0-lx,x0+lx,y0-ly,y0+ly,**intkwargs)-flux
                return intres*intres
            
            if np.isscalar(v0):
                v0 = (v0,v0)
        else:
            raise ValueError('unrecognized mode')
                
        if minfunc!='brent':
            res = fmin(f,v0,full_output=1,**kwargs)
        else:
            res = fmin(f,full_output=1,**kwargs)
        self.lastfluxsize = res
        val = res[0]
        
        return val.ravel()[0] if val.size == 1 else tuple(val)
    
    def plot(self,datarange=None,nx=100,ny=100,clf=True,cb=True,data='auto',
                  log=False,**kwargs):
        """
        Plots the model over the range requested using the
        :func:`matplotlib.pyplot.imshow` function.
        
        :param datarange: 
            Specifies the range to plot in the form (xl,xu,yl,yu). If None, it
            will be inferred from the data.
        :type datarange: 4-tuple or None
        :param nx: determines the number of pixels in the x-direction
        :type nx: int
        :param ny: determines the number of pixels in the y-direction
        :type ny: int
        :param clf: If True, the figure will be cleared before plotting.
        :type clf: bool
        :param data: 
            Data to be plotted along with the model.  If None, no data will be
            plotted, or if 'auto', it will be taken from :attr:`model.data` if 
            present.
        :type data: array-like, None, or 'auto'
        :param log:  
            Can be False, True (natural log), '10' (base-10), 'mag', or 'mag##'
            (pogson magnitudes with ## as zeropoint)
        :type log: bool or string
        
        kwargs are passed into :func:`matplotlib.pyplot.imshow`
        
        :except ValueError: if data is not present and `datarange` is None
        
        """
        from matplotlib import pyplot as plt
        from operator import isMappingType
        
        isinter = plt.isinteractive()
        try:
            plt.ioff()
            if data == 'auto':
                data = self.data or None
                
            maxmind = None
            if data: #TODO:correct coord conv
                xd,yd = data[0][0],data[0][1]
                dataval = data[1]
                if datarange is None:
                    datarange = (np.min(xd),np.max(xd),np.min(yd),np.max(yd))
                maxmind = (np.max(dataval),np.min(dataval))
            elif datarange is None:
                raise ValueError("Can't choose limits for plotting without data or a range hint")
                
            
            grid = np.mgrid[datarange[0]:datarange[1]:1j*nx,datarange[2]:datarange[3]:1j*ny]
            res = self(grid)
            
            if log:
                if 'mag' in log:
                    lognomag = log.replace('mag','')
                    zpt = float(lognomag) if lognomag.strip() != '' else 0
                    logfunc = lambda x:zpt-2.5*np.log10(x)
                elif log == '10':
                    logfunc = np.log10
                else:
                    logfunc = np.log 
                res = logfunc(res)
                if data:
                    dataval = logfunc(dataval)
                if maxmind is not None:
                    maxmind = logfunc(maxmind)
            
            if maxmind:
                norm = plt.normalize(min(np.min(res),maxmind[1]),max(np.max(res),maxmind[0]))
            else:
                norm = plt.normalize(np.min(res),np.max(res))
            
            if clf:
                plt.clf()
                
            kwargs.setdefault('aspect','auto')
            #plt.imshow(res[::-1].T,norm=norm,extent=datarange,origin='lower',**kwargs)
            plt.imshow(res.T,norm=norm,extent=datarange,origin='lower',**kwargs)
            
            if cb:
                if isMappingType(cb):
                    plt.colorbar(**cb)
                else:
                    plt.colorbar()
            
            if data:
                if isMappingType(data):
                    kwscat = dict(data)
                else:
                    kwscat = {}
                kwscat.setdefault('norm',norm)
                kwscat.setdefault('c',dataval)
                plt.scatter(xd,yd,**kwscat)
                    
            plt.xlim(datarange[0],datarange[1])
            plt.ylim(datarange[2],datarange[3])
            
            if isinter:
                plt.draw()
                plt.show()
        finally:
            plt.interactive(isinter)
    def plot3d(self,datarange=None,nx=100,ny=100,clf=True,cb=True,data='auto',**kwargs):
        """
        Generate a 3D plot of the model using mayavi.
        
        See :meth:`plot` for meaning of the arguments
        """
        from enthought.mayavi import mlab as M
        from operator import isMappingType
        
        if data == 'auto':
            if self.data:
                data = self.data[:2]
            else:
                data = None
        
        if data: #TODO:correct coord conv
            xd,yd = data[0][0],data[0][1]
            if datarange is None:
                datarange = (np.min(xd),np.max(xd),np.min(yd),np.max(yd))
            maxmind = (np.max(data[1]),np.min(data[1]))
        elif datarange is None:
            raise ValueError("Can't choose limits for plotting without data or a range hint")
            maxmind = None
        
        grid = np.mgrid[datarange[0]:datarange[1]:1j*nx,datarange[2]:datarange[3]:1j*ny]
        res = self(grid)
        
#        if maxmind:
#            norm = plt.normalize(min(np.min(res),maxmind[1]),max(np.max(res),maxmind[0]))
#        else:
#            norm = plt.normalize(np.min(res),np.max(res))
        
        if clf:
            M.clf()
            
        M.mesh(grid[0],grid[1],res)
        
        if cb:
            if isMappingType(cb):
                M.colorbar(**cb)
            else:
                M.colorbar()
        
        if data:
            if isMappingType(data):
                kwscat = dict(data)
            else:
                kwscat = {}
            kwscat.setdefault('c',data[1])
            M.points3d(xd,yd,**kwscat)
                
        #M.xlim(datarange[0],datarange[1])
        #M.ylim(datarange[2],datarange[3])
            
    def plotResiduals(self,x=None,y=None,fmt='or',xaxis='r',clf=True,relative=False):
        """
        Plots residuals between data and the model for a variety of projections.
        If x and y are None, the residuals will be inferred from the last-fit
        data. For details on how residuals are computed see
        :meth:`FunctionModel.residuals`
        
        :param x: input data
        :type x: array-like
        :param y: output data
        :type y: array-like
        :param fmt: fmt argmument for :func:`matplotlib.pyplot.ploy`
        :param xaxis: 
            Sets the plane in which to plot residuals. Can be
            'r','theta','x', or 'y'.
        :type xaxis: string
        :param clf: If True, the plot will be cleared before plotting.
        :type clf: bool
        :param relative:
            If true, the plot will be relative residuals (e.g. residual/model)
            instead of the absolute residuals.
        
        """
        from matplotlib import pyplot as plt
        
        isinter = plt.isinteractive()
        try:
            x,y,res = self.residuals(x,y,retdata=True)
            if relative:
                res = res/self(x)
            
            if xaxis=='r':
                outcoordsys = 'polar'
                cind = 0
            elif xaxis=='theta':
                outcoordsys = 'polar'
                cind = 1
            elif xaxis=='x':
                outcoordsys = 'cartesian'
                cind = 0
            elif xaxis=='y':
                outcoordsys = 'cartesian'
                cind = 1
            else:
                raise ValueError('unrecognized x-axis type')
            if clf:
                plt.clf()
            xt = self.transformCoordinates(x,outcoordsys=outcoordsys)[cind]
            plt.plot(xt,res,fmt)
            plt.axhline(0,color='k',ls=':')
            plt.xlabel(xaxis)
            
            if isinter:
                    plt.draw()
                    plt.show()
        finally:
            plt.interactive(isinter)
    
class FunctionModel2DScalarAuto(FunctionModel2DScalar):
    """
    A :class:`FunctionModel2DScalar` that has its parameters automatically 
    determined by the :meth:`f` method.
    
    Equivalent to setting the metaclass to :class:`AutoParamsMeta`.
    """
    __metaclass__ = AutoParamsMeta
    
class CompositeModel2DScalar(FunctionModel2DScalar,CompositeModel):
    """
    This model is a composite model of :class:`FunctionModel2DScalar` models.
    """
    def __init__(self,*args,**kwargs):
        """
        See :class:`CompositeModel` for initializer arguments
        """
        super(CompositeModel2DScalar,self).__init__(*args,**kwargs)
        for m in self._models:
            if not isinstance(m,FunctionModel2DScalar):
                raise ModelTypeError('Input model %s is not a 2D->scalar model'%m)
        
        self.incoordsys = 'cartesian'
        
    def _getIncoordsys(self):
        return self._incoordsys
    def _setIncoordsys(self,val):
        for m in self._models:
            m.incoordsys = val
        #this ensures the conversions never happen at the composite level
        self._fcoordsys = self._incoordsys = val 
    incoordsys = property(_getIncoordsys,_setIncoordsys,doc='Input coordinate system')
    
    
class FunctionModel2DScalarSeperable(FunctionModel2DScalar):
    """
    A :class:`FunctionModel2DScalar` that is seperable and follows a radial and 
    polar function that are seperable - e.g. 
    :math:`F(r,\\theta) = R(r)*T(\\theta)`
    """
    def __init__(self,rmodel,thetamodel=None):
        """
        :param rmodel: 
            A FunctionModel1D that is to be used as the radial function.
        :type rmodel: A valid input to :func:`get_model_instance`
        :param rmodel: 
            A FunctionModel1D that is to be used as the polar function.
        :type rmodel: A valid input to :func:`get_model_instance`
        """
        if rmodel is None:
            self.__dict__['_rmodel'] = None
        elif isinstance(rmodel,FunctionModel1D):
            self.__dict__['_rmodel'] = rmodel
        else:
            self.__dict__['_rmodel'] = get_model_instance(rmodel,FunctionModel1D)
            
        if thetamodel is None:
            self.__dict__['_thetamodel'] = None    
        elif isinstance(thetamodel,FunctionModel1D):
            self.__dict__['_thetamodel'] = thetamodel
        else:
            self.__dict__['_thetamodel'] = get_model_instance(thetamodel,FunctionModel1D)
            
        self._fixParams()
    
    _fcoordsys = 'polar'
    def f(self,x,*params):
        self.parvals = params
        R = self._rmodel(x[0]) if self._rmodel is not None else 1
        T = self._thetamodel(x[1]) if self._thetamodel is not None else 1
        return R*T
    
    def _gettmod(self):
        return self._thetamodel
    def _settmod(self,val):
        if not (val is None or isinstance(val,FunctionModel1D)):
            raise TypeError('input value is not a FunctionModel1D')
        self._thetamodel = val
        self._fixParams()
        
    def _getrmod(self):
        return self._rmodel
    def _setrmod(self,val):
        if not (val is None or isinstance(val,FunctionModel1D)):
            raise TypeError('input value is not a FunctionModel1D')
        self._rmodel = val
        self._fixParams()
        
    def _fixParams(self):
        ps = list()
        if self._rmodel is not None:
            ps.extend(self._rmodel.params)
        if self._thetamodel is not None:
            ps.extend(self._thetamodel.params)
        self._pars = ps
        
    def __getattr__(self,name):
        if name in self.params:
            if hasattr(self._rmodel,name):
                return getattr(self._rmodel,name)
            elif hasattr(self._thetamodel,name):
                return getattr(self._thetamodel,name)
        else:
            raise AttributeError(name)
        
    def __setattr__(self,name,val):
        if name in self.params:
            if hasattr(self._rmodel,name):
                return setattr(self._rmodel,name,val)
            elif hasattr(self._thetamodel,name):
                return setattr(self._thetamodel,name,val)
        else:
            super(FunctionModel2DScalarSeperable,self).__setattr__(name,val)
            
    @property
    def rangehint(self):
        if self._rmodel is None:
            rh = (1,1)
        else:
            rh = self._rmodel.rangehint
            
        if rh is None:
            return None
        else:
            return -rh[1],rh[1],-rh[1],rh[1] #assume theta is not localized
            
class FunctionModel2DScalarDeformedRadial(FunctionModel2DScalar):
    """
    A `FunctionModel2DScalar` with a profile that is flattened along an axis 
    but has a radial form if unflattened
    
    atob is a ratio of the large to small axis, and pa is the
    position angle measured in radians from the x-axis to the y-axis 
    (e.g. west to north)
    """
    def __init__(self,rmodel,atob=1,pa=0):
        if isinstance(rmodel,FunctionModel1D):
            self.__dict__['_rmodel'] = rmodel
        else:
            self.__dict__['_rmodel'] = get_model_instance(rmodel,FunctionModel1D)
        self.__dict__['_rmodel'].incoordsys = 'cartesian'
        self._fixPars()
        self.atob = atob
        self.pa = pa
        
    _fcoordsys = 'cartesian'    
    def f(self,arrin,*params):
        #self.parvals = params
        params = list(params)
        pa = params.pop()
        atob = params.pop()
        
        #xr=(-sinpa*x+cospa*y)/cosi
        #yr=(x*cospa+y*sinpa)
        x,y = arrin
        sinpa,cospa = np.sin(pa),np.cos(pa)
        yr=(-sinpa*x+cospa*y)*atob
        xr=(x*cospa+y*sinpa)
        return self._rmodel((xr*xr+yr*yr)**0.5)
        
        
    def _fixPars(self):
        ps = list(self._rmodel.params)
        ps.append('atob')
        ps.append('pa')
        self.__dict__['_pars'] = tuple(ps)
        
    def _getrmodel(self):
        return self._rmodel
    def _setrmodel(self,val):
        if not isinstance(val,FunctionModel1D):
            raise TypeError('input value is not a FunctionModel1D')
        self._rmodel = val
        self._rmodel.incoordsys = 'cartesian'
        self._fixPars()
        
    def __getattr__(self,name):
        if hasattr(self._rmodel,name):
            return getattr(self._rmodel,name)
        else:
            raise AttributeError(name)
        
    def __setattr__(self,name,val):
        if name in self._rmodel.params:
            return setattr(self._rmodel,name,val)
        else:
            super(FunctionModel2DScalarDeformedRadial,self).__setattr__(name,val)
            
    def _getPaDeg(self):
        return np.degrees(self.pa)
    def _setPaDeg(self,val):
        self.pa = np.radians(val)     
    padeg=property(_getPaDeg,_setPaDeg,"""position angle (from +x-axis
            to +y-axis) in degrees
            """)
    
    def _getInc(self):
        return np.arccos(1/self.atob)
    def _setInc(self,val):
        self.atob = 1/np.cos(val)
    inc=property(_getInc,_setInc,doc="""
        inclination angle of object in radians - maps onto atob
        """)
        
    def _getIncDeg(self):
        return np.degrees(np.arccos(1/self.atob))
    def _setIncDeg(self,val):
        self.atob = 1/np.cos(np.radians(val))
    incdeg=property(_getIncDeg,_setIncDeg,doc="""
        inclination angle of object in degrees - maps onto atob
        """)
        
    def _getE(self):
        return (1-(1/self.atob)**2)**0.5
    def _setE(self,val):
        self.atob = (1- val**2)**-0.5
    e = property(_getE,_setE,doc='ellipticity given by e^2 = 1-(b/a)^2')
            
    @property
    def rangehint(self):
        rh = self._rmodel.rangehint
        if rh is None:
            return None
        else:
            return -rh[1],rh[1],-rh[1],rh[1] #assume theta is not localized
        
    
        
#<-------------------------------Module functions ----------------------------->  
__model_registry={}
def register_model(model,name=None,overwrite=False,stripmodel=True):
    """
    Registers a model at the module package level for :func:`get_model` and
    :func`list_model`.
    
    :param model: The model to register
    :type model: a class object that is a subclass of class:`ParametricModel`
    :param name: 
        The name to assign to this model (lower case class name will be used if
        this is None)
    :type name: string or None
    :param overwrite:
        If True, if a model already exists with the provided name, it will be
        silently overwritten. Otherwise a KeyError will be raised.
    :type overwrite: bool
    :param stripmodel:
        If True, the characters 'model' will be stripped from the name before
        the model is assigned.
    :type stripmodel: bool
    
    :except KeyError: 
        If `overwrite` is False and a model already exists with the name
        provided.
        
    **Example**
    
    .. testsetup::
    
        from astropysics.models.core import register_model,list_models
    
    .. doctest::
    
        >>> register_model(MyFavoriteModel,name='mymodel',stripmode=True)
        >>> list_models(include=[MyFavoriteModel])
        ['mymodel']
        >>> register_model(MyOtherFavoriteModel,name=None,stripmode=True)
        >>> list_models(include=[MyFavoriteModel])
        ['myotherfavorite']
        >>> register_model(MyMostFavoriteModel,name=None,stripmode=False)
        >>> list_models(include=[MyFavoriteModel])
        ['mymostfavoritemodel']
    """
    if not issubclass(model,ParametricModel):
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
    
def get_model(model,baseclass=None):
    """
    Returns the class object for the requested model in the model registry
    
    :param model: can be a name, instance, or class object
    :param baseclass:
        If not None, specifies that the requested model me a subclass of the
        provided base class. If not, a TypeError will be raised.
    
    :returns: A class object for the requested model.
    
    :except KeyError: If a model name is not present in the registry
    :except TypeError: For all other invalid inputs
    """
    from inspect import isclass
    
    if isinstance(model,basestring):    
        res = __model_registry[model]
    elif isclass(model):
        if issubclass(model,ParametricModel):
            res = model
        else:
            raise TypeError('class object is not a Model')
    elif issubclass(model.__class__,ParametricModel):
        res = model.__class__
    else:
        raise TypeError('attempted to get invalid model')
    
    if baseclass:
        if not issubclass(res,baseclass):
            raise TypeError('%s is not a subclass of %s'%(res,baseclass))
        
    return res

def get_model_instance(model,baseclass=None,**kwargs):
    """
    Returns an instance of the supplied model - if the input is actually an 
    instance of a model, the same instance will be passed in - otherwise, 
    a new instance will be created.
    
    :param model: A model name, model class, or model object.
    :param baseclass: 
        If not None, ensures that the instance is a subclass of the provided
        base class.
    :type baseclass: class or None
    
    kwargs will either be passed into the new model constructor or applied as
    attributes to an already-existing model object.
    
    :returns: An instance of the request model type
    
    :except ValueError: if the model is not a subclass of the `baseclass`
    """
    if isinstance(model,ParametricModel if baseclass is None else baseclass):
        for k,v in kwargs.iteritems():
            setattr(model,k,v)
        return model
    else:
        return get_model(model,baseclass)(**kwargs)
        

def list_models(include=None,exclude=None,baseclass=None):
    """
    Lists the registered model objects in the package, possibly subject to
    constraints.
    
    :param include:
        A sequence of model names to include in the list (e.g. the function will
        just validate that they are valid models and return the strings)
    :type include: sequence of strings or None
    :param exclude:
        A sequence of model names that should be excluded.
    :type exclude: sequence of strings or None
    :param baseclass:  
        If not None, all models that are not subclasses of this class will be
        filtered out of the results.
    :type baseclass: a class or None
    
    :returns: 
        A list of strings for the models that can be used with :func:`get_model`
        or :func:`get_model_instance`.
    
    :except ValueError: if any provided model strings are not in the registry
    :except TypeError: if both include and exclude are not None
    
    """
    from operator import isSequenceType
    
    if baseclass is not None:
        res = [k for k,c in __model_registry.iteritems() if issubclass(c,baseclass)]
    else:
        res = __model_registry.keys()
        
    if include is not None:
        if exclude is not None:
            raise TypeError("can't specify both included models and excluded models")
        if isinstance(include,basestring):
            include = include.split(',')
        elif not isSequenceType(include):
            include = [include]
        for m in include:
            if m not in res:
                raise ValueError('modelname to include %s not a valid model name'%m)
            res = include
    elif exclude is not None:
        if isinstance(exclude,basestring):
            exclude = exclude.split(',')
        elif not isSequenceType(exclude):
            exclude = [exclude]
        for m in exclude:
            if m not in res:
                raise ValueError('modelname to exclude %s not a valid model name'%m)
            res.remove(m)
    
    return res


def intersect_models(m1,m2,bounds=None,nsample=1024,full_output=False,**kwargs):
    """
    Determine the points where two models intersect.
    
    :param m1: the first model
    :type m1: a :class:`FunctionModel1D` object
    :param m2: the second model
    :type m2: a :class:`FunctionModel1D` object
    :param bounds:
        If None, the bounds will be determined from the model data if any is
        saved. Otherwise, a 2-tuple (min,max) defining the region to search for
        intersections.
    :type bound: tuple or None
    :param nsample: number of points within the bounds to sample for intersections.
        
    :returns:
        A sorted array of points where the two models intersect on the interval
        (up to a resolution of (max-min)/nsample), or if full_output is True, returns
        (array,scipy.optimize.zeros.RootResults)
    
    kwargs are passed into :func:`scipy.optimize.brentq`
    """
    from scipy.optimize import brentq
    
    if bounds is None:
        data1 = m1.data[:2] if hasattr(m1,'data') and m1.data is not None else None
        data2 = m2.data[:2] if hasattr(m2,'data') and m2.data is not None else None
        
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

def scale_model(model,scaleparname='A'):
    """
    A convinience function to generate a CompositeModel with a scaling 
    factor multiplying the wrapped model.
    
    A convinience to generate a CompositeModel with a scaling factor multiplying
    the wrapped model (e.g. :math:`A m(x)`).
    
    :param model: the model to wrap (or a name of a model that will be created)
    :type model: :class`FunctionModel` object or string
    :param scaleparname:
        name for the parameter controlling the scaling factor
    :type scaleparname: string
    
    :returns: a :class:`CompositeModel` (or subclass) object
    """
    model = get_model_instance(model)
    if scaleparname in model.params:
        scaleparname += '1'
    if isinstance(model,FunctionModel1D):
        compclass = CompositeModel1D
    else:
        compclass = CompositeModel
    return compclass((model,'constant'),operation='*',
                     parnames={'C1':scaleparname})

def offset_model(model,offsetparname='C'):
    """
    A convinience to generate a CompositeModel with an additive offset 
    on the wrapped model (e.g. :math:`m(x)+C`). 
    
    :param model: the model to wrap (or a name of a model that will be created)
    :type model: :class`FunctionModel` object or string
    :param offsetparname:
        name for the parameter controlling the additive offset value
    :type offsetparname: string
    
    :returns: a :class:`CompositeModel` (or subclass) object
    """
    model = get_model_instance(model)
    if offsetparname in model.params:
        offsetparname += '1'
    if isinstance(model,FunctionModel1D):
        compclass = CompositeModel1D
    else:
        compclass = CompositeModel
    return compclass((model,'constant'),operation='+',
                     parnames={'C1':offsetparname})
                     
def scale_and_offset_model(model,scaleparname='A',offsetparname='C'):
    """
    A convinience function to generate a CompositeModel with a multiplicative
    scale and an additive offset on the wrapped model. (e.g. :math:`A m(x)+C`)
    
    :param model: the model to wrap (or a name of a model that will be created)
    :type model: :class`FunctionModel` object or string
    :param scaleparname:
        name for the parameter controlling the scaling factor
    :type scaleparname: string
    :param offsetparname:
        name for the parameter controlling the additive offset value
    :type offsetparname: string
    
    :returns: a :class:`CompositeModel` (or subclass) object
    """
    model = get_model_instance(model)
    if scaleparname in model.params:
        scaleparname += '1'
    if offsetparname in model.params:
        offsetparname += '2'
    if isinstance(model,FunctionModel1D):
        compclass = CompositeModel1D
    else:
        compclass = CompositeModel
    return compclass((model,'constant','constant'),operation=['*','+'],
                     parnames={'C1':scaleparname,'C2':offsetparname})
    
def binned_weights(values,n,log=False):
    """
    Produces an array of weights that are generated by subdividing the values
    into n bins such that each bin has an equal share of the total number of
    values.
    
    :param values: the input values
    :type values: array-like
    :param n: number of bins
    :type n: int
    :param log: 
        If True, the values are evenly-spaced on logarithmic intervals,
        otherwise, linear.
    :type log: bool
    
    :returns: An array of weights on [0,1] with shape matching `values`
    
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

def intrinsic_to_observed_ellipticity(ei,i,degrees=True):
    """
    Converts intrinsic ellipticity to observed where :math:`e^2 = 1-(b/a)^2`
    
    :param ei: intrinsic ellipticity
    :type ei: float or array-like
    :param i: inclination angle
    :type i: float or array-like
    :param degrees: if True, the inclination is assumed to be in degrees
    
    :returns: observed ellipticity
    """
    if degrees:
        i = np.radians(i)
        
    return 1 - ((1-ei)**2*np.sin(i)**2+np.cos(i)**2)**0.5
    
def observed_to_intrinsic_ellipticity(eo,i,degrees=True):
    """
    Converts observed ellipticity to intrinsic where :math:`e^2 = 1-(b/a)^2`
    
    :param eo: observed ellipticity
    :type eo: float or array-like
    :param i: inclination angle
    :type i: float or array-like
    :param degrees: if True, the inclination is assumed to be in degrees
    
    :returns: intrinsic ellipticity
    """
    if degrees:
        i = np.radians(i)
        
    return 1 - ((1-eo)-np.cos(i)**2)**0.5/np.sin(i)
