from __future__ import division
from math import pi
import numpy as np

try:
    import spylot
except ImportError:
    print 'Spylot not found - some spec package functions may not work correctly'
    
class Spectrum(object):
    """
    Represents a flux/luminosity/intensity
    
    Note that operations are performed in-place, and properties retrieve the 
    same versions that are changed (except ivar)
    """
    def __init__(self,x,flux,err=None,ivar=None,unit='wl',copy=True):
        """
        sets the x-axis values and the flux.  Optionally, an error can be
        given by supplying an error or an inverse variance (bot not both)
        
        copy determines if the inputs will be copied if they are already arrays
        """
        
        x = np.array(x,copy=copy)
        flux = np.array(flux,copy=copy)
        if x.shape != flux.shape:
            raise ValueError("x and flux don't match shapes")
        
        if ivar is None and err is None:
            err = np.zeros_like(flux)
        elif ivar is not None:
            err = np.array(ivar,copy=False)**-0.5
            if err.shape != flux.shape:
                raise ValueError("ivar and flux don't match shapes")
            
        elif err is not None:
            err=np.array(err,copy=copy)
            err[err<0]=-1*err[err<0]
            if err.shape != flux.shape:
                raise ValueError("err and flux don't match shapes")
            
        else:
            raise ValueError("can't set both err and ivar at the same time")
        
        self._phystype,self._unit,self._scaling = self._strToUnit(unit)
        
        self._x = x
        self._flux = flux
        self._err = err
        
    def _strToUnit(self,typestr):
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
        elif u == 'm' or u == 'wavelength-m':
            val =  'wavelength-m'
            scaling=1e2
        elif u == 'cm' or u == 'wavelength-cm':
            val =  'wavelength-cm'
            scaling=1
        elif u == 'f' or u == 'nu' or u == 'hz' or u == 'frequency-hz':
            val =  'frequency-Hz'
            scaling=1
        elif u == 'thz' or u == 'frequency-THz':
            val =  'frequency-THz'
            scaling=1e12
        elif u == 'e' or u == 'en' or u == 'energy' or u == 'energy-eV':
            from .constants import ergperev
            val =  'energy-eV'
            scaling=ergperev
        elif u == 'erg' or u == 'energy-erg':
            val =  'energy-erg'  
            scaling=1
        elif u == 'J' or u == 'energy-J':
            val =  'energy-J'   
            scaling=1e-7
        else:
            raise ValueError('unrecognized unit')
        
        phystype,unit = val.split('-')
        return phystype,unit,scaling
        
    def _convertType(self,oldtype,newtype,oldscale,newscale):
        #TODO:fix wl to freq
        from .constants import c,h
        
        x,flux,err = self._x/oldscale,self._flux*oldscale,self._err*oldscale # convert to cgs
    
        if newtype == 'energy': #TODO:check
            if oldtype == 'energy': 
                fluxscale = 1
            elif oldtype == 'wavelength':
                x = h*c/x
                fluxscale = x*x/c/h #convert to fnu and divide by h
                
            elif oldtype == 'frequency':
                x = h*x
                fluxscale = 1/h
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'wavelength':
            if oldtype == 'energy': #TODO:check
                x = h*c/x
                fluxscale = x*x*h/c #convert to fnu and then to  by h
            elif oldtype == 'wavelength':
                fluxscale = 1
            elif oldtype == 'frequency':
                x = c/x
                fluxscale = x*x/c #flambda = fnu*nu^2/c
            else:
                raise ValueError('unrecognized oldtype')
        elif newtype == 'frequency':
            if oldtype == 'energy': #TODO:check
                x = x/h
                fluxscale = h
            elif oldtype == 'wavelength':
                x = c/x
                fluxscale = x*x/c #fnu = flambda*lambda^2/c
            elif oldtype == 'frequency':
                fluxscale = 1
            else:
                raise ValueError('unrecognized oldtype')
        else:
            raise ValueError('unrecognized newtype')
        
        return x*newscale,flux*fluxscale/newscale,err*fluxscale/newscale
        
    
    
    
    #------------------------Properties--------------------------------->
    @property
    def npix(self):
        return self._x.shape[-1]
    
    @property
    def shape(self):
        return self._x.shape
    
    def _getFlux(self):
        return self._flux
    def _setFlux(self,flux):
        flux = np.array(flux,copy=False)
        if flux.shape != self._flux.shape:
            raise ValueError("new flux doesn't match old flux shape")
        self._flux[:] = np.array(flux)
    flux = property(_getFlux,_setFlux)
    
    def _getNFlux(self):
        from .constants import h
        
        oldtype = self.unit
        try:
            self.unit = 'hz'
            return self._flux/h/self._x
        finally:
            self.type = oldtype
    def _setNFlux(self,nflux):
        nflux = np.array(nflux)
        if nflux.shape != self._flux.shape:
            raise ValueError("new flux doesn't match flux shape")
        raise NotImplementedError
        self._flux = 42
    nflux = property(_getNFlux,_setNFlux)
    
    def _getErr(self):
        return self._err
    def _setErr(self,err):
        if err.shape != self._err.shape:
            raise ValueError("new err doesn't match old err shape")
        self._err[:]=np.array(err)
    err = property(_getErr,_setErr)
    
    def _getIvar(self):
        return 1/self._err/self._err
    def _setIvar(self,ivar):
        if ivar.shape != self._flux.shape:
            raise ValueError("newivar doesn't match flux shape")
        self._err=np.array(ivar)**-0.5
    ivar = property(_getIvar,_setIvar)
    
    def _getX(self):
        return self._x
    def _setX(self,x):
        x = np.array(x)
        if x.shape != self._flux.shape:
            raise ValueError("new x doesn't match flux shape")
        self._x = np.array(x)
    x = property(_getX,_setX)
    
    def _getUnit(self):
        return self._phystype+'-'+self._unit
    def _setUnit(self,typestr):
        newtype,newunit,newscaling = self._strToUnit(typestr)
        
        x,flux,err = self._convertType(self._phystype,newtype,self._scaling,newscaling)
        
        self._x[:] = x
        self._flux[:] = flux
        self._err[:] = err
        self._phystype = newtype
        self._unit = newunit
        self._scaling = newscaling
    unit = property(_getUnit,_setUnit)
    
    #<----------------------Tests/Info--------------------------------->
    def __len__(self):
        return len(self._x)
    
    def getDx(self,mean=True):
        """
        get the spacing of the x-axis, which is always 1 element shorter than x
        
        if mean, returns the mean of the spacing
        """
        dx = np.convolve(self._x,(1,-1),mode='valid')
        if mean:
            return dx.mean()
        else:
            return dx
        
    def getDlogx(self,mean=True,logbase=10):
        """
        get the logarithmic spacing of the x-axis, which is always 1 element
        shorter than x
        
        if mean, returns the mean of the spacing
        """
        x = np.log(self._x)/np.log(logbase)
        dlogx = np.convolve(x,(1,-1),mode='valid')
        if mean:
            return dlogx.mean()
        else:
            return dlogx
    
    def isLinear(self,eps=1e-10):
        """
        Determines if the x-spacing is linear (e.g. evenly spaced so that 
        dx ~ constant with a standard deviation of eps)
        """
        return np.std(self.getDx(False)) < eps
        
    def isLogarithmic(self,eps=1e-10):
        """
        Determines if the x-spacing is logarithmic (e.g. evenly spaced in 
        logarithmic bins so that dx ~ x to tolerance of eps)
        """
        return np.std(self.getDlogx(False)) < eps
    
    #<----------------------Operations---------------------------->
    
    def smooth(self,width=1,filtertype='gaussian',replace=True):
        """
        smooths the flux in this object by a filter of the given filtertype 
        (can be either 'gaussian' or 'boxcar')
        
        if replace is True, the flux in this object is replaced by the smoothed
        flux and the error is smoothed in the same fashion
        
        returns smoothedflux
        """
        import scipy.ndimage as ndi 
        
        if filtertype == 'gaussian':
            filter = ndi.gaussian_filter1d
        elif filtertype == 'boxcar' or type == 'uniform':
            filter = ndi.uniform_filter1d
        else:
            raise ValueError('unrecognized filter type %s'%filtertype)
        
        smoothed = filter(self._flux,width)
        
        if replace:
            self.flux = smoothed
            self.err = filter(self._err,width)
        
        return smoothed
    
    def resample(self,newx,interpolation='linear',replace=True,**kwargs):
        """
        this interpolates the flux to populate the supplies x-axis
        
        kwargs go into the interpolation routine as described below
        
        interpolations can be:
        'linear': simple linear interpolation
        'spline': a k-order spline with smoothing factor s is used, where s and 
        k are set by kwargs.  if the 'save' kwarg is True, the spline is saved
        and will be used for subsequent resamplings.  if 'clear' is True, the 
        existing spline will be cleared and a new spline will be calculated.
        
        WARNING: this does not treat the errors properly yet - currently just interpolates
        
        returns newx,newflux,newerr
        """
        if interpolation == 'linear':
            newflux = np.interp(newx,self._x,self._flux)
            #TODO: fix errors
            newerr = np.interp(newx,self._x,self._err)
        elif 'spline' in interpolation:
            from scipy.interpolate import UnivariateSpline
            
            s = kwargs.pop('s',None)
            k = kwargs.pop('k',3)
            save = kwargs.pop('save',False)
            clear = kwargs.pop('clear',False)
            
            fspline = espline = None
            if hasattr(self,'_spline') and not (clear or save):
                fspline,espline = self._spline
                
                    
            if fspline is None:
                fspline = UnivariateSpline(self._x,self._flux,k=k,s=s)
            if espline is None:
                espline = UnivariateSpline(self._x,self._flux,k=k,s=s)
                
            if save:
                self._spline = (fspline,espline)
            
            newflux = fspline(newx)
            #TODO: fix errors
            newerr = espline(newx)
        else:
            raise ValueError('unrecognized interpolation technique')
        
        if len(kwargs)>0:
            raise TypeError('unexpected keyword %s'%kwargs.keys()[0])
        
        if replace:
            if newx.shape != self._x.shape:
                raise ValueError("can't replace if new x-axis shape doesn't match old")
            
            self._x[:] = newx 
            self._flux[:] = newflux
            self._err[:] = newerr   
        return newx,newflux,newerr
    
    def linearize(self,lower=None,upper=None,**kwargs):
        """
        convinience function for resampling to an equally-spaced linear x-axis
        
        if lower or upper are None, the upper and lower x values are used
        kwargs go into Spectrum.resample
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        
        newx = np.linspace(lower,upper,self.npix)
        return self.resample(newx,**kwargs)
    
    def logify(self,lower=None,upper=None,**kwargs):
        """
        convinience function for resampling to an x-axis that is equally spaced
        in logarithmic bins.  Note that lower and upper are the x-axis values 
        themselves, NOT log(xvalue)
        
        if lower or upper are None, the upper and lower x values are used
        """
        if lower is None:
            lower = np.min(self._x)
        if upper is None:
            upper = np.max(self._x)
        
        newx = np.logspace(np.log10(lower),np.log10(upper),self.npix)
        return self.resample(newx,**kwargs)
    
    def plot(self,fmt=None,ploterrs=True,smoothing=None,clf=True,**kwargs):
        """
        uses matplotlib to plot the Spectrum object
        
        kwargs go into 
        """
        import matplotlib.pyplot as plt
        
        x,y,e = self.x,self.flux,self.err
        
        if fmt is None:
            res = [plt.plot(x,y,**kwargs)]
        else:
            res = [plt.plot(x,y,fmt,**kwargs)]
            
        if ploterrs and np.any(e):
            from operator import isMappingType
            if not isMappingType(ploterrs):
                ploterrs = {}
            ploterrs.setdefault('ls','--')
            
            res.append(plt.plot(x,e,**ploterrs))
        plt.xlim(np.min(x),np.max(x))
            
        return res
        
    
    
def align_spectra(specs,ressample='super',interpolation='linear',copy=False):
    """
    resample the spectra in the sequence specs so that they all have the same 
    x-axis.  
    
    ressample can be 'super' lower-resolution spectra to the highest or 'sub' 
    to interolate the highest-resolution spectrum to the lowest
    alternateively, 'logsuper' or 'logsub' will use the logarithmic resolution
    to determine
    
    copy makes new copies of the Spectrum objects 
    
    returns specs, or new copies if copy is True
    """
    from operator import isSequenceType
    if not isSequenceType(specs):
        raise ValueError('specs must be a sequence')
    
    for s in specs:
        if s.__class__.__name__ != 'Spectrum':
            raise ValueError(str(s)+' is not a spectrum')
    
    if copy:
        from copy import deepcopy
        specs= [deepcopy(s) for s in specs]
    
    if 'super' in ressample:
        super = True
    elif 'sub' in ressample:
        super = False
    else:
        raise ValueError('unrecognized ressample value')
    
    logres = 'log' in ressample
        
    
    
    
    if logres:
        reses=[s.getDlogx() for s in specs]
    else:
        reses=[s.getDx() for s in specs]
    mins=[np.min(s) for s in specs]
    maxs=[np.max(s) for s in specs]
    
    if super:
        templi = np.where(reses == max(reses))[0][0]
    else:
        templi = np.where(reses == max(reses))[0][0]
    
    x = specs[templi].x
    
    
    for s in specs:
        s.resample(x,interpolation)
    
    return specs

def zfind(object,templates,lags=True,dochecks = True, interpolation='linear',verbose = False):
    """
    Computes the lag with the best-fitting chi-squared over all the templates
    
    templates should be nstars X npix 
    
    lags specify the possible pix-lags to try - if True, it will try all 
    possibilities, if a sequence of numbers, will try all supplied lags, and 
    if False/None, does no lag
    
    dochecks makes sure the spectra are all logarithmically spaced and interpolates
    if necessary
    
    interpolation determines the type of interpolation to use whenever it is necessary
    
    returns bestz,besti,chi2s,zs,coeffs(nzs X nstars),dofs
    """
    from operator import isSequenceType,isMappingType
    
    if dochecks:
        from warnings import warn
        
        if hasattr(object,'flux'): #allows for type-matching
            pass
        elif isSequenceType(object):
            object = Spectrum(*object)
        elif isMappingType(object):
            object = Spectrum(**object)
        else:
            raise TypeError("couldn't interpret input object")
        
        if not object.isLogarithmic():
            warn('object is not logarithmic - interpolating to logarithmic spacing')
            object.logify()
        
        
        if isinstance(templates,Spectrum):
            templates=[templates]
        elif isSequenceType(templates):
            temptemplates=[]
            for t in templates:
                if hasattr(t,'flux'):#assume this is a flux-object
                    temptemplates.append(t)
                elif isSequenceType(t):
                    if len(t) == len(object): #assume this is a flux array
                        temptemplates.append(Spectrum(object.x,t))
                    else:
                        temptemplates.append(Spectrum(*t))
                elif isMappingType(t):
                    temptemplates.append(Spectrum(**t))
                else:
                    raise TypeError("couldn't interpret input object")
            templates = temptemplates
            del temptemplates
        else:
            raise TypeError("input templates not a sequence")
        
        #search for first logarithmically spaced template and use that one as the basis
        tgoodi = None
        for i,t in enumerate(templates):
            if t.isLogarithmic():
                tgoodi = i
                break
        
        if tgoodi is None:
            warn('No logarithmic templates found - interpolating to logarithmic spacing for template 0')
            templates[0].logify()
            tgoodi = 0
            
        tgoods = templates[tgoodi].shape
        tgoodx = templates[tgoodi].x
            
        for i,t in enumerate(templates):
            if t.shape != tgoods or np.any(t.x != tgoodx):
                warn('Template %i does not match base template %i - interpolating'%(i,tgoodi))
                t.resample(tgoodx,interpolation=interpolation)
        
                #TODO:resample smarter so that there aren't 0-set areas
        if object.shape != tgoods or np.any(object.x != tgoodx):
            warn('object x-axis does not match templates - interpolating')
            newt = object.resample(tgoodx,interpolation=interpolation,replace=False)
            object = Spectrum(*newt)#TODO:make this more robust?
    
    dlogx = object.getDlogx()
        
    objflux = object.flux
    ivar = object.ivar
    infi=np.isinf(ivar)
    if all(ivar):
        ivar[infi]=1 #equally weight
    else:
        ivar[infi]= 0 #ignore 0 ivars
    tempfluxes=np.array([t.flux for t in templates]).T
    
    npix,nstars = tempfluxes.shape
    if npix != objflux.shape[0] != ivar.shape[0]:
        raise ValueError("pixels don't match")
    
    
    if lags is True:
        lags = np.arange(2*npix-3)-npix+2
    elif not np.any(lags):
        lags=np.array((0,))
    elif isSequenceType(lags):
        lags = np.array(lags,dtype=int)
    else:
        raise ValueError("couldn't understand lags")
    
    dofs=np.ndarray(len(lags),int)
    rchi2s=np.ndarray(len(lags))
    coeffs=np.ndarray((len(lags),nstars))
    fits= np.zeros((len(lags),npix))
    #go over all possible lags, compute chi-squared
    for i,l in enumerate(lags): #TODO: arrayify this?
        if verbose:
            print 'lag=',l
        b = np.matrix(objflux).T
        A = np.matrix(tempfluxes)
        W = np.array(np.tile(ivar,nstars).reshape(A.shape))
        
        if l > 0:
            objrange=slice(l,None)
            matrange=slice(None,-l)
            
            
        elif l < 0:
            objrange=slice(None,l)
            matrange=slice(-l,None)
        else:
            objrange=slice(None,None)
            matrange=slice(None,None)
        
        b = b[objrange]
        A = A[matrange]
        W = W[matrange]
        dof = np.sum(ivar[objrange] > 0) - nstars
        
            
        
        #cos = np.linalg.solve(A,b) #does not work for non-square mats, so:
        #AAi = np.linalg.inv(A.T * A) #often mats are singular
        AAi = np.linalg.pinv(A.T * np.multiply(W,A)) #TODO:speed tests?
        #AAi = np.linalg.pinv(A.T * A) #TODO:speed tests?
        
        coeff = AAi *  ( A.T * np.multiply(W[:,0:1],b) )
        #coeff = AAi *  ( A.T * b )
        
        fit = np.asarray(A * coeff).flatten()
        chis = fit - b.A.flatten()
        chis = chis
        csq = chis*chis*ivar[objrange] #faster than cs**2
        
        coeffs[i] = np.asarray(coeff).flatten()
        fits[i][objrange] = fit
        #fits[i][round((npix-len(fit))/2.0):round((npix+len(fit))/2.0)] = fit
        dofs[i] = dof
        rchi2s[i]=sum(csq)/dof
    minvalid = rchi2s[np.isfinite(rchi2s)].min()
    imin = np.where(rchi2s==minvalid)
    if imin[0].size > 1:
        raise ValueError('found more than 1 chi-squared minimum')
    zs = 10**(lags*dlogx)-1    
        
    #TODO:go bzck to zs
    #return lags[imin],imin,rchi2s,zs,coeffs,dofs,fits
    return zs[imin],imin,rchi2s,zs,coeffs,dofs,fits

    
#<---------------------Data loading functions---TODO:move to plugin-like?------>

def load_deimos_spectrum(fn,plot=True,extraction='horne',retdata=False,smoothing=None):
    """
    extraction type can 'horne' or 'boxcar'
    
    returns Spectrum object with ivar, [bdata,rdata]
    """
    import pyfits
    if 'spec1d' not in fn or 'fits' not in fn:
        raise ValueError('loaded file must be a 1d spectrum from deep pipeline')
    
    if extraction == 'horne':
        extname = 'HORNE'
    elif extraction == 'boxcar':
        extname = 'BXSPF'
    else:
        raise ValueError('unrecgnized extraction type %s'%extraction)
    
    f=pyfits.open(fn)
    try:
        extd = dict([(f[i].header['EXTNAME'],i) for i in range(1,len(f))])
        bi,ri = extd[extname+'-B'],extd[extname+'-R']
        bd,rd=f[bi].data,f[ri].data
        x=np.concatenate((bd.LAMBDA[0],rd.LAMBDA[0]))
        flux=np.concatenate((bd.SPEC[0],rd.SPEC[0]))
        ivar=np.concatenate((bd.IVAR[0],rd.IVAR[0]))
        sky=np.concatenate((bd.SKYSPEC[0],rd.SKYSPEC[0]))
        
        changei = len(bd.LAMBDA[0])
        
        fobj = Spectrum(x,flux,ivar)
        fobj.sky = sky
        
        if smoothing:
            fobj.smooth(smoothing,replace=True)
        
        if plot:
#            from operator import isMappingType
#            if not isMappingType(plot):
#                plot={}
            from matplotlib import pyplot as plt
            plt.plot(fobj.x[:changei],fobj.flux[:changei],'-b')
            plt.plot(fobj.x[changei:],fobj.flux[changei:],'-r')
        
        if retdata:
            return fobj,bd,rd
        else:
            return fobj
    finally:
        f.close()
        
def load_deimos_templates(fns):
    """
    This will generate a dictionary of Spectrum objects from the specified  
    templates 
    fn can be a single file, a list of files, or a unix-style pattern
    """
    import pyfits
    from operator import isSequenceType
    from warnings import warn
    
    if type(fns) == str:
        from glob import glob
        fns = glob(fns)
        
    if not isSequenceType(fns):
        raise ValueError('improper filename format')
    
    tempd={}
    xd={}
    for fn in fns:
        f = pyfits.open(fn)
        try:
            h=f[0].header
            d=f[0].data
            if len(d.shape) == 1:
                d = (d,)
            if len(d) == 1:
                if h.has_key('NAME0'):
                    k=h['NAME0'].strip()
                elif h.has_key('NAME'):
                    k=h['NAME'].strip()
                elif h.has_key('OBJECT'):
                    k=h['OBJECT'].strip()
                else:
                    k = fn.strip()
                    warn('could not infer spectrum name - using filename %s as key name'%k)
                    
                xd[k] = 10**(h['COEFF0']+h['COEFF1']*np.arange(d[0].shape[-1]))
                tempd[k] = d[0].copy()
            else:
                x = 10**(h['COEFF0']+h['COEFF1']*np.arange(d.shape[-1]))
                for i,spec in enumerate(d):
                    if  h.has_key('NAME%i'%i):
                        k=h['NAME%i'%i].strip()
                    elif h.has_key('OBJECT') and h.has_key('EIGEN%i'%i):
                        k = '%s%i'%(h['OBJECT'].strip(),h['EIGEN%i'%i])
                    elif h.has_key('NAME'):
                        k = '%s-%i'%(h['NAME'].strip(),i)
                    elif h.has_key('OBJECT'):
                        k = '%s-%i'%(h['OBJECT'].strip(),i)
                    else:
                        k = '%s-%i'%(fn.strip(),i)
                        warn('could not infer spectrum name - using filename %s as key name'%k)
                        
                    xd[k] = x
                    tempd[k]=spec.copy() #TODO: see if its better to not copy here
                    
        finally:
            f.close()
    
    return dict(((k,Spectrum(xd[k],tempd[k])) for k in tempd))
    
    
         