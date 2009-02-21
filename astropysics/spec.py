#Copyright (c)2008 Erik Tollerud (etolleru@uci.edu) 
"""
This module is for observed or synthetic spectra (Spectrum object) and 
related operations.

Methods/functions tend to be oriented towards optical, but will work in other
bands.
"""

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
            if np.isscalar(ivar):
                err = ivar**-0.5*np.ones_like(flux)
            else:
                err = np.array(ivar,copy=False)**-0.5
            if err.shape != flux.shape:
                raise ValueError("ivar and flux don't match shapes")
            
        elif err is not None:
            if np.isscalar(err):
                err = err*np.ones_like(flux)
            else:
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
        
    def copy(self):
        """
        Generates a deep copy of this Spectrum
        """
        from copy import deepcopy
        return deepcopy(self)
        
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
        
        print oldtype,newtype,oldscale,newscale #DB
        
        x,flux,err = self._x/oldscale,self._flux/oldscale,self._err/oldscale # convert to cgs
    
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
        
        return x*newscale,flux*fluxscale*newscale,err*fluxscale*newscale
        
    
    
    
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
        print newtype,newunit,newscaling #DB
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
    
    def __eq__(self,other):
        try:
            xeq = self._x == other._x
            feq = self._flux == other._flux
            eeq = self._err == other._err
            return np.all(xeq & feq & eeq) 
        except AttributeError:
            return False
    
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
        
    def isXMatched(self,other,tol=1e-10):
        """
        tests if the x-axis of this Spectrum matches that of another Spectrum
        or equal length array, with an average deviation less than tol
        """
        from operator import isSequenceType
        if isSequenceType(other):
            ox = other
        else:
            ox = other.x
            
        try:
            return np.std(self.x - ox) < tol
        except (TypeError,ValueError):
            return False
    
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
        (can be either 'gaussian' or 'boxcar'/'uniform')
        
        if replace is True, the flux in this object is replaced by the smoothed
        flux and the error is smoothed in the same fashion
        
        width is in pixels
        
        returns smoothedflux,smoothederr
        """
        import scipy.ndimage as ndi 
        
        if filtertype == 'gaussian':
            filter = ndi.gaussian_filter1d
        elif filtertype == 'boxcar' or type == 'uniform':
            filter = ndi.uniform_filter1d
        else:
            raise ValueError('unrecognized filter type %s'%filtertype)
        
        smoothedflux = filter(self._flux,width)
        smoothederr = filter(self._err,width)
        
        if replace:
            self.flux = smoothedflux
            self.err = smoothederr
        
        return smoothedflux,smoothederr
    
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
        
        if smoothing:
            x,(y,e) = self.x,self.smooth(smoothing,replace=False)
        else:
            x,y,e = self.x,self.flux,self.err
        
        if 'ls' not in kwargs and 'linestyle' not in kwargs:
            kwargs['ls']='steps'
            
        if clf:
            plt.clf()
        
        if fmt is None:
            res = [plt.plot(x,y,**kwargs)]
        else:
            res = [plt.plot(x,y,fmt,**kwargs)]
            
        if ploterrs and np.any(e):
            from operator import isMappingType
            if not isMappingType(ploterrs):
                ploterrs = {}
            ploterrs.setdefault('ls','-')
            
            res.append(plt.plot(x,e,**ploterrs))
        plt.xlim(np.min(x),np.max(x))
        
        xl=self.unit
        xl=xl.replace('wavelength','\\lambda')
        xl=xl.replace('frequency','\\nu')
        xl=xl.replace('energy','E')
        xl=xl.replace('angstrom','\\AA')
        xl=xl.replace('micron','\\mu')
        xl=tuple(xl.split('-'))
        plt.xlabel('$%s/{\\rm %s}$'%xl)
        
        plt.ylabel('$ {\\rm Flux}/({\\rm erg}\\, {\\rm s}^{-1}\\, {\\rm cm}^{-2} {\\rm %s}^{-1})$'%xl[1])
            
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


def zfind(specobj,templates,lags=(0,200),checkspec=True,checktemplates=True,verbose=True,interpolation = None):
    """
    computes the best fit by linear least-squares fitting of templates to the 
    spectrum for each possible pixel offset.  Weighted fits will be done if 
    the spectrum ivars are different.
    
    lags can either be a sequence of lags or a 2-tuple specifying the lower and 
    upper possible lags
    
    specobj must be a Spectrum object or a sequence of (flux,[x],[ivar]) or 
    flux 
    
    templates can be either a sequence of Spectrum objects or an array with at
    least one dimension matching the pixel dimension.  Long templates are not
    yet supported.
    
    interpolation is the technique for interpolating for sub-pixel lags.  If 
    None, no interpolation is used, so lags must be integers (but this method is
    much faster)
    
    returns besti,lags,coeffs,xs,fitfluxes,rchi2s
    """
    if interpolation is not None:
        raise NotImplementedError('interpolation not implemented yet')
    
    from operator import isSequenceType
    if not isinstance(specobj,Spectrum) and isSequenceType(specobj):
        if len(specobj) > 3:
            specobj = Spectrum(np.logspace(0,1,len(specobj)),specobj)
        else:
            flux = specobj[0]
            if len(specobj) < 3:
                ivar = None
            else:
                ivar = specobj[2]
            if len(specobj) < 2:
                x = np.arange(len(flux))
            else:
                x = specobj[1]
            specobj = Spectrum(x,flux,ivar=ivar)
        
    x = specobj.x
    if checkspec:
        if not specobj.isLogarithmic():
            newx = np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),len(x))
            specobj = specobj.copy()
            specobj.resample(newx)
            x = specobj.x
    flux = specobj.flux
    npix = specobj.shape[0]
    ivar = specobj.ivar
    
    if checktemplates:
        templates=[t if isinstance(t,Spectrum) else Spectrum(x,t) for t in templates]
        for i,t in enumerate(templates):
            if not t.isXMatched(x):
                if verbose:
                    print 'template',i,'does not match spectrum -- resampling'
                t.resample(x)
        templates = [t.flux for t in templates]
    templates=np.atleast_1d(templates)
    
    if templates.shape[0] == npix:
        tm = np.matrix(templates)
    elif templates.shape[1] == npix:
        tm = np.matrix(templates).T
        
    y = np.matrix(flux).T
    ws = np.matrix(ivar).T
    
    
    if type(lags) is tuple and len(lags) == 2:
        lags = np.arange(*lags)
    llags = lags[lags<0]
    ulags = lags[lags>0]
    
    ls,slices=[],[]
    for l in llags:
        ls.append(l)
        #tml = tm[-l:,:]
        #yl = y[:l]
        slices.append(((slice(-l,None),slice(None)),slice(None,l)))
    if 0 in lags:
        ls.append(0)
        slices.append(((slice(None),slice(None)),slice(None)))
    for l in ulags:
        ls.append(l)
        #tml = tm[:-l,:]
        #yl=y[l:]
        slices.append(((slice(None,-l),slice(None)),slice(l,None)))
        
        
    cs,dsq,fitfluxes=[],[],[]
    
    #don't do weighting if all of the errors are identical -- matrix becomes singular
    useweights = np.any(ivar-ivar[0]) and not np.all(~np.isfinite(ivar)) 
    for l,s in zip(ls,slices):
        if verbose:
            print 'doing lag',l
        A = tm[s[0]]
        v = y[s[1]]
        w = ivar[s[1]]
        
        if useweights:
            try:
                AT = np.multiply(A.T,w)
                cs.append(np.linalg.inv(AT*A)*AT*v)
            except np.linalg.LinAlgError,e:
                if verbose:
                    print 'Error inverting matrix in lag',l,':',e
                cs.append(np.linalg.pinv(A)*v)
        else:
            cs.append(np.linalg.pinv(A)*v)
            #TODO: faster inversion schemes?
        
        fitspec = A*cs[-1]
        fitfluxes.append(fitspec)
        fitdiff = (v-fitspec).A
        dsq.append(np.sum(fitdiff*fitdiff))
    ls=np.array(ls)
    cs=np.array(cs)
    dsq=np.array(dsq)
    dofs = np.sum(ivar!=0) - np.abs(ls)
    rchi2s = dsq/dofs
    
    mins = np.where(rchi2s==min(rchi2s))[0]
    if len(mins) == 1:
        besti = mins[0]
    else:
        besti = mins
        
    xs = [x[s[1]] for s in slices]
    fitspecs = [s.A[:,0] for s in fitspecs]
    
    return besti,ls,cs,xs,fitfluxes,rchi2s
        
            
    
    
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
    
def load_spylot_spectrum(s,bandi):
    x=s.getCurrentXAxis()
    f=s.getCurrentData(bandi=bandi)
    if s.isContinuumSubtracted():
        e=s.getContinuumError()
    else:
        e=s.getWindowedRMSError(bandi=bandi)
    return Spectrum(x,f,e)
         