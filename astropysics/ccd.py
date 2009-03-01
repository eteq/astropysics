#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains objects and functions for viewing and/or reducing CCD 
images.

This package deals mostly with the raw CCD images themselves - for science tools
see the phot and spec packages.
"""

from __future__ import division
import numpy as np

class CCDImage(object):
    """
    Represents a CCD image.  
    
    Currently intended as an abstract base class -- subclasses should follow the
    guidelines below to support various formats/file types.
    
    MUST implement:
    *_extractArray(self,range) -- should return the 2d (x,y) array 
    corresponding to the range (xlower,xupper,ylower,yupper) or None for entire 
    image.  Raise IndexError if the range is not valid
    
    *_applyRange(self,range,data) -- should take the provided range and apply 
    the data to whatever the base store is for that range.  Raise IndexError if
    the range is not valid, or AttributeError if the data should be read-only.
    
    MAY implement:    
    * close(self) do any necessary file cleanup (default does nothing)
    """
    def __init__(self,range=None,scaling=None):
        #informational property attributes
        self._pixscale = None
        self._zpt = None
        
        #internal variables
        self.__examcid = None
        self._changed = False
        self._fstatd = {}
        self._lstatd = {}
        self._scaling = ''
        self.applyChangesOnActivate = False#True #TODO:check this out for consistency
        
        self._rng = range
        self._scalefunc = self._invscalefunc = None
        self.setScaling(scaling) #implicitly calls activateRange and will NotImplementedError if not overridden
        
    def __del__(self):
        self.close()
    
    def activateRange(self,range):
        """
        chooses the range to activate for further operations
        
        can either be None (whole file), a 4-tuple of
        the form (xl,xu,yl,yu) or a 3-tuple of the form (xcen,ycen,radius)
        
        Note that the 3-tuple form will be range-checked and auto-shrunk if the 
        radius is too big, but the 4-tuple will not
        """
            
        if self.applyChangesOnActivate and self._changed:
            try:
                self.applyChanges()
            except:
                from warnings import warn
                warn("No inverse function available - can't save",category = RuntimeWarning)
                
        if range is None:
            pass
        elif len(range) == 3:
            xcen,ycen,rad = range
            nx,ny = self.shape
            if xcen < 0 or xcen >= nx or ycen < 0 or ycen >= ny:
                raise IndexError('center of range %i,%i is not in image'%(xcen,ycen))
            range=[xcen-rad,xcen+rad,ycen-rad,ycen+rad]
            
            xl = xcen-rad if xcen-rad > 0 else 0
            xu = xcen+rad if xcen+rad <= nx else nx
            yl = ycen-rad if ycen-rad > 0 else 0
            yu = ycen+rad if ycen+rad <= ny else ny

            range=(xl,xu,yl,yu)
            
        elif len(range) == 4: 
            range=tuple(range) #no range checking
        else:
            raise ValueError('Unregonized form for range')
            
        im = self._extractArray(range)
        
        self._rng = range
        self._active = self._scalefunc(im)
        
    def applyChanges(self):
        if not self._invscalefunc:
            raise Exception('No inverse function available')
        if self._invscalefunc is self._scalefunc: #linear mapping so changes already applied
            pass
        else:
            altim = self._invscalefunc(self._active)
            range = self._rng
            
            self._applyArray(range,altim)
            
        self._changed = False
        self._lstatd.clear() 
        self._gstatd.clear()
        
    def _extractArray(self,range):
        raise NotImplementedError
        
    def _applyArray(self,range,data):
        raise NotImplementedError
    
    def close(self):
        """
        default close method does nothing, override if file cleanup is necessary
        """
        pass
        
    def setScaling(self,scalefunc,invscalefunc=None):
        """
        scalefunc must be a callable or evaluate to False (linear mapping)
        
        alternatively, it may be a string specifying one of the builtin types
        of scaling
        
        if inverse function is None, changes cannot be saved
        """
        if scalefunc:
            if isinstance(scalefunc,basestring):
                if self._scaling == scalefunc: 
                    return #if same scaling as before, do nothing
                else:
                    sfl=scalefunc.lower()
                    if sfl=='linear':
                        self.setScalingLinear()
                    elif sfl=='log':
                        self.setScalingLog()
                    elif sfl=='sb':
                        self.setScalingSurfaceBrightness()
                    elif sfl=='exp':
                        self.setScalingExp()
                    elif sfl=='power':
                        self.setScalingPower()
                    elif sfl=='asinh' or sfl=='arcsinh':
                        self.setScalingASinh()
                    else:
                        raise ValueError('Unrecognized scaling string')
                    
            elif not callable(scalefunc):
                raise ValueError('scalefunc is not a callable')
            newf=scalefunc
            newif=invscalefunc
            self._scaling = 'custom'
        else:
            newf = lambda x:x
            newif = newf
            self._scaling = 'linear'
            
        if self._scalefunc != newf or self._invscalefunc != newif :
            self._scalefunc = newf
            self._invscalefunc = newif    
            self.activateRange(self._rng)
        
    def setScalingLinear(self):
        self.setScaling(None)
        
    def setScalingLog(self,lower=0,upper=3,base=10.0):
        """
        linearly maps the values to a range from 10^lower to 10^upper and applies
        logarithmic scaling of the appropriate base to the result
        
        if upper is None, just offsets to the lower range
        if lower is None, logarithmic scaling is applied directly
        """
        from math import e,log
        if base == 10.0: 
            logfunc=np.log10
            denom=1
        elif base == e:
            logfun=np.log
            denom=1
        else:
            logfunc=np.log
            denom=log(base)
        
        if lower is None:
            self._scalefunc=lambda x: logfunc(x)/denom
            self._invscalefunc=lambda x:base**x
        elif upper is None:
            lower = base**lower
            xmin=self._linearStats['min']
            self._scalefunc=lambda x: logfunc(x-xmin+lower)/denom
            self._invscalefunc=lambda x: base**(x)+xmin-lower
        else:
            upper = base**upper
            lower = base**lower
            xmin,xmax=self._linearStats['min'],self._linearStats['max']
            range = (upper-lower)/(xmax-xmin)
            
            self._scalefunc=lambda x: logfunc((x-xmin)*range+lower)/denom
            self._invscalefunc=lambda x:(((base**x)-lower)/range)+xmin
        
        self.activateRange(self._rng)
        self._scaling = 'log'
        
    def setScalingSurfaceBrightness(self):
        """
        set the scaling so that output is magnitudes per square arcsecond.  Note
        that the zeropoint and pixelscale properties must be set for this scaling
        to work.
        
        note that if any part of the image is negative, errors can occur
        """        
        if self.zeropoint is None:
            raise ValueError('no zero point set')
        if self.pixelscale is None:
            raise ValueError('No pixel scale set')
        
        zpt,A=self.zeropoint,self.pixelscale[0]*self.pixelscale[1]
        offset = 2.5*np.log10(A) - zpt
            
        self._scalefunc = lambda im:-2.5*np.log10(im)+offset
        self._invscalefunc = lambda imsb:10**((imsb-offset)/-2.5)
        
        self.activateRange(self._rng)
        self._scaling = 'sb'
        
        
    def setScalingExp(self,base=10):
        from math import log
        logbase=log(base)
        self._scalefunc=lambda x:base**x
        self._invscalefunc=lambda x:np.log(x)/logbase
        
        self.activateRange(self._rng)
        self._scaling = 'exp'
        
    def setScalingPower(self,power=2):
        self._scalefunc=lambda x:x**power
        self._invscalefunc=lambda x:x**(1/power)
        
        self.activateRange(self._rng)
        self._scaling = 'power'
        
    def setScalingASinh(self):
        #TOOD:test
        self._scalefunc = np.arcsinh
        self._invscalefunc = np.sinh
        
        self.activateRange(self._rng)
        self._scaling = 'asinh'
    
    def _getScaling(self):
        return self._scaling
        
    scaling = property(_getScaling,setScaling)
        
    def clipOutliers(self,limits=(1,99),percentage=True,action='noclipmedian'):
        """
        clip points lying outside the range specified in the 2-sequence limits
        
        if percentage is outside True, the limits are percentages, otherwise 
        fractions
        
        action is a replacement action (see clipSigma)
        """
        if np.isscalar(limits):
            limits=(-limits,limits)
        elif len(limits)==2:
            pass
        else:
            raise ValueError('unrecognized limit form')
        
        im=self._active
        
        if percentage:
            mi,ma=np.min(im),np.max(im)
            rng=ma-mi
            limits=(mi+limits[0]*rng/100,mi+limits[1]*rng/100)
            
        clipcond=np.logical_or(im>limits[1],im<limits[0])
        self._active,nclip=self._repl_inds(action,clipcond)
        
        self._changed = True
        return nclip
        
    def clipSigma(self,sigma=12,fullsig=True,action='noclipmedian'):
        """
        if fullsig is true, the whole image is used for computing the standard
        deviation instead of just the active region
        
        action can be a value to set the outliers to, 'mean','median',
        'noclipmean','noclipmedian', 'fullmean', or 'fullmedian', or 'mask'
        
        returns number of pixels replaced
        """
        im=self._active
            
        if fullsig:
            statdat=self._scalefunc(self._extractArray(None))
        else:
            statdat=self._scalefunc(self._active)
            
        slim=sigma*self._fullStats['std']
        mean=self._fullStats['mean']
        
        wcond=np.logical_or(im>(mean+slim),im<(mean-slim))
        self._active,nclip=self._repl_inds(action,wcond)
        
        self._changed = True
        return nclip
    
    def clipInvalids(self,action='mask'):
        """
        clip points that are non-finite (i.e. NaN or infinite)
        
        action is a replacement action (see clipSigma)
        """
        wcond=np.logical_not(np.isfinite(self._active))
        self._active,nclip= self._repl_inds(action,wcond)
        
        self._changed = True
        return nclip
    
    def clipThreshold(self,thresh=0,action='mask',comp='<'):
        """
        clip data points that are above/below a threshold
        
        thresh is the threshold
        
        comp can be '<', '>','<=', or '>=' for clipping points that are 
        less that, greater than, less than or equal to, or greater than or
        equal to the threshold.
        
        action is a replacement action (see clipSigma)
        """
        
        if comp == '<':
            m = self._active < thresh
        elif comp == '>':
            m = self._active > thresh
        elif comp == '>=':
            m = self._active >= thresh
        elif comp == '<=':
            m = self._active <= thresh
        else:
            raise ValueError('unrecognized comp operator')
        self._active,nclip = self._repl_inds(action,m)
        
        self._changed = True
        return nclip
    
    def offsetData(self,offset='minp1'):
        """
        apply a global offset to all the current data
        
        offset can either by a value or 'gmin' to offset the data by the global 
        minimum or 'gminp1' to offset the lowest value to 1, or 'min' or 'minp1'
        to use the selected area minimum
        """
        
        
        if offset == 'gmin':
            offset = self._fullStats['min']
        elif offset == 'gminp1':
            offset = self._fullStats['min'] + 1
        elif offset == 'min':
            offset = np.min(self._active)
        elif offset == 'minp1':
            offset = np.min(self._active) + 1
        self._active-=offset
        self._changed = True
        
    def _repl_inds(self,action,cliparr):
        im=self._active
        
        if action == 'mask':
            im = np.ma.masked_array(im,cliparr,copy=False)
        elif action == 'median':
            im[cliparr]=np.median(im,None)
        elif action == 'fullmedian':
            im[cliparr]=self._fullStats['median']
        elif action == 'noclipmedian':
            im[cliparr]=np.median(im[~cliparr],None)
        elif action == 'mean':
            im[cliparr]=np.mean(im,None)
        elif action == 'fullmean':
            im[cliparr]=self._fullStats['mean']
        elif action == 'noclipmean':
            im[cliparr]=np.mean(im[~cliparr],None)
        elif type(action) != str and np.isscalar(action):
            im[cliparr]=action
        else:
            raise ValueError('unrecognized action')
        
        return im,sum(cliparr)
    
    @property
    def _fullStats(self):
        """
        statistics for the full image 
        """
        if not self._fstatd:
            v = self._scalefunc(self._extractArray(None).ravel())
            #TODO:deal with clipping/masks
            
            self._fstatd['median']=np.median(v)
            self._fstatd['mean']=np.mean(v)
            self._fstatd['std']=np.std(v)
            self._fstatd['min']=np.min(v)
            self._fstatd['max']=np.max(v)
            
        return self._fstatd
    
    @property
    def _linearStats(self):
        """
        statistics for the full image with linear scaling
        """
        if not self._lstatd:
            v = self._extractArray(None).ravel()
            #TODO:deal with clipping/masks
            
            self._lstatd['median']=np.median(v)
            self._lstatd['mean']=np.mean(v)
            self._lstatd['std']=np.std(v)
            self._lstatd['min']=np.min(v)
            self._lstatd['max']=np.max(v)
            
        return self._lstatd
        
    def plotImage(self,valrange='p99',flipaxis='y',invert=False,cb=True,clickinspect=True,clf=True,**kwargs):
        """
        This plots the associated active range of the image using matplotlib
        
        valrange can be:
        a 2-tuple with (lower,upper) range to display
        's#'/'sigma#'/'sig#' to use the specified S.D.s from the median
        'n#' to ignore the highest and lowest n values
        'p##.#' to use the central ## percent of the values for ranging
        (p or n can be ##,## to give lower,upper bounds)
        (if 'g' appears in any of the above codes, the global value is used)
        'i#,#,#,#' to ignore the specified values in the calculation of the range
        
        
        if clickinspect is True, will return the cid of the matplotlib event
        
        kwargs are passed into matplotlib imshow
        """
        import matplotlib.pyplot as plt
        
        vals=self._active
            
        if self.__examcid is not None:
            plt.gcf().canvas.mpl_disconnect(self.__examcid)
        if clf:
            plt.clf()
        if valrange:
            #if 'sig' in valrange or 'sigma' in valrange:
            if 's' in valrange:
                if 'g' in valrange:
                    #valrange=float(valrange.replace('g','').replace('sigma','').replace('sig',''))*self._fullStats['std']
                    valrange=float(valrange.replace('g','').replace('s',''))*self._fullStats['std']
                    vm=self._fullStats['median']
                else:
                    #valrange=float(valrange.replace('sigma','').replace('sig',''))*np.std(vals)
                    valrange=float(valrange.replace('s',''))*np.std(vals)
                    vm=np.median(vals)
                valrange=(vm-valrange,vm+valrange)
                
            elif 'n' in valrange or 'p' in valrange:
                if 'g' in valrange and valrange.replace('g','') in self._fstatd: #use cached value for this valrange if present
                    valrange = self._fstatd[valrange.replace('g','')]
                else:
                    if 'g' in valrange:
                        statvals = self._scalefunc(self._extractArray(None).ravel())
                    else:
                        statvals = vals.ravel()
                    vrspl=valrange.split(',')
                    if len(vrspl)==2:
                        if 'p' in valrange:
                            niglow=(100-float(vrspl[0].replace('g','').replace('p','')))*statvals.size/200
                            nigup=(100-float(vrspl[1].replace('g','').replace('p','')))*statvals.size/200
                        else:
                            niglow=float(vrspl[0].replace('g','').replace('n',''))
                            nigup=float(vrspl[1].replace('g','').replace('n',''))
                        niglow,nigup=int(round(niglow)),int(round(nigup))
                    elif len(vrspl)==1:
                        if 'p' in valrange:
                            nignore=(100-float(valrange.replace('g','').replace('p','')))*statvals.size/200
                        else:
                            nignore=float(valrange.replace('g','').replace('n',''))
                            
                        niglow=nigup=int(round(nignore))
                    else:
                        raise ValueError('unrecognized valrange w/p or n')
                    
                    sortval=statvals[np.argsort(statvals)]
                    
                    if niglow < 0:
                        niglow=0
                    if nigup < 0:
                        nigup=0
                        
                    if niglow+nigup >= sortval.size:
                        from warnings import warn
                        warn('ignored all of the values - displaying all instead')
                        valrange = (sortval[0],sortval[-1])
                    else:
                        if 'g' in valrange:
                            self._fstatd[valrange.replace('g','')] =  valrange = (sortval[niglow],sortval[-1 if nigup == 0 else -nigup])
                        else:
                            valrange = (sortval[niglow],sortval[-1 if nigup == 0 else -nigup])
                    
                    
                    
            elif 'i' in valrange:
                igvals=[float(vstr) for vstr in valrange.replace('i','').split(',')]
                m=np.ones(vals.shape,bool)
                for igval in igvals:
                    m = m & (vals != igval)
                valrange = (np.min(vals[m]),np.max(vals[m]))
            elif len(valrange) == 2:
                pass #correct form already
            elif np.isscalar(valrange):
                valrange=(0,valrange)
            else:
                raise ValueError('unrecognized valrange')
            from matplotlib.colors import Normalize
            nrm=Normalize(valrange[0],valrange[1])
            if invert:
                nrm.vmax,nrm.vmin = -1*nrm.vmin,-1*nrm.vmax
        else:
            nrm=None
        
        preint=plt.isinteractive()    
        try:
            plt.ioff()
            plt.imshow(-vals.T if invert else vals.T,norm=nrm,**kwargs)
            if cb:
                plt.colorbar()
            if 'x' in flipaxis:
                plt.xlim(*(plt.xlim()[::-1]))
            if 'y' in flipaxis:
                plt.ylim(*(plt.ylim()[::-1]))
        finally:
            if preint:
                plt.ion()
                plt.draw()
                
        if clickinspect:
            def onclick(event):
                if event.inaxes and event.button>2:
                    x,y=event.xdata,event.ydata
                    print x,y,vals[y,x] #TODO:check if there is a 1-pixel offset
                    
            self.__examcid=plt.gcf().canvas.mpl_connect('button_release_event',onclick)
        
        
    def plotHist(self,perc=None,gauss=None,clf=True,**kwargs):
        """
        kwargs go into hist
        reset defaults: 'log'->True,'bins'->100 or 25,'histtype'->'step'
        """
        import matplotlib.pyplot as plt
        
        vals=self._active
        
        if clf:
            plt.clf()
            
        if 'log' not in kwargs:
            kwargs['log'] = True
        if 'bins' not in kwargs:
            kwargs['bins'] = 100 if vals.size > 1e5 else 10
        if 'histtype' not in kwargs:
            kwargs['histtype'] = 'step'
            
        res=plt.hist(vals.ravel(),**kwargs)
        if perc:
            if isinstance(perc,basestring):
                perc=[float(p) for p in perc.split(',')]
            if np.isscalar(perc):
                perc=(perc,perc)
            elif len(perc) !=2:
                raise ValueError('unrecognized perc')
            mi,ma=np.min(vals),np.max(vals)
            rng=perc[0]*(ma-mi)/100
            yl=plt.ylim()
            plt.vlines([ma-rng,mi+rng],1 if kwargs['log'] else 0,np.max(res[0]),color='k',linestyle='--')
            plt.ylim(*yl)
            #plt.axvline(ma-rng,c='k',ls='--')
            #plt.axvline(mi+rng,c='k',ls='--')
        if gauss:
            gf=np.std(vals.ravel())*np.random.randn(vals.size)+np.median(vals.ravel())
            kwargs['ec']='g'
            plt.hist(gf,**kwargs)
            
        
    def getStats(self):
        """
        return mean,median,stddev,min, and max for the scaled active region
        """
        im=self._active
        v=im.ravel()
        
        return {'mean':np.mean(v),'median':np.median(v),'std':np.std(v),'min':np.min(v),'max':np.max(v)}
    
    @property
    def shape(self):
        """
        tuple with dimensions of the currently selected image
        (the internal representation is flipped for FITS - this is the correct orientation)
        """
        return self._extractArray(None).shape
    
    @property
    def size(self):
        """
        number of pixels in the currently selected image
        """
        return self._extractArray(None).size
    
    def _getPixScale(self):
        return self._pixscale
    def _setPixScale(self,val):
        if np.isscalar(val):
            self._pixscale = (val,val)
        else:
            if len(val) != 2:
                raise ValueError('set pixel scale as (xscale,yscale)')
            self._pixscale = tuple(val)
    pixelscale = property(_getPixScale,_setPixScale,doc="""
    Pixel scale of the image in 
    """)
            
    def _getZpt(self):
        return self._zpt
    def _setZpt(self,val):
        if not np.isscalar(val):
            return ValueError('zero points must be scalars')
        self._zpt = val
    zeropoint = property(_getZpt,_setZpt,doc="""
    The photometric zero-point for this CCDImage, to be used in other 
    processing tasks (see astropysics.phot)
    """)
    
    def _getData(self):
        return self._active
    data = property(_getData,doc="""
    direct access to the array holding the current data.  If this data is 
    subsequently edited, applyChanges() should be called to ensure consistancy
    """)
    
    def _getRange(self):
        return self._rng
    def _setRange(self,value):
        self.activateRange(value)
    range = property(_getRange,_setRange,doc="""
    returns or sets the range (setting is syntactic sugar for activateRange)
    """)
    
            
class FitsImage(CCDImage):
    def __init__(self,fn,range=None,scaling=None,hdu=0,memmap=0):
        import pyfits
        
        fnl = fn.lower()
        if fnl.endswith('fit') or fnl.endswith('fits'):
            self.fitsfile = pyfits.open(fn,memmap=memmap)
        else:
            raise IOError('unrecognized file type')
        
        self._chdu = None
        self._newhdu = hdu 
        
        CCDImage.__init__(self,range=range,scaling=scaling)
        #self._updateFromHeader()#called implicity in CCDImage.__init__ through activateRange
        
    def close(self):
        self.fitsfile.close()
    
    def __del__(self):
        if self.fitsfile is not None:
            self.fitsfile.close()
    
    def _extractArray(self,range):
        """
        kwarg hdu chooses which hdu to use for this image
        """
        if self._newhdu != self._chdu:
            
            oldhdu,oldfstatd,oldlstatd = self._chdu,self._fstatd,self._lstatd
            self._fstatd.clear()
            self._lstatd.clear()
            try:
                self._chdu = self._newhdu
                res = self._extractArray(range)
            except IndexError:
                print 'New HDU Has incompatible range - using whole image'
                res = self._extractArray(None,hdu=None)
                self._rng = None
            except:
                self._fstatd = oldfstatd
                self._lstatd = oldlstatd
                self._chdu = oldhdu
                raise
            
            self._updateFromHeader()
            return res
                
        if range is None:
            im = self.fitsfile[self._chdu].data.T
        else:            
            ny,nx = self.fitsfile[self._chdu]._dimShape()
            xl,xu,yl,yu = range
            if xl < 0 or xu > nx or yl < 0 or yu > ny:
                raise IndexError('Attempted range %i,%i;%i,%i on image of size %i,%i!'%(xl,xu,yl,yu,nx,ny))
            im = self.fitsfile[self._chdu].data.T[xl:xu,yl:yu]
            
        return im
    
    
    def _applyArray(self,range,imdata):
        if self._invscalefunc is self._scalefunc: #linear mapping so changes already applied
            pass
        else:
            if range is None:
                self.fitsfile[self._chdu].data = imdata.T
            elif len(range) == 4:
                xl,xu,yl,yu = range
                ny,nx = self.fitsfile[self._chdu].shape #fits files are flipped
                if xl < 0 or xu > nx or yl < 0 or yu > ny:
                    raise IndexError('Attempted range %i,%i;%i,%i on image of size %i,%i!'%(xl,xu,yl,yu,nx,ny))
                self.fitsfile[self._chdu].data[yl:yu,xl:xu] = imdata.T
            else:
                raise ValueError('Unregonized form for range')
        
    def _updateFromHeader(self):
        d = dict(self.fitsfile[self._chdu].header.items())
        
        d.pop('',None)
        
        if 'PIXSCALE' in d:
            self.pixelscale = d['PIXSCALE']
        elif 'PIXSCAL1' in d:
            if 'PIXSCAL2' in d:
                self.pixelscale = (d['PIXSCAL1'],d['PIXSCAL2'])
            else:
                self.pixelscale = d['PIXSCAL1']
                
        if 'ZEROPOINT' in d:
            self.zeropoint=d['ZEROPOINT']
        elif 'ZEROPT' in d:
            self.zeropoint=d['ZEROPT']
    
    def _getHdu(self):
        return self._chdu
    def _setHdu(self,value):
        self._newhdu = value
        currrng = self._rng
        self.activateRange(currrng)
    hdu = property(_getHdu,_setHdu,doc="""
    the hdu of the current fits file
    """)       
    
    
def load_image_file(fn,**kwargs):
    """
    Factory function for generating CCDImage objects from files - the exact 
    class type will be inferred from the extension.  kwargs will be
    passed into the appropriate constructor
    
    Supports:
    *FITS files
    """
    from os import path
    ext = path.splitext(fn)[-1].lower()[1:]
    if ext =='fits' or ext == 'fit':
        return FitsImage(fn,**kwargs)
    else:
        raise ValueError('Unrecognized file type for file '+fn)
    
        
def mosaic_objects(xcens,ycens,radii,images,row=None,titles=None,noticks=True,
                   clf=True,logify=False,imdict=None,**kwargs):
    """
    generate mosaics of objects from different images
    xcens,ycens,and radii are the x,y of the center and the radius to plot
    images are the image filenames
    row determines which row of the mosaic the object is placed in.  If None,
    an optimal square/rectangular mosaic will be used
    titles is a sequence of strings to title each object
    if noticks is True, x and y ticks are not drawn
    clf determines if the figure is cleared before the mosaic is generated
    imdict is a dictionary with values that are CCDImage objects.  Any image 
    filenames will be placed into the provided dictionary keys with values
    generated by calling load_image_files on the filenames
    """
    import matplotlib.pyplot as plt
    from operator import isSequenceType
    
    if isSequenceType(xcens):
        onev=np.ones(len(xcens))
    elif isSequenceType(ycens):
        onev=np.ones(len(ycens))
    elif isSequenceType(radii):
        onev=np.ones(len(radii))
    elif isSequenceType(images) and type(images) is not str:
        onev=np.ones(len(images))
    elif isSequenceType(row):
        onev=np.ones(len(row))
    else:
        onev=np.ones(1)
        
    if not isSequenceType(xcens):
        xcens=onev*xcens
    if not isSequenceType(ycens):
        ycens=onev*ycens
    if not isSequenceType(radii):
        radii=onev*radii
    if row is not None and not isSequenceType(row):
        row=onev*row
    if isinstance(images,basestring) or not isSequenceType(images):
        images=[images for i in onev]
        
    if not titles:
        titles=['' for i in onev]
    
    if row is None:
        n=len(onev)
        cols=int(np.ceil(n**0.5))
        rws,row=0,[]
        while rws*cols<n:
            for i in range(cols):
                row.append(rws+1)
            rws+=1
        row=row[:n]
            
    else:
        rws=np.max(row)
        cols=dict(((r,0) for r in row))
        for r in row:
            cols[r]+=1
        cols=max(cols.values())
        row = row.astype(int)
    if imdict is None:
        imdict={}
    
    subps=[]
    rowcounts=[0 for i in range(rws)]
    for x,y,r,im,rw,ti in zip(xcens,ycens,radii,images,row,titles):
        #determine which subplot index to use ensuring that everything ends up in the right row
        i=rowcounts[rw-1]+(rw-1)*cols+1
        rowcounts[rw-1]+=1
        subps.append((x,y,r,ti,i,im))
    
    for img in images:
        if img not in imdict:
            imdict[img]=load_image_file(img)    
        
    preint=plt.isinteractive()    
    try:
        plt.ioff()
        if clf:
            plt.clf()
            
        for x,y,r,ti,i,im in subps:
            im=imdict[im]
            if logify:
                im.setScalingLog()
            plt.subplot(rws,cols,i)
            im.activateRange((x,y,r))
            kwargs['cb']=False
            kwargs['clf']=False
            im.plotImage(**kwargs)
            if noticks:
                plt.xticks([])
                plt.yticks([])
            plt.title(ti)
        plt.show()
        plt.draw()
    finally:
        if preint:
            plt.ion()
        
    
def kcorrect_images(images,bands,z,range=None,zeropoints=None,pixelareas=None,retdict=False,offset=None,**kckwargs):
    """
    This function performs kcorrections pixel-by-pixel on a matched set of 
    images.  Note that one must be carefule to ensure the images are matched
    (e.g. they must be degraded to the largest PSF image), something that this 
    function does nothing about.  
    See astropysics.phot.kcorrect for more details on the k-correction algorithm
    
    images are a sequence of CCDImage objects, or a list of filenames that will 
    be loaded as images (using the 0th hdu if FITS).
    
    bands should be the photometric bands of each of the images as accepted by
    kcorrect
    
    z is the redshift to correct this image to
    
    range is the range from each of the images to select
    
    zeropoints and pixelareas (sq asec) will be deduced from the CCDImage 
    properties zeropoint and pixelscale if they are not provided (if they are,
    they will everwrite the existing properties)
    
    offset is the type of offset to be used (passed into offsetData)
    
    extra kwargs will be passed into astropysics.phot.kcorrect
    
    returns absmag,kcorrection,chi2 as arrays shaped like the input range
    (for absmag and kcorrection, the first dimension is the bans)
    
    if True, retdict means the absmag and kcorrection will be returned as
    dictionaries of arrays with the band names as the keys
    """
    from .phot import kcorrect
    nbands = len(bands)
    if  len(images) != nbands:
        raise ValueError("images and bands don't match!")
    if zeropoints is not None and nbands != len(zeropoints):
        raise ValueError("zeropoints and # of bands don't match!")
    if pixelareas is not None and nbands != len(pixelareas):
        raise ValueError("pixelscale and # of bands don't match!")
    
    imobj,imdata = [],[]
    for i,im in enumerate(images):
        if isinstance(im,CCDImage):
            imobj.append(im)
        elif isinstance(im,basestring):
            imobj.append(load_image_file(im))
        else:
            raise ValueError('image #%i is not a CCDImage or string'%i)
        
    if zeropoints is not None:
        for zpt,im in zip(zeropoints,imobj):
            im.zeropoint = zpt
    if pixelareas is None:
        pixelareas=[]
        for pa,im in zip(pixelareas,imobj):
            im.pixelscale = pa**0.5
            
    for i,im in enumerate(imobj):
        im.offsetData(offset)
        im.activateRange(range)
        im.setScaling('sb')
        imdata.append(im.data)
        
    dshape = imdata[0].shape
        
    sbs = [da.ravel() for da in imdata]
    
    #TODO: implement errors
    #TODO: shape matching tests?
    ams,kcs,chi2s = kcorrect(sbs,z*np.ones(imdata[0].size),filterlist=bands) 
    
    targetshape = tuple(np.r_[nbands,dshape]) #all should match
    
    ams = ams.reshape(targetshape)
    kcs = kcs.reshape(targetshape)
    chi2s = chi2s.reshape(dshape)
    
    if retdict:
        dams = dict([(b,am) for am,b in zip(ams,bands)])
        dkcs =  dict([(b,kc) for kc,b in zip(kcs,bands)])
        return dams,dkcs,chi2s
    else:
        return ams,kcs,chi2s