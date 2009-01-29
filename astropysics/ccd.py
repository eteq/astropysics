#Erik Tollerud (etolleru@uci.edu) 2008

"""
This module contains objects and functions for viewing and/or reducing CCD 
images.

This package deals mostly with the raw CCD images themselves - for science tools
see the phot and spec packages.
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

class CCDImage(object):
    """
    Represents a CCD image.  Currently only fits format is supported.
    """
    def __init__(self,fn,irange=None,scaling=None,ihdu=0,memmap=0):
        import pyfits
        fnl = fn.lower()
        if fnl.endswith('fit') or fnl.endswith('fits'):
            self.file = pyfits.open(fn,memmap=memmap)
        else:
            self.file = None
            raise IOError('unrecognized file type')
        
        
        self.chdu = ihdu
        self._rng = None
        self.__examcid = None
        self._changed = False
        self.applyChangesOnActivate = False#True #TODO:check this out for consistency
        self.setScaling(scaling)
        self._calcGlobalStats()
        
        self.activateRange(irange) 
        
        
        
    def __del__(self):
        if self.file is not None:
            self.file.close()
        
    def close(self):
        self.file.close()
        
    def setHDU(self,hdu):
        self.file[hdu]
        self.chdu=hdu
        
    def setScaling(self,scalefunc,invscalefunc=None):
        """
        scalefunc must be a callable or evaluate to False (linear mapping)
        
        alternatively, it may be a string specifying one of the builtin types
        of scaling
        
        if inverse function is None, changes cannot be saved
        """
        if scalefunc:
            if type(scalefunc) == str:
                sfl=scalefunc.lower()
                if sfl=='linear':
                    return self.setScalingLinear()
                elif sfl=='log':
                    return self.setScalingLog()
                elif sfl=='exp':
                    return self.setScalingExp()
                elif sfl=='power':
                    return self.setScalingPower()
                elif sfl=='asinh':
                    return self.setScalingASinh()
                else:
                    raise ValueError('Unrecognized scaling string')
            elif not callable(scalefunc):
                raise ValueError('scalefunc is not a callable')
            self._scalefunc=scalefunc
            self._invscalefunc=invscalefunc
        else:
            self._scalefunc=lambda x:x
            self._invscalefunc=self._scalefunc
        
            
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
            xmin=self._global_min
            self._scalefunc=lambda x: logfunc(x-xmin+lower)/denom
            self._invscalefunc=lambda x: base**(x)+xmin-lower
        else:
            upper = base**upper
            lower = base**lower
            xmin,xmax=self._global_min,self._global_max
            range = (upper-lower)/(xmax-xmin)
            
            self._scalefunc=lambda x: logfunc((x-xmin)*range+lower)/denom
            self._invscalefunc=lambda x:(((base**x)-lower)/range)+xmin
        
        self.activateRange(self._rng)
        
    def setScalingExp(self,base=10):
        from math import log
        logbase=log(base)
        self._scalefunc=lambda x:base**x
        self._invscalefunc=lambda x:np.log(x)/logbase
        
    def setScalingPower(self,power=2):
        self._scalefunc=lambda x:x**power
        self._invscalefunc=lambda x:x**(1/power)
        
    def setScalingASinh(self):
        #TOOD:test
        self._scalefunc=np.asinh
        self._invscalefunc=np.sinh
        
        
    def activateRange(self,range):
        """
        chooses the range to activate for further operations
        
        can either be None (whole file), a 4-tuple of
        the form (xl,xu,yl,yu) or a 3-tuple of the form (xcen,ycen,radius)
        """
        if self.applyChangesOnActivate and self._changed:
            try:
                self.applyChanges()
            except:
                from warnings import warn
                warn("No inverse function available - can't save",category = RuntimeWarning)
                
        if range is None:
            im = self.file[self.chdu].data
        else:
            if len(range) == 3:
                xcen,ycen,rad = range
                range=(ycen-rad,ycen+rad,xcen-rad,xcen+rad)
            elif len(range) == 4: 
                xl,xu,yl,yu=range
                range=(yl,yu,xl,xu)
            else:
                raise ValueError('Unregonized form for range')
            
            xl,xu,yl,yu = range
            im = self.file[self.chdu].data
            xr,yr=im.shape
            
            #do checks
            if xu <= xl:
                if xu == xl:
                    raise ValueError('tried to activate size 0 range')
                xu,xl=xl,xu
            if yu <= yl:
                if yu == yl:
                    raise ValueError('tried to activate size 0 range')
                yu,yl=yl,yu
            #ensure range does not go off edge - recenter    
            if xu-xl > xr:
                xl,xu=0,xr
            elif xl < 0:
                xl,xu=0,xu-xl
            elif xu > xr:
                xl,xu=xl-xu,xr
            if yu-yl > yr:
                yl,yu=0,yr
            elif yl < 0:
                yl,yu=0,yu-yl
            elif yu > yr:
                yl,yu=yl-yu,yr
                
                
            im = im[xl:xu,yl:yu]
        
        self._rng = range
        self._active = self._scalefunc(im)
        
    def getRange(self,editing = True):
        """
        editing = whether or not you intend to change the data
        """
        self._changed = editing
        return self._rng
    
    def applyChanges(self):
        if not self._invscalefunc:
            raise Exception('No inverse function available')
        if self._invscalefunc is self._scalefunc: #linear mapping so changes already applied
            pass
        else:
            altim = self._invscalefunc(self._active)
            range = self._rng
            
            if range is None:
                self.file[self.chdu].data = altim
            elif len(range) == 3:
                xcen,ycen,rad = range
                self.file[self.chdu].data[xcen-rad:xcen+rad,ycen-rad:ycen+rad]  = altim
            elif len(range) == 4:
                xl,xu,yl,yu = range
                self.file[self.chdu].data[xl:xu,yl:yu] = altim
            else:
                raise ValueError('Unregonized form for range')
            
        self._calcGlobalStats()
        self._changed = False
    
    def clipOutliers(self,limits=(1,99),percentage=True,action='noclipmedian'):
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
        clipi=self._repl_inds(action,clipcond)
        
        self._changed = True
        return len(clipi[0])
        
    def clipSigma(self,sigma=12,fullsig=True,action='noclipmedian'):
        """
        if fullsig is true, the whole image is used for computing the standard
        deviation instead of just the active region
        
        action can be a value to set the outliers to, 'mean','median',
        'noclipmean','noclipmedian', 'fullmean', or 'fullmedian'
        
        returns number of pixels replaced
        """
        im=self._active
            
        if fullsig:
            statdat=self._scalefunc(self.file[self.chdu].data)
        else:
            statdat=self._scalefunc(self._active)
            
        slim=sigma*statdat.std()
        mean=statdat.mean()
        
        wcond=np.logical_or(im>(mean+slim),im<(mean-slim))
        clipi=self._repl_inds(action,wcond)
        
        self._changed = True
        return len(clipi[0])
    
    def clipInvalids(self,action):
        wcond=np.logical_not(np.isfinite(self._active))
        clipi=self._repl_inds(action,wcond)
        
        self._changed = True
        return len(clipi[0])
        
    def _repl_inds(self,action,cliparr):
        clipi=np.where(cliparr)
        nclipi=np.where(np.logical_not(cliparr))
        im=self._active
        
        if action == 'median':
            im[clipi]=np.median(im,None)
        elif action == 'fullmedian':
            im[clipi]=self._global_median
        elif action == 'noclipmedian':
            im[clipi]=np.median(im[nclipi],None)
        elif action == 'mean':
            im[clipi]=np.mean(im,None)
        elif action == 'fullmean':
            im[clipi]=self._global_mean
        elif action == 'noclipmean':
            im[clipi]=np.mean(im[nclipi],None)
        elif type(action) != str and np.isscalar(action):
            im[clipi]=action
        else:
            raise ValueError('unrecognized action')
        
        return clipi
    
    def _getGlobalData(self):
        return self.file[self.chdu].data.ravel()
    
    def _calcGlobalStats(self):
        v = self.file[self.chdu].data.ravel()
        
        self._global_median=np.median(v)
        self._global_mean=np.mean(v)
        self._global_std=np.std(v)
        self._global_min=np.min(v)
        self._global_max=np.max(v)
        
    def showImage(self,valrange='p99',flipaxis='y',invert=False,cb=True,clickinspect=True,clf=True,**kwargs):
        """
        valrange can be:
        a 2-tuple with (lower,upper) range to display
        'sigma#'/'sig#' to use the specified S.D.s from the median
        'n#' to ignore the highest and lowest n values
        'p##.#' to use the central ## percent of the values for ranging
        (p or n can be ##,## to give lower,upper bounds)
        'i#,#,#,#' to ignore the specified values in the calculation of the range
        
        if clickinspect is True, will return the cid of the matplotlib event
        
        kwargs are passed into imshow
        """
        
        vals=self._active
            
        if self.__examcid is not None:
            plt.gcf().canvas.mpl_disconnect(self.__examcid)
        if clf:
            plt.clf()
        if valrange:
            if 'sig' in valrange or 'sigma' in valrange:
                valrange=float(valrange.replace('sigma','').replace('sig',''))*self._global_std
                vm=np.median(vals,None)
                valrange=(vm-valrange,vm+valrange)
                
            if 'n' in valrange or 'p' in valrange:
                vrspl=valrange.split(',')
                if len(vrspl)==2:
                    if 'p' in valrange:
                        niglow=(100-float(vrspl[0].replace('p','')))*vals.size/200
                        nigup=(100-float(vrspl[1].replace('p','')))*vals.size/200
                    else:
                        niglow=float(vrspl[0].replace('n',''))
                        nigup=float(vrspl[1].replace('n',''))
                    niglow,nigup=int(round(niglow)),int(round(nigup))
                elif len(vrspl)==1:
                    if 'p' in valrange:
                        nignore=(100-float(valrange.replace('p','')))*vals.size/200
                    else:
                        nignore=float(valrange.replace('n',''))
                        
                    niglow=nigup=int(round(nignore))
                else:
                    raise ValueError('unrecognized valrange w/p or n')
                
                
                
                sortval=vals.copy().ravel()
                sortval.sort()
                
                if niglow < 0:
                    niglow=0
                if nigup < 0:
                    nigup=0
                    
                if niglow+nigup >= sortval.size:
                    from warnings import warn
                    warn('ignored all of the values - displaying all instead')
                    valrange = (sortval[0],sortval[-1])
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
            plt.imshow(-vals if invert else vals,norm=nrm,**kwargs)
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
        
        
    def showHist(self,perc=None,gauss=None,clf=True,**kwargs):
        """
        kwargs go into hist
        reset defaults: 'log'->True,'bins'->100 or 25,'histtype'->'step'
        """
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
            if type(perc) == str:
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
    
    def getShape(self):
        return self.file[self.chdu]._dimShape()
    
    def getSize(self):
        return self.file[self.chdu].size()
        
        
def mosaic_objects(xcens,ycens,radii,images,row=None,titles=None,noticks=True,clf=True,logify=False,**kwargs):
    """
    generate mosaics of objects from different images
    xcens,ycens,and radii are the x,y of the center and the radius to plot
    images are the image filenames
    row determines which row of the mosaic the object is placed in.  If None,
    an optimal square/rectangular mosaic will be used
    titles is a sequence of strings to title each object
    if noticks is True, x and y ticks are not drawn
    clf determines if the figure is cleared before the mosaic is generated
    """
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
    if type(images) is str or not isSequenceType(images):
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
       
    imd=dict(((i,[]) for i in images))
    rowcounts=[0 for i in range(rws)]
    for x,y,r,im,rw,ti in zip(xcens,ycens,radii,images,row,titles):
        #determine which subplot index to use ensuring that everything ends up in the right row
        i=rowcounts[rw-1]+(rw-1)*cols+1
        rowcounts[rw-1]+=1
        imd[im].append((x,y,r,ti,i))
    
    preint=plt.isinteractive()    
    try:
        plt.ioff()
        if clf:
            plt.clf()
        for img,ts in imd.iteritems():
            im=CCDImage(img)
            if logify:
                im.setScalingLog()
            for x,y,r,ti,i in ts:
                plt.subplot(rws,cols,i)
                im.activateRange((x,y,r))
                kwargs['cb']=False
                kwargs['clf']=False
                im.showImage(**kwargs)
                if noticks:
                    plt.xticks([])
                    plt.yticks([])
                plt.title(ti)
            im.close()
        plt.show()
        plt.draw()
    finally:
        if preint:
            plt.ion()
        
    
    