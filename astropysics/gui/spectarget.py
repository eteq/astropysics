#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This module contains the internals for the SpecTarget gui.
"""
from __future__ import division,with_statement
import numpy as np

from enthought.traits.api import HasTraits,on_trait_change,Instance,Float,\
                                 String,ListStr,ListFloat,Button,Bool,\
                                 Property,Event,Array,Enum,TraitError
from enthought.traits.ui.api import View,Group,HGroup,HFlow,Item,Handler
from enthought.traits.ui.menu import ModalButtons
from enthought.chaco.api import Plot,ArrayPlotData,HPlotContainer,ColorBar,\
                            LinearMapper,jet,LassoOverlay,ColormappedSelectionOverlay
from enthought.chaco.tools.api import PanTool, ZoomTool,LassoSelection
from enthought.enable.component_editor import ComponentEditor

from ..phot import CMDAnalyzer

class MaskMaker(object):
    """
    base class for objects that make masks - subclasses should:
    
    * replace name with a name for the mask-making operation
    * shortlen: maximum length of shortname
    * override makeMask(self,filename)
    * size: approximate dimensions of a mask as (len,width) in degrees
    
    expectations for special priority groups:
    
    * -1: alignment stars
    * -2: guide stars
    * -9: remove from list unless guide/align star
    
    """
    name = 'Default mask-maker'
    shortlen = None
    size = None
    def makeMask(self,filename):
        raise NotImplementedError
    
    def __init__(self,spectargobj):
        self.sto = spectargobj

class ODHandler(Handler):
    def apply(self,info):
        o = info.object
        o.cmda.offsetbands = o.offsetbands
        o.cmda.offsetweights = 1.0/o.offsetscales
        o.cmda.locweight = 1.0/o.locscale if o.locscale != 0 else 0
        o.st.updatepri = True
        
    def closed(self,info, is_ok):
        if is_ok:
            self.apply(info)
        Handler.closed(self,info,is_ok)

class OffsetDialog(HasTraits):
    offsetbands = ListStr
    offsetscales = Array(float)
    locscale = Float
    
    view = View(Group(Item('offsetbands'),
                      Item('offsetscales'),
                      Item('locscale'),
                      layout='normal'),
                resizable=True,title='Offset Settings',handler=ODHandler(),
                buttons = ModalButtons,kind='modal')
    
    def __init__(self,cmda,st):
        self.cmda = cmda
        self.st = st
        self.offsetbands = cmda.offsetbands if cmda.offsetbands is not None else []
        self.offsetscales = 1.0/cmda.offsetweights
        self.locscale = 1.0/cmda.locweight if cmda.locweight != 0 else 0
            
class SpecTarget(HasTraits):
    distance = Float(10)
    distmod = Float(0)
    cmda = Instance(CMDAnalyzer)
    plotcontainer = HPlotContainer
    locplot = Plot
    cmplot = Plot
    cb = ColorBar
    data = ArrayPlotData
    locs = ArrayPlotData
    xax = String
    yax = String
    priband = String('')
    pribandw = Float(0)
    pripower = Float(2)
    groupcutoff = Float(100)
    agband = String
    astarcut = Float(18.5)
    astarcut2 = Float(20)
    gstarcut = Float(16)
    agigpri = Float(0.5)
    hideg0 = Bool(True)
    agigg1 = Bool(False)
    offsetset = Button
    priorities = Property(depends_on='distmod,offsetset')
    
    masktype = Enum(['DEIMOS'])
    maskmaker = Instance(MaskMaker) 
    maskfn = String('spectarg1')
    maskshortname = String('Nonam1')
    masklongname = String('')
    makemaskbutton = Button
    
    updatepri = Event
    updatedata = Event
    updatelocs = Event
    updategastars = Event
    
    view = View(Group(Item('plotcontainer',editor=ComponentEditor(size=(800,400)),show_label=False),
                      Group(Group(Item('distance',label='Distance (kpc)'),
                                  Item('distmod',label='Distance Modulus'),
                                  Item('xax',label='CMD x-axis'),
                                  Item('yax',label='CMD y-axis'),
                                  columns=2),
                            Group(Item('agband',label='Alignment/Guide Star Band'),
                                  Item('gstarcut',label='Guide Star Cutoff'),
                                  Item('astarcut',label='Alignment Star Cutoff'),
                                  Item('astarcut2',label='Alignment Star Cutoff 2'),
                                  Item('agigpri',label='A/G Ignore priority'),
                                  Item('agigg1',label='A/G Ignore if group 1?'),
                                  columns=6),
                            HGroup(Item('priband',label='Priority Band'),
                                   Item('pribandw',label='Band Priority Weight'),
                                   Item('offsetset',label='Offset Settings'),
                                   Item('pripower',label='Power for converting offsets into priorities'),
                                   Item('groupcutoff',label='Last Group Cutoff'),
                                   Item('hideg0',label='Set Group 0 priority to 0?')
                                   ),
                            Group(Item('masktype',label='Mask Type'),
                                  Item('maskfn',label='Mask File(base)'),
                                  Item('makemaskbutton',label='Make Mask'),
                                  Item('maskshortname',label='Mask Name(short)'),
                                  Item('masklongname',label='Mask Name(long)'),
                                  columns=3)),
                      layout='normal'),
                resizable=True,title='Spectra Targeter')
    
    
    def __init__(self,*args,**kwargs):
        """
        generates the SpecTarget GUI - args can be either:
        
        * SpecTarget(cmdanalyzer) 
            provide a ready-made CMDAnalyzer
        * SpecTarget(fiducialdict,data) 
            provide a dictionary of fiducial data and a set of data (array or
            dict)
        
        kwargs get applied to the CMDAnalyzer
        """
        if len(args) == 1:
            cmda = args[0]
        elif len(args) == 2:
            cmda = CMDAnalyzer(args[0].values(),args[0].keys())
            cmda.setData(args[1])
        #elif len(args) == 3:
        #    cmda = CMDAnalyzer(args[0].values(),args[0].keys())
        #    cmda.setData(args[1])
        #    cmda.locs = args[2]
        else:
            raise TypeError('SpecTarget initializer takes 1 or 2 arguments')
        for k,v in kwargs.iteritems():
            setattr(cmda,k,v)
        self.cmda = cmda
        
        self._distinner = False
        
        ob = cmda.offsetbands
        
        if ob is not None and len(ob) > 1:
            self.xax,self.yax = ob[0],ob[1]
        else:
            bi1,bi2 = cmda.validbandinds[0],cmda.validbandinds[1]
            self.xax,self.yax = cmda.bandnames[bi1],cmda.bandnames[bi2]
        self.agband = self.yax

        self.data = ArrayPlotData()
        self.locs = ArrayPlotData()
        self._update_x_data()
        self._update_y_data()
        self.updatelocs = True
        self.updatepri = True
        self.updategastars = True
                                  
        cmdplot = Plot(self.data,resizable='hv')
        objplotcmd = cmdplot.plot(('x','y','pri'),name='target_mags',type='cmap_scatter',marker='dot',marker_size=1,color_mapper=jet)[0]
        cmdplot.plot(('fidx','fidy'),name='fid_mags',type='line',color='black')
        cmdplot.plot(('gx','gy'),name='guide_mags',type='scatter',color='gray',marker='inverted_triangle',outline_color='red',marker_size=5)  
        cmdplot.plot(('ax','ay'),name='align_mags',type='scatter',color='black',marker='diamond',outline_color='red',marker_size=4)
        cmdplot.tools.append(PanTool(cmdplot,drag_button='right'))
        cmdplot.tools.append(ZoomTool(cmdplot)) 
        
        
        
        
        locplot = Plot(self.locs,resizable='hv')
        objplotloc = locplot.plot(('ra','dec','pri'),name='target_locs',type='cmap_scatter',marker_size=2,color_mapper=jet)[0]
        locplot.plot(('centerx','centery'),name='cen_locs',type='scatter',color='black',marker='cross',marker_size=10,line_width=4)
        locplot.plot(('boxx','boxy'),name='box_locs',type='line',color='black',line_width=2)
        locplot.plot(('gra','gdec'),name='guide_locs',type='scatter',color='gray',marker='inverted_triangle',outline_color='red',marker_size=5)
        locplot.plot(('ara','adec'),name='align_locs',type='scatter',color='black',marker='diamond',outline_color='red',marker_size=5)
        locplot.tools.append(PanTool(locplot,drag_button='right'))
        locplot.tools.append(ZoomTool(locplot))  
        
        cb = ColorBar(index_mapper=LinearMapper(range=cmdplot.color_mapper.range),color_mapper=cmdplot.color_mapper,resizable='v',width=30,padding=5)
        
        ls1 = LassoSelection(objplotcmd,selection_datasource=objplotcmd.color_data)
        objplotcmd.tools.append(ls1)
        objplotcmd.overlays.append(LassoOverlay(objplotcmd,lasso_selection=ls1))
        objplotcmd.active_tool = ls1
        ls1.on_trait_change(self._selection_changed,'selection_changed')
        ls1.on_trait_change(self._selection_completed,'selection_completed')
        objplotcmd.overlays.append(ColormappedSelectionOverlay(objplotcmd,selection_type='mask',selection_datasource=objplotcmd.color_data))
        
        ls2 = LassoSelection(objplotloc,selection_datasource=objplotloc.index)
        objplotloc.tools.append(ls2)
        objplotloc.overlays.append(LassoOverlay(objplotloc,lasso_selection=ls2))
        objplotloc.active_tool = ls2
        ls2.on_trait_change(self._selection_changed,'selection_changed')
        ls2.on_trait_change(self._selection_completed,'selection_completed')
        
        self.plotcontainer = HPlotContainer(use_backbuffer=True)
        self.locplot = locplot
        self.cmplot = cmdplot
        self.cb = cb
        self.plotcontainer.add(locplot)
        self.plotcontainer.add(cmdplot)
        self.plotcontainer.add(cb)
        
        self.on_trait_change(self._cmda_late_changed,'cmda')
        self.on_trait_change(self._xax_late_changed,'xax')
        self.on_trait_change(self._yax_late_changed,'yax')
        self.on_trait_change(self._late_dogastars,'agband,gstarcut,astarcut,agigpri')
        
        
        self.priband = self.yax #must be done after data is set
        self._cmda_late_changed() #clean things up in case random stuff changes
        self._masktype_changed(self.masktype) #call to intialize default mask maker
        
    def _selection_changed(self,obj,name,new):
        cds = (self.cmplot.plots['target_mags'][0].color_data,self.locplot.plots['target_locs'][0].color_data)
        datasource = obj.selection_datasource   
        targds = cds[1] if cds[0] == datasource else cds[0]
        
        mask = datasource.metadata['selection'].astype(bool)
        targds.set_mask(mask)
        
    def _selection_completed(self,obj,name,new):
        cds = (self.cmplot.plots['target_mags'][0].color_data,self.locplot.plots['target_locs'][0].color_data)
        datasource = obj.selection_datasource   
        targds = cds[1] if cds[0] == datasource else cds[0]
        
        targds.set_mask(None)
        
    def _maskfn_changed(self):
        self.dsimready = False

    def _distance_changed(self,old,new):
        if self._distinner:
            self._distinner = False
        else:
            self.cmda.distance = self.distance
            self._distinner = True
            self.distmod = self.cmda.distmod
            self.updatedata = True
        
    def _distmod_changed(self,old,new):
        if self._distinner:
            self._distinner = False
        else:
            self.cmda.distmod = self.distmod
            self._distinner = True
            self.distance = self.cmda.distance
            self.updatedata = True
            
    def _priband_changed(self,old,new):
        if new not in self.cmda.validbandnames:
            print 'band',new,'not valid'
            self.priband = old if new != old else ''
        else:
            self.updatepri = True
            
    def _pribandw_changed(self):
        if self.priband:
            self.updatepri = True
            
    def _agigpri_changed(self):
        self.updatepri = True
            
    def _pripower_changed(self):
        self.updatepri = True
        
    def _hideg0_changed(self):
        self.updatedata = True
        #m = self._get_g0m()
        #for v in  self.locplot.plots.itervalues():
        #    v = v[0]
        #    v.index.set_mask(m)
            
            
#    def _get_distance(self):
#        return self.cmda.distance
#    def _set_distance(self,val):
#        self.cmda.distance = val
#        self.updatedata = True
#        self.test = 1
        
#    def _get_distmod(self):
#        return self.cmda.distmod
#    def _set_distmod(self,val):
#        self.cmda.distmod = val
#        self.updatedata = True


    
    def _cmda_late_changed(self):
        self.distance = self.cmda.distance
        self.updatedata = True
        
    
    def _updatedata_fired(self):
        try:
            self._xax_late_changed()
            self._yax_late_changed()
        except ValueError:
            bi1,bi2 = self.cmda.validbandinds[0],cmda.validbandinds[1]
            self.xax = self.cmda.bandnames[bi1]
            self.yax = self.cmda.bandnames[bi2]
            self._xax_late_changed()
            self._yax_late_changed()
        self.updatepri = True
        self.updategastars = True
        
    def _xax_late_changed(self):
        self._update_x_data()
        self.cmplot.x_axis.title = self.xax
        
    def _get_g0m(self):
        if self.hideg0 and self.cmda.datagroups is not None:
            return self.cmda.datagroups !=0
        else:
            return np.ones(self.cmda._nd,dtype=bool)
        
            
    def _update_x_data(self):
        try:
            x = self.cmda.getData(self.xax)
            fx = self.cmda.getFiducialData(self.xax)
        except ValueError,e:
            if e.message == ' is not a valid fiducial band' or e.message == ' is not a valid band':
                return
        
        self._fix_axes()
            
        self.data.set_data('x',x)
        self.data.set_data('fidx',fx)
        self.updategastars = 'xd'
        
            
    def _yax_late_changed(self):
        self._update_y_data()
        self.cmplot.y_axis.title = self.yax
        
    def _update_y_data(self):
        try:
            y = self.cmda.getData(self.yax)
            fy = self.cmda.getFiducialData(self.yax)
        except ValueError,e:
            if e.message == ' is not a valid fiducial band' or e.message == ' is not a valid band':
                return
            
        self._fix_axes()
                
        self.data.set_data('y',y)
        self.data.set_data('fidy',fy)
        self.updategastars = 'yd'
        
    def _fix_axes(self):
        if self.cmplot is not Plot:
            oristr=[]
            if '-' not in self.yax:
                oristr.append('top')
            else:
                oristr.append('bottom')
            if '-' not in self.xax:
                oristr.append('right')
            else:
                oristr.append('left')
            oristr = ' '.join(oristr)
            self.cmplot.default_origin = oristr
            for pl in self.cmplot.plots.values():
                pl[0].origin = oristr
        
        
    def _updatelocs_fired(self):
        if self.cmda.locs is None:
            ra,dec = np.zeros(self.cmda._nd),np.zeros(self.cmda._nd)
        else:
            
            ra = self.cmda.locs[0]
            dec = self.cmda.locs[1]
            
        cen = self.cmda.center
        
        self.locs.set_data('ra',ra)
        self.locs.set_data('dec',dec)
        self.locs.set_data('centerx',(cen[0],))
        self.locs.set_data('centery',(cen[1],))
        if self.maskmaker is not None:
            l,w = self.maskmaker.size
            l2,w2=l/2.0,w/2.0
            self.locs.set_data('boxx',np.array([-w2,w2,w2,-w2,-w2])+cen[0])
            self.locs.set_data('boxy',np.array([-l2,-l2,l2,l2,-l2])+cen[1])
        else:
            self.locs.set_data('boxx',np.array(tuple()))
            self.locs.set_data('boxy',np.array(tuple()))
        
        self.updategastars = 'loc'
    
    def _late_dogastars(self):
        self.updategastars = True
    
    def _updategastars_fired(self,val):
        if self.agband:
            mg,ma = self._guidealignmasks()
            if val is True or val == 'loc':
                if self.cmda.locs is None:
                    ra,dec = np.zeros(self.cmda._nd),np.zeros(self.cmda._nd)
                else:
                    ra = self.cmda.locs[0]
                    dec = self.cmda.locs[1]
                self.locs.set_data('gra',ra[mg])
                self.locs.set_data('gdec',dec[mg])
                self.locs.set_data('ara',ra[ma])
                self.locs.set_data('adec',dec[ma])
                
            if val is True or val == 'xd':
                x = self.cmda.getData(self.xax)
                self.data.set_data('gx',x[mg])
                self.data.set_data('ax',x[ma])
                
            if val is True or val == 'yd':
                y = self.cmda.getData(self.yax)
                self.data.set_data('gy',y[mg])
                self.data.set_data('ay',y[ma])
        else:
            self.locs.set_data('gra',[])
            self.locs.set_data('gdec',[])
            self.locs.set_data('ara',[])
            self.locs.set_data('adec',[])
            self.locs.set_data('gx',[])
            self.locs.set_data('gy',[])
            self.locs.set_data('ax',[])
            self.locs.set_data('ay',[])
    def _guidealignmasks(self):
        mdat = self.cmda.getData(self.agband)
        pricut = self.agigpri > self.priorities
        mg = (mdat < self.gstarcut) & pricut
        ma = (mdat < self.astarcut) & pricut & ~mg
        return mg,ma
        
    def _offsetset_fired(self):
        od = OffsetDialog(self.cmda,self)
        od.edit_traits()
        
    def _masktype_changed(self,newval):
        if newval == 'DEIMOS':
            self.maskmaker = DEIMOSMaskMaker(self)
        else:
            raise TraitError('Unrecognized Mask Type')
        
        self.updatelocs = True
        
    def _maskshortname_changed(self,new,old):
        maxn = self.maskmaker.shortlen
        if len(new) > maxn:
            self.maskshortname = new[:maxn]
        
    def _makemaskbutton_fired(self):
        self.maskmaker.makeMask(self.maskfn)
        
        #increment number if there is one in the name
        name = self.maskshortname
        for cut in range(len(name)+1):
            try:
                num = int(name[cut:])
                num+=1
                self.maskshortname = name[:cut]+str(num)
                self.masklongname = self.maskshortname
                break
            except ValueError:
                pass
        #do the same for the file name
        name = self.maskfn
        for cut in range(len(name)+1):
            try:
                num = int(name[cut:])
                num+=1
                self.maskfn = name[:cut]+str(num)
                break
            except ValueError:
                pass
            
        
        
    def _get_priorities(self):
        
        off = self.cmda.getOffsets()
        pri = 1-off/off.max()
        
        if self.priband and self.pribandw:
            data = self.cmda.getData(self.priband)
            if self.pribandw < 0:
                data = (data-data.min())/(data.max()-data.min())
            else:
                data = (data.max()-data)/(data.max()-data.min())
            pri = pri+self.pribandw*data
            pri = (pri - pri.min())/(pri.max()-pri.min())
        
        pri = pri**self.pripower
        if self.hideg0:
#            if self.hideg0 and self.cmda.datagroups is not None:
#                m = self.cmda.datagroups ==1
#            else:
#                m =  slice(None)
#            pri[m] = 1
            pri[~self._get_g0m()] = 0
            
        return pri
        
    def _updatepri_fired(self):
        #pri = np.arange(len(self.data.get_data('x')))
        #pri = self.cmda.getOffsets()
        self.data.set_data('pri',self.priorities)
        self.locs.set_data('pri',self.priorities)
        if self.agigpri < 1:
            self.updategastars = True
        
class DEIMOSMaskMaker(MaskMaker):
    #import pyraf #do this to make sure pyraf is installed
    
    name = 'Dsim'
    shortlen = 6
    size = (5.0/60,16.0/60)
    
    exporttoiraf = False
    def makeMask(self,fn):
        from pyraf import iraf
        import os
        from os import path
        from warnings import warn
        
        
        fn = fn if fn.endswith('.in') else fn+'.in'
        basefn = path.splitext(fn)[0]
        print 'making',fn
        if path.exists(fn):
            warn('path %s exists, overwriting'%fn)
            os.remove(fn)
        self.writeInFile(fn)
        
        
        dskw={}
        dskw['output'] = basefn+'.dsim'
        dskw['mdf'] = basefn+'.fits'
        dskw['plotfile'] = basefn+'.plot'
        dskw['ra0'] = self.sto.cmda.center[0]*24.0/360.0
        dskw['dec0'] = self.sto.cmda.center[1]
        dskw['equinox'] = 2000
        dskw['guiname'] = self.sto.maskshortname
        dskw['maskid'] =  self.sto.masklongname.replace(' ','_') if self.sto.masklongname else self.sto.maskshortname
        print dskw['ra0'],dskw['dec0'],iraf.dsimulator.PA0
        for k in ('output','mdf','plotfile'):
            kfn = dskw[k]
            print 'deleteing',kfn
            if path.exists(kfn):
                os.remove(kfn)
        print 'running dsimulator'
        iraf.dsimulator(fn,**dskw)
        print 'dsimulator 1 complete!'
        fn2 = fn.replace('.in','.in2')
        dskw['ra0'],dskw['dec0'],dskw['PA0'] = self.secondFile(dskw['output'],fn2)
        print 'deleting intermediate files'
        for k in ('output','mdf','plotfile'):
            kfn = dskw[k]
            print 'deleteing',kfn
            if path.exists(kfn):
                os.remove(kfn)
        print 'Running with offset pa'
        iraf.dsimulator(fn2,**dskw)
        print 'setting selected objects to 0 group'
        self.doneToG0(dskw['output'],self.sto.cmda)
        print 'saving iraf parameters'
        with open(basefn+'.irafpars','w') as f:
            f.write('input:%s\n'%fn2)
            for t in dskw.iteritems():
                f.write('%s\t=\t%s\n'%t)
        if self.exporttoiraf:
            print 'exporting settings to iraf'
            for k,v in dskw.iteritems():
                setattr(iraf.dsimulator,k,v)
                #iraf.dsimulator.output = ''
                #iraf.dsimulator.mdf = ''
                #iraf.dsimulator.plotfile = ''
                iraf.dsimulator.saveParList()
        print 'all done!'
        
    def writeInFile(self,fn,pa=None):
        from astropysics.coords import AngularCoordinate
        
        pri = 1000*self.sto.priorities
        ni = len(pri)
        mg,ma = self.sto._guidealignmasks()
        
        
        ra,dec = self.sto.cmda.locs
        magband = self.sto.priband if self.sto.priband else self.sto.agband
        
        mags = self.sto.cmda.getData(magband)
        self.sto.groupcutoff
        if self.sto.cmda.datagroups is None:
            samp = np.ones(ni)
            samp[mags > self.sto.groupcutoff] = 2
            pri[mg] = -1
            pri[ma] = -2
        else:
            samp = self.sto.cmda.datagroups.copy()
            ssamp = set(samp)
            keepbright = samp==1 if self.sto.agigg1 else np.zeros(samp.shape,dtype=bool)
            pri[(samp==-1) & ~ma] = -1
            pri[(samp==-2) & ~mg] = -2
            #pri[samp==0] = 0 #do this automatically or not from the GUI 
            pri[mg & ~keepbright] = -1
            pri[ma & ~keepbright] = -2
            samp[(samp==-1) | (pri==-1)] = 6
            samp[(samp==-2) | (pri==-2)] = 6
            samp[(mags > self.sto.groupcutoff) & (samp>0)] = 5
            #samp[pri==0] = 5
        if pa is not None:
            pa = np.ones(ni)*pa
        
        #pri[pri==0]=1 #avoid priority-0 cases due to dsim weirdness
    
        presel = np.zeros(ni)
        with open(fn,'w') as f:
            for i,n in enumerate(self.sto.cmda.datanames):
                rao,deco = AngularCoordinate(ra[i]),AngularCoordinate(dec[i])
                ras,decs = rao.getHmsStr(canonical = True),deco.getDmsStr(canonical = True)
                if samp[i] != -9 or pri[i]==-1 or pri[i] == -2: #-9 is a code to not write at all
                    if pa is None:
                        t = (n[:16],ras,decs,2000,mags[i],magband,int(pri[i]),samp[i],presel[i])
                        f.write('%s\t%s %s %i %2.2f %s %4i %i %i INDEF\n'%t)
                    else:
                        t = (n[:16],ras,decs,2000.0,mags[i],magband,int(pri[i]),samp[i],presel[i],pa[i])
                        f.write('%s\t%s\t%s\t%06.1f\t%2.2f\t%s\t%i\t%i\t%i\t%i\n'%t)
        
#    def fixPaFile(self,fn1,fn2,degoffset=5):
#        #fn1 is th dsim output of the first dsim run, fn2 is the target file name
#        from astropysics.coords import AngularPosition
        
#        with open(fn1,'r') as fr:
#            with open(fn2,'w') as fw:
#                for l in fr:
#                    if 'PA' in l:
#                        ras,decs,epochs,pas = l.split()[-5:-1]
#                        panum = float(pas.replace('PA=',''))
#                        fw.write(l)
#                    else:
#                        ls = l.strip()
#                        if ls.startswith('#') or ls == '':
#                            fw.write(l)
#                        else:
#                            fw.write(ls)
#                            fw.write('  %i\n'%(panum+degoffset))
#        ap = AngularPosition(ras,decs)
#        return ap.ra.d,ap.dec.d,panum
    
    def secondFile(self,fn1,fn2,degoffset=5):
        #fn1 is th dsim output of the first dsim run, fn2 is the target file name
        from astropysics.coords import AngularPosition,AngularCoordinate
        
        with open(fn1,'r') as fr:
            with open(fn2,'w') as fw:
                for l in fr:
                    if 'PA' in l:
                        ls = l.split()
                        for i in range(len(ls)-1):
                            try:
                                raac = AngularCoordinate(ls[i],sghms=True)
                                decac = AngularCoordinate(ls[i+1],sghms=False)
                                break
                            except ValueError:
                                pass
                        pas = ls[-2]
                        panum = float(pas.replace('PA=',''))
                        fw.write(l)
                    else:
                        ls = l.strip()
                        if ls.startswith('#') or ls == '':
                            fw.write(l)
                        else:
                            #fw.write(ls)
                            #fw.write('  %i\n'%(panum+degoffset))
                            name,ra,dec,epoch,mag,band,pri,grp,presel = ls.split()[:9]
                            if int(presel):
                                pri = int(pri)
                            else:
                                pri = -2 if (float(mag) < self.sto.astarcut2) else int(pri)
                            t=(name,ra,dec,epoch,mag,band,pri,grp,presel,panum+degoffset)
                            fw.write('%s\t%s  %s  %s %s %s %4.1i %s %s %i\n'%t)
                            
        ap = AngularPosition(raac,decac)
        return ap.ra.d,ap.dec.d,panum
    
    @staticmethod
    def doneToG0(fn,cmda,forceradec=False):
        objs=[]
        ras,decs=[],[]
        
        with open(fn) as f:
            watch = False
            for l in f:
                if watch and l.strip() and not l.lstrip().startswith('#'):
                    ls=l.split()
                    if int(ls[6]) >= 0 and int(ls[8]) >= 0:
                        objs.append(int(ls[0]))
                        ras.append(ls[1])
                        decs.append(ls[2])
                if 'Selected Objects' in l:
                    watch = True
                if 'Selected Guide' in l:
                    #watch = False
                    break
        objs = np.array(objs)-1
        
        if cmda.datagroups is not None:
            g = cmda.datagroups
        else:
            g = np.ones(cmda._nd)
            
        
        if forceradec or objs.max() >= len(g):
            #if necessary, do ra/dec matching instead of direct indexing
            from ..coords import match_coords,radec_str_to_decimal
            rao,deco = cmda.locs
            ras,decs = radec_str_to_decimal(ras,decs)
            mask = match_coords(rao,deco,ras,decs,0.5/3600)[0]
            g[mask] = 0
        else:
            if len(objs) > 0:
                g[objs] = 0
        if not all(g):
            cmda.datagroups = g
                
def spec_target(*args):
    """
    generates and shows a SpecTarget gui interface - accepts two forms:
    
    * spec_target(cmdanalyzer) - provide a ready-made CMDAnalyzer
    * spec_target(fiducialdict,data) - provide a dictionary of fiducial 
        data and a set of data (array or dict)
        
    A GUI application instance must already exist (e.g. interactive mode of 
    ipython)    
        
    returns the SpecTarget instance 
    """
    st = SpecTarget(*args)
    st.edit_traits()
    return st
    