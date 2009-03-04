#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the SpecTarget gui.
"""
from __future__ import division
import numpy as np

from enthought.traits.api import HasTraits,on_trait_change,Instance,Float,\
                                 String,ListStr,ListFloat,Button,Bool,\
                                 Property,Event,Array,Enum,TraitError
from enthought.traits.ui.api import View,Group,HGroup,HFlow,Item,Handler
from enthought.traits.ui.menu import ModalButtons
from enthought.chaco.api import Plot,ArrayPlotData,HPlotContainer,ColorBar,LinearMapper,jet
from enthought.chaco.tools.api import PanTool, ZoomTool
from enthought.enable.component_editor import ComponentEditor

from ..phot import CMDAnalyzer

class MaskMaker(object):
    """
    base class for objects that make masks - subclasses should:
    *replace name with a name for the mask-making operation
    *override makeMask(self,filename,spectargobj)
    """
    name = 'Default mask-maker'
    def makeMask(self,filename,spectargobj):
        raise NotImplementedError

class ODHandler(Handler):
    def apply(self,info):
        o = info.object
        o.cmda.offsetbands = o.offsetbands
        o.cmda.offsetweights = o.offsetweights
        o.st.updatepri = True
        
    def closed(self,info, is_ok):
        if is_ok:
            self.apply(info)
        Handler.closed(self,info,is_ok)

class OffsetDialog(HasTraits):
    offsetbands = ListStr
    offsetweights = Array(float)
    
    view = View(Group(Item('offsetbands'),
                      Item('offsetweights'),
                      layout='normal'),
                resizable=True,title='Offset Settings',handler=ODHandler(),
                buttons = ModalButtons,kind='modal')
    
    def __init__(self,cmda,st):
        self.cmda = cmda
        self.st = st
        self.offsetbands = cmda.offsetbands if cmda.offsetbands is not None else []
        self.offsetweights = cmda.offsetweights
            
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
    agband = String
    astarcut = Float(18)
    gstarcut = Float(18)
    offsetset = Button
    priorities = Property(depends_on='distmod,offsetset')
    
    masktype = Enum(['DEIMOS'])
    maskmaker = Instance(MaskMaker) 
    maskfn = String('spectarg-1')
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
                                  columns=3),
                            Item('offsetset',label='Offset Settings'),
                            Group(Item('masktype',label='Mask Type'),
                                  Item('maskfn',label='Mask File(base)'),
                                  Item('makemaskbutton'),#,label='Make Mask'),
                                  columns=3)),
                      layout='normal'),
                resizable=True,title='Spectra Targeter')
    
    
    def __init__(self,*args,**kwargs):
        """
        generates the SpecTarget GUI - args can be either:
        
        *SpecTarget(cmdanalyzer) - provide a ready-made CMDAnalyzer
        *SpecTarget(fiducialdict,data) - provide a dictionary of fiducial 
        data and a set of data (array or dict)
        
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
        
        if len(ob) > 1:
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
        cmdplot.plot(('x','y','pri'),name='target_mags',type='cmap_scatter',marker='dot',marker_size=2,color_mapper=jet)
        cmdplot.plot(('fidx','fidy'),name='fid_mags',type='line',color='black')
        cmdplot.plot(('gx','gy'),name='guide_mags',type='scatter',color='gray',marker='inverted_triangle',outline_color='red',marker_size=3)  
        cmdplot.plot(('ax','ay'),name='align_mags',type='scatter',color='black',marker='diamond',outline_color='red',marker_size=3)
        cmdplot.tools.append(PanTool(cmdplot))
        cmdplot.tools.append(ZoomTool(cmdplot))  
        
        locplot = Plot(self.locs,resizable='hv')
        locplot.plot(('ra','dec','pri'),name='target_locs',type='cmap_scatter',marker_size=3,color_mapper=jet)
        locplot.plot(('centerx','centery'),name='fid_locs',type='scatter',color='black',marker='cross',marker_size=10,line_width=3)
        locplot.plot(('gra','gdec'),name='guide_locs',type='scatter',color='gray',marker='inverted_triangle',outline_color='red',marker_size=6)
        locplot.plot(('ara','adec'),name='align_locs',type='scatter',color='black',marker='diamond',outline_color='red',marker_size=6)
        locplot.tools.append(PanTool(locplot))
        locplot.tools.append(ZoomTool(locplot))  
        
        cb = ColorBar(index_mapper=LinearMapper(range=cmdplot.color_mapper.range),color_mapper=cmdplot.color_mapper,resizable='v',width=30,padding=5)
        
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
        self.on_trait_change(self._late_dogastars,'agband,gstarcut,astarcut')
        
        self._cmda_late_changed() #clean things up in case random stuff changes
        self._masktype_changed(self.masktype) #call to intialize default mask maker
        
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
        self.dist = self.cmda.distance
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
        
        self.updategastars = 'loc'
    
    def _late_dogastars(self):
        self.updategastars = True
    
    def _updategastars_fired(self,val):
        if self.agband:
            mdat = self.cmda.getData(self.agband)
            mg = mdat < self.gstarcut
            ma = (mdat < self.astarcut) & ~mg
            
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
        
    def _offsetset_fired(self):
        od = OffsetDialog(self.cmda,self)
        od.edit_traits()
        
    def _masktype_changed(self,newval):
        if newval == 'DEIMOS':
            self.maskmaker = DEIMOSMaskMaker()
        else:
            raise TraitError('Unrecognized Mask Type')
        
    def _makemaskbutton_fired(self):
        self.maskmaker.makeMask(self.maskfn,self)
        
        
    def _get_priorities(self):
        off = self.cmda.getOffsets()
        off = 1-off/off.max()
        return off**4
        
    def _updatepri_fired(self):
        #pri = np.arange(len(self.data.get_data('x')))
        #pri = self.cmda.getOffsets()
        self.data.set_data('pri',self.priorities)
        self.locs.set_data('pri',self.priorities)
        
class DEIMOSMaskMaker(MaskMaker):
    import pyraf #do this to make sure pyraf is installed
    
    name = 'Dsim'
    def makeMask(self,fn,st):
        from pyraf import iraf
        from os import path
        
        
        fn = fn if fn.endswith('.in') else fn+'.in'
        print 'making',fn
        
        
    def writeInFile(self,fn)
    