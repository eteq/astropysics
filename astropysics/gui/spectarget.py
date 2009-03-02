#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the SpecTarget gui.
"""
from __future__ import division
import numpy as np

from enthought.traits.api import HasTraits,on_trait_change,Instance,Float,String,ListStr
from enthought.traits.ui.api import View,Group,HGroup,HFlow,Item
from enthought.chaco.api import Plot,ArrayPlotData
from enthought.enable.component_editor import ComponentEditor

from ..phot import CMDAnalyzer

class SpecTarget(HasTraits):
    distance = Float(10)
    distmod = Float(0)
    cmda = Instance(CMDAnalyzer)
    locplot = Plot
    cmplot = Plot
    data = ArrayPlotData
    locs = ArrayPlotData
    xax = String
    yax = String
    offsetbands = ListStr
    
    view = View(Group(HGroup(Item('locplot',editor=ComponentEditor(),show_label=False),
                             Item('cmplot',editor=ComponentEditor(),show_label=False)),
                      Group(Item('distance',label='Distance (kpc)'),
                            Item('distmod',label='Distance Modulus'),
                            Item('xax',label='CMD x-axis'),
                            Item('yax',label='CMD y-axis'),
                            Item('offsetbands',label='Bands for Calculating Offsets'),
                            layout='normal'),
                      ),
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
        
        offsetbands = cmda.offsetbands if cmda.offsetbands is not None else []
        
        if len(offsetbands) > 1:
            self.xax,self.yax = offsetbands[0],offsetbands[1]
        else:
            bi1,bi2 = cmda.validbandinds[0],cmda.validbandinds[1]
            self.xax,self.yax = cmda.bandnames[bi1],cmda.bandnames[bi2]

        self.data = ArrayPlotData()
        self._update_x_data()
        self._update_y_data()
                                  
        cmdplot = Plot(self.data,resizable='v')
        self.cmdrender1 = cmdplot.plot(('x','y'),type='scatter',color='blue',marker='dot')
        self.cmdrender2 = cmdplot.plot(('fidx','fidy'),type='line',color='red')   
        
            
        self.locs = ArrayPlotData()
        self.update_locs()
        
        locplot = Plot(self.locs,resizable='v')
        self.locrender1 = locplot.plot(('ra','dec'),type='scatter',color='blue')
        self.locrender2 = locplot.plot(('centerx','centery'),type='scatter',color='black',marker='cross',marker_size=10,line_width=3)
        
        #self.data = ArrayPlotData()
        self.locplot = locplot
        self.cmplot = cmdplot
        
        self.on_trait_change(self._cmda_late_changed,'cmda')
        self.on_trait_change(self._xax_late_changed,'xax')
        self.on_trait_change(self._yax_late_changed,'yax')
        
        self._cmda_late_changed() #clean things up in case random stuff changes
    
    def _dist_changed(self):
        self.cmda.distance = dist
        self.distmod = self.cmda.distmod
        
    def _distmod_changed(self):
        self.cmda.distmod = distmod
        self.distance = self.cmda.distance
        
    def _offsetbands_changed(self):
        self.cmda.offsetbands = self.offsetbands
    
    def _cmda_late_changed(self):
        self.dist = self.cmda.distance
        self._update_data()
        self.offsetbands = self.cmda.offsetbands if self.cmda.offsetbands is not None else []
    
    def _update_data(self):
        try:
            self._xax_late_changed()
            self._yax_late_changed()
        except ValueError:
            bi1,bi2 = self.cmda.validbandinds[0],cmda.validbandinds[1]
            self.xax = self.cmda.bandnames[bi1]
            self.yax = self.cmda.bandnames[bi2]
            self._xax_late_changed()
            self._yax_late_changed()
        
    def _xax_late_changed(self):
        self._update_x_data()
        self.cmplot.x_axis.title = ('-' if '-' not in self.xax else '')+self.xax
       
        
    def _update_x_data(self):
        try:
            x = self.cmda.getData(self.xax)
            fx = self.cmda.getFiducialData(self.xax)
        except ValueError,e:
            if e.message == ' is not a valid fiducial band' or e.message == ' is not a valid band':
                return
        #TODO: figure out how to ivert the axes:
        if '-' not in self.xax:
            x , fx = -1*x,fx*-1
            
        self.data.set_data('x',x)
        self.data.set_data('fidx',fx)
        
            
    def _yax_late_changed(self):
        self._update_y_data()
        self.cmplot.y_axis.title = ('-' if '-' not in self.yax else '')+self.yax
        
    def _update_y_data(self):
        try:
            y = self.cmda.getData(self.yax)
            fy = self.cmda.getFiducialData(self.yax)
        except ValueError,e:
            if e.message == ' is not a valid fiducial band' or e.message == ' is not a valid band':
                return
        #TODO: figure out how to ivert the axes:
        if '-' not in self.yax:
            y , fy = -1*y,fy*-1
            
        self.data.set_data('y',y)
        self.data.set_data('fidy',fy)
        
        
    def update_locs(self):
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
        
    
def do_dsim():
    from pyraf import iraf
    iraf.dsimulator(infile,**kwargs)