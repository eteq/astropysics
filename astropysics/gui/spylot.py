#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for spylot, the Python Spectrum
Plotter.
"""
from __future__ import division,with_statement
import numpy as np

from enthought.traits.api import HasTraits,List,Instance,Property,Int,Range, \
                                 Float,Event,Bool,DelegatesTo,Dict,Str,List, \
                                 Button,cached_property,on_trait_change
from enthought.traits.ui.api import View,Item,Label,RangeEditor,Group,HGroup, \
                                    VGroup,HFlow,SetEditor,spring

from enthought.chaco.api import Plot,ArrayPlotData,PlotAxis,LinearMapper, \
                                DataRange1D,DataLabel,LinePlot,ArrayDataSource
from enthought.chaco.tools.api import PanTool, ZoomTool,BroadcasterTool
from enthought.enable.api import ComponentEditor

from .. import spec

def _hist_sample_x(x):
    if len(x)==1:
        return x
    x = np.convolve(x,[.5,.5],mode='valid')
    l = (2*x[0]-x[1],)
    r = (2*x[-1]-x[-2],)
    return np.concatenate((l,np.repeat(x,2),r))

def _hist_sample_y(y):
    if len(y)==1:
        return y
    else:
        return np.repeat(y,2)
    
    
class LineListEditor(HasTraits):
    candidates = List
    candidatenames = Property(depends_on='candidates')
    selectednames = List
    selectedobjs = Property(depends_on='selectednames')
    selectedlocs = Property(depends_on='selectednames')
    lineplot = Instance(LinePlot)
    
    major_view = View(Item('selectednames',editor=SetEditor(name='candidatenames',
                                    left_column_title='Candidate Lines',
                                    right_column_title='Selected Lines'),
                                    show_label=False),
                        Item('lineplot',show_label=False,style='custom'),
                        title='Major Lines')
                                        
    minor_view = View(Item('selectednames',editor=SetEditor(name='candidatenames',
                                    left_column_title='Candidate Lines',
                                    right_column_title='Selected Lines'),
                                    show_label=False),
                        Item('lineplot',show_label=False,style='custom'),
                        title='Minor Lines')

    @cached_property
    def _get_candidatenames(self):
        return tuple([line.name for line in self.candidates])
    
    @cached_property
    def _get_selectedlocs(self):
        return tuple([o.loc for o in self.selectedobjs])
    
    @cached_property
    def _get_selectedobjs(self):
        return tuple([self.candidates[self.candidatenames.index(n)] for n in self.selectednames])
    
    
    def _set_selectedobjs(self,val):
        cs,cns,cls = [],[],[]
        for c in self.candidates:
            cs.append(c)
            cns.append(c.name)
            cls.append(c.loc)
        
        names = []
        for v in val:
            if isinstance(v,basestring):
                if v not in cns:
                    raise ValueError('line name '+v+' not found')
                names.append(v)
            elif np.isscalar(v):
                if v not in cls:
                    raise ValueError('line location '+str(v)+' not found')
                names.append(cns[cls.index(v)])
            else:
                if v not in cs:
                    raise ValueError('line object '+str(v)+' not found')
                names.append(cns[cs.index(v)])
        self.names = names

class Spylot(HasTraits):
    defaultlines = 'galaxy' #can be 'galaxy' or 'stellar' or None
    
    plot = Instance(Plot)
    histplot = Bool(True)
    
    majorlineeditor = LineListEditor
    minorlineeditor = LineListEditor
    showmajorlines = Bool(True)
    showminorlines = Bool(False)
    
    #majorlabels = List(DataLabel)
    #showmajorlabels = Bool(False)
    editmajor = Button('Edit Major')
    #minorlabels = List(DataLabel)
    #showminorlabels = Bool(False)
    editminor = Button('Edit Minor')
    
    specs = List(Instance(spec.Spectrum))
    currspeci = Int
    currspec = Property
    z = Float
    lowerZ = Float(0.0)
    upperZ = Float(1.0)
    _zql,_zqh = min(spec.Spectrum._zqmap.keys()),max(spec.Spectrum._zqmap.keys())
    zqual = Range(_zql,_zqh,-1)
    
    spechanged = Event
    specleft = Button('<')
    specright = Button('>')
    
    scaleerr = Bool(False)
    scaleerrfrac = Float(0.7)
    fluxformat = Button('Flux Line Format')
    errformat = Button('Error Line Format')
    showgrid = Bool(True)
    showzero = Bool(False)
    
    _titlestr = Str('Spectrum 0/0')
    _oldunit = Str('')
    
    traits_view = View(VGroup(HGroup(Item('specleft',show_label=False,enabled_when='currspeci>0'),spring,
                                     Item('_titlestr',style='readonly',show_label=False),spring,
                                     Item('specright',show_label=False,enabled_when='currspeci<(len(specs)-1)')),
                               HGroup(spring,Item('showmajorlines',label='Show major lines?'),
                                     Item('editmajor',show_label=False),
                                     Item('showmajorlines',label='Show major lines?'),
                                     Item('editmajor',show_label=False),spring),     
#                              HGroup(spring,Item('showmajorlines',label='Show major lines?'),
#                                     Item('showmajorlabels',label='Show major labels?'),
#                                     Item('editmajor',show_label=False),spring),
#                              HGroup(spring,Item('showminorlines',label='Show minor lines?'),
#                                     Item('showminorlabels',label='Show minor labels?'),
#                                     Item('editminor',show_label=False),spring),
                              HGroup(spring,
                                     Item('scaleerr',label='Scale Error to flux?'),
                                     Item('scaleerrfrac',show_label=False,enabled_when='scaleerr'),
                                     Item('fluxformat',show_label=False),
                                     Item('errformat',show_label=False),
                                     Item('showgrid',label='Grid?'),
                                     Item('showzero',label='Zero line?'),spring),
                              Item('plot',editor=ComponentEditor(),show_label=False,width=768),
                              Item('z',editor=RangeEditor(low_name='lowerZ',high_name='upperZ',format='%.3f',auto_set=True)),
                              HGroup(Item('lowerZ'),spring,
                                     Item('zqual',style='custom',label='Z quality',editor=RangeEditor(cols=_zqh-_zql+1,low=_zql,high=_zqh)),
                                     spring,Item('upperZ'))
                      ),
                       resizable=True, 
                       title='Spectrum Plotter',
                       )
                      
    def __init__(self,specs):
        #pd = ArrayPlotData(x=[1],x0=[1],flux=[1],err=[1]) #reset by spechanged event
        pd = ArrayPlotData(x=[1],flux=[1],err=[1]) #reset by spechanged event
        pd.set_data('majorx',[1,1])#reset by majorlines change
        pd.set_data('majory',[0,1])#reset by majorlines change
        pd.set_data('minorx',[1,1])#reset by minorlines change
        pd.set_data('minory',[0,1])#reset by minorlines change
        
        self.plot = plot = Plot(pd,resizeable='hv')
        plot.plot(('x','flux'),name='flux',type='line',line_style='solid',color='blue')
        plot.plot(('x','err'),name='err',type='line',line_style='dash',color='green')
        
        topmapper = LinearMapper(range=DataRange1D())
        plot.x_mapper.range.on_trait_change(self._update_upperaxis_range,'updated')
        plot.x_mapper.on_trait_change(self._update_upperaxis_screen,'updated')
        self.upperaxis = PlotAxis(plot,orientation='top',mapper=topmapper)
        plot.overlays.append(self.upperaxis)
        
        plot.padding_top = 30 #default is a bit much
        
        majorlineplot = plot.plot(('majorx','majory'),name='majorlineplot',type='line',line_style='dash',color='red')[0]
        #majorlineplot.value_mapper = LinearMapper(range=DataRange1D(plot._get_or_create_datasource('majory')))
        majorlineplot.value_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
        majorlineplot.visible = self.showmajorlines
        del plot.x_mapper.range.sources[-1]  #remove the line plot from the x_mapper sources so scaling is only on the spectrum
        self.majorlineeditor = LineListEditor(lineplot=majorlineplot)
        
        minorlineplot = plot.plot(('minorx','minory'),name='minorlineplot',type='line',line_style='dot',color='red')[0]
        #minorlineplot.value_mapper = LinearMapper(range=DataRange1D(plot._get_or_create_datasource('minory')))
        minorlineplot.value_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
        minorlineplot.visible = self.showminorlines
        del plot.x_mapper.range.sources[-1]  #remove the line plot from the x_mapper sources so scaling is only on the spectrum
        self.minorlineeditor = LineListEditor(lineplot=minorlineplot)
        
        idat = ArrayDataSource((0.0,1.0))
        vdat = ArrayDataSource((0.0,0.0))
        self.zeroline = LinePlot(index=idat,value=vdat,line_style='solid',color='black')
        self.zeroline.index_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
        self.zeroline.value_mapper = self.plot.y_mapper
        self.zeroline.visible = self.showzero
        self.plot.add(self.zeroline)
        
        if Spylot.defaultlines:
            defaultlines = _get_default_lines(Spylot.defaultlines)
            self.majorlineeditor.candidates = defaultlines[0]
            self.majorlineeditor.selectednames = defaultlines[1]
            self.minorlineeditor.candidates = defaultlines[0]
            self.minorlineeditor.selectednames = defaultlines[2]
        
        plot.tools.append(PanTool(plot))
        plot.overlays.append(ZoomTool(plot))
        
        if specs is None:
            specs = []
        elif isinstance(specs,spec.Spectrum):
            specs = [specs]
        self.specs = specs
        
        self.spechanged = True
        
        self.on_trait_change(self._majorlines_changed,'majorlineeditor.selectedobjs')
        self.on_trait_change(self._minorlines_changed,'minorlineeditor.selectedobjs')
        
    def _update_upperaxis_range(self):
        newlow = self.plot.x_mapper.range.low/(self.z+1)
        newhigh = self.plot.x_mapper.range.high/(self.z+1)
        self.upperaxis.mapper.range.set_bounds(newlow,newhigh)
        
    def _update_upperaxis_screen(self):
        self.upperaxis.mapper.screen_bounds = self.plot.x_mapper.screen_bounds
        
    def _get_currspec(self):
        return self.specs[self.currspeci]
    def _set_currspec(self,val):
        self.currspeci = self.specs.index(val)
    
    def _currspeci_changed(self,new,old):
        nsp = len(self.specs)
        if nsp == 0:
            self.currspeci = 0
            if new != 0:
                raise IndexError("No spectra present - can't choose current")
        elif not (0 <= new < nsp):
            self.currspeci = old
            raise IndexError('Spectrum index %i invalid'%new)
        self.spechanged = True
    
    def _specs_changed(self):
        #TODO: check more carefully to make sure changes don't invalidate currspeci
        self.spechanged = True
        
    def _histplot_changed(self):
        self.spechanged = True
    
    def _z_changed(self):
        self.currspec.z = self.z
        #self.plot.data.set_data('x0',_hist_sample_x(self.currspec.x0) if self.histplot else self.currspec.x0)
        self._majorlines_changed()
        self._minorlines_changed()
        self._update_upperaxis_range()
        
    def _zqual_changed(self):
        self.currspec.zqual = self.zqual
        
    def _scaleerr_changed(self):
        self.spechanged = True #TODO: control more finely
        
    def _scaleerrfrac_changed(self):
        self.spechanged = True #TODO: control more finely
        
    def _spechanged_fired(self):
        p = self.plot
        s = self.currspec
        
        if s.unit != self._oldunit:
            self._minorlines_changed()
            self._majorlines_changed()
            self._oldunit = s.unit
        
        if self.histplot:
            x = _hist_sample_x(s.x)
            #x0 = _hist_sample_x(s.x0)
            flux = _hist_sample_y(s.flux)
            err = _hist_sample_y(s.err)
        else:
            #x0 = s.x0
            x = s.x
            flux = s.flux
            err = s.err
        
        if self.scaleerr:
            err = err*self.scaleerrfrac*flux.max()/err.max()
            
        p.data.set_data('x',x)
        #p.data.set_data('x0',x0)
        p.data.set_data('flux',flux)
        p.data.set_data('err',err)
        
        self.z  = s.z
        self.zqual = s.zqual
        
        units = tuple(s.unit.split('-'))
        p.x_axis.title = '%s/%s'%units
        self.upperaxis.title = 'rest %s/%s'%units
        p.y_axis.title ='Flux/(erg s^-1 cm^-2 %s^-1)'%units[1]
        p.request_redraw()
        
        self._titlestr = 'Spectrum %i/%i: %s'%(self.currspeci+1,len(self.specs),s.name)
        
    def _specleft_fired(self):
        self.currspeci -= 1
    
    def _specright_fired(self):
        self.currspeci += 1

    def _majorlines_changed(self):
        majorlines = self.majorlineeditor.selectedobjs
        if len(majorlines)>0:
            u = self.currspec.unit
            linelocs = [l.getUnitLoc(u) for l in majorlines]
            linex = np.repeat(linelocs,2)
            linex.sort()
            liney = np.tile((0,1,1,0),np.ceil(len(majorlines)/2))[:len(majorlines)*2]
            self.plot.data.set_data('majorx',linex*(self.z+1))
            self.plot.data.set_data('majory',liney)
            self.plot.plots['majorlineplot'][0].request_redraw()
        else:
            self.plot.data.set_data('majorx',[0])
            self.plot.data.set_data('majory',[0])
    def _showmajorlines_changed(self):
        majorlineplot = self.plot.plots['majorlineplot'][0]
        majorlineplot.visible = self.showmajorlines 
        majorlineplot.request_redraw()
        
    def _minorlines_changed(self):
        minorlines = self.minorlineeditor.selectedobjs
        if len(minorlines)>0:
            u = self.currspec.unit
            linelocs = [l.getUnitLoc(u) for l in minorlines]
            linex = np.repeat(linelocs,2)
            linex.sort()
            liney = np.tile((0,1,1,0),np.ceil(len(minorlines)/2))[:len(minorlines)*2]
            self.plot.data.set_data('minorx',linex*(self.z+1))
            self.plot.data.set_data('minory',liney)
            self.plot.plots['minorlineplot'][0].request_redraw()
        else:
            self.plot.data.set_data('minorx',[0])
            self.plot.data.set_data('minory',[0])
    
    def _showminorlines_changed(self):
        minorlineplot = self.plot.plots['minorlineplot'][0]
        minorlineplot.visible = self.showminorlines 
        minorlineplot.request_redraw()
        
    def _editmajor_fired(self):
        self.majorlineeditor.edit_traits(view='major_view')
    
    def _editminor_fired(self):
        self.minorlineeditor.edit_traits(view='minor_view')
        
    def _fluxformat_fired(self):
        from copy import copy
        plot = self.plot.plots['flux'][0]
        v = copy(plot.trait_view())
        v.title = 'Flux Line Format'
        plot.edit_traits(view=v)
    
    def _errformat_fired(self):
        from copy import copy
        plot = self.plot.plots['err'][0]
        v = copy(plot.trait_view())
        v.title = 'Error Line Format'
        plot.edit_traits(view=v)
        
    def _showgrid_changed(self,new):
        p = self.plot
        p.x_grid.visible = new
        p.y_grid.visible = new
        p.request_redraw()
        
    def _showzero_changed(self,new):
        self.zeroline.visible = new
        self.plot.request_redraw()
            
def _get_default_lines(linetypes):
    candidates = spec.load_line_list(linetypes)
    if linetypes == 'galaxy':
        major=['Ly_alpha','H_alpha','H_beta','H_gamma','[OIII]_5007','[OIII]_4959','[OII]_3727','Balmer_limit','[SII]_6717','[SII]_6731']
        return candidates,major,[c.name for c in candidates if c.name not in major]
    elif linetypes == 'stellar':
        major = ['Halpha']
        return candidates,major,[c.name for c in candidates if c not in major]
    else:
        return candidates,[],[]
    
def spylotSpecs(specs):
    sp = Spylot(specs)
    sp.configure_traits()
    return sp