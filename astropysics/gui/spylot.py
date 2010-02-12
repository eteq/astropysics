#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This module contains the internals for spylot, the Python Spectrum
Plotter.
"""
#TODO: clean up upper interface
#TODO: smarter scaling/edge cleaning, inc scaling history ("reset" option?)
#TODO: objcat support
#TODO: remember line lists (better placed in spec)
#TODO: add diagnostics like SNR

from __future__ import division,with_statement
import numpy as np
import os

from enthought.traits.api import HasTraits,List,Instance,Property,Int,Range, \
                                 Float,Event,Bool,DelegatesTo,Dict,Str,List, \
                                 Button,cached_property,on_trait_change,Enum, \
                                 Tuple,Array,Trait,File
from enthought.traits.ui.api import View,Item,Label,RangeEditor,Group,HGroup, \
                                    VGroup,HFlow,SetEditor,spring,Handler, \
                                    TextEditor,ColorTrait,TabularEditor, \
                                    Include
from enthought.traits.ui.key_bindings import KeyBinding,KeyBindings
from enthought.traits.ui.tabular_adapter import TabularAdapter
from enthought.chaco.api import Plot,ArrayPlotData,PlotAxis,LinePlot, \
                                DataRange1D,DataLabel,ArrayDataSource, \
                                AbstractOverlay,LinearMapper,TextBoxOverlay
from enthought.chaco.tools.api import PanTool, ZoomTool,BroadcasterTool, \
                                      RangeSelection,RangeSelectionOverlay, \
                                      DragTool    
from enthought.enable.api import ComponentEditor,BaseTool,Interactor,Line, \
                                 Component 

from .. import spec

def _hist_sample_x(x):
    if len(x)<3:
        return x
    x = np.convolve(x,[.5,.5],mode='valid')
    l = (2*x[0]-x[1],)
    r = (2*x[-1]-x[-2],)
    return np.concatenate((l,np.repeat(x,2),r))

def _hist_sample_y(y):
    if len(y)<3:
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
        return tuple([str(line) for line in self.candidates])
    
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
            cns.append(str(c))
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


#overlays and tools
class FeatureClickTool(BaseTool):
    
    event_state = Enum("normal", "mousedown")
    plot = Instance(Plot)
    mouse_clicked = Event

    def normal_left_down(self, event):
        self.event_state = 'mousedown'

    def mousedown_left_up(self, event):
        datax,datay = self.plot.map_data((event.x, event.y))
        self.mouse_clicked = (datax,datay)
        
class MouseMoveReporter(BaseTool):
    overlay = Instance(TextBoxOverlay)
    plot = Instance(Plot)
    def normal_mouse_move(self, event):
        if self.plot is not None:
            datax,datay = self.plot.map_data((event.x, event.y))
        else:
            x,y = event.x,event.y
        if self.overlay is not None and self.overlay.visible:
            self.plot.request_redraw()
            self.overlay.text = '%.3f,%.3f'%(datax,datay)

class RangeSelectionBugfix(RangeSelection):
    """
    this subclasses RangeSelection to fix a bug where the selection_completed
    event is not fired for the left mouse button
    """
    def selecting_left_up(self, event):
        self.event_state = "selected"
        self.selection_completed = self.selection
        return        
    def normal_left_down(self,event):
        self.downwindow = event.window
        return RangeSelection.normal_left_down(self,event)
    
    def normal_right_down(self,event):
        self.downwindow = event.window
        if self.left_button_selects and event.right_down:
            return
        else:
            return RangeSelection.normal_right_down(self,event)
        
    def selecting_mouse_move(self, event):
        """ Handles the mouse being moved when the tool is in the 'selecting'
        state.
        
        Expands the selection range at the appropriate end, based on the new
        mouse position.
        """
        if self.selection is not None:
            axis_index = self.axis_index
            low = self.component.position[axis_index]
            high = low + self.component.bounds[axis_index] - 1
            tmp = self._get_axis_coord(event)
            if tmp >= low and tmp <= high:
                new_edge = self.mapper.map_data(self._get_axis_coord(event))
                #new_edge = self._map_data(self._get_axis_coord(event))
                if self._drag_edge == "high":
                    low_val = self.selection[0]
                    if new_edge >= low_val:
                        self.selection = (low_val, new_edge)
                    else:
                        self.selection = (new_edge, low_val)
                        self._drag_edge = "low"
                else:
                    high_val = self.selection[1]
                    if new_edge <= high_val:
                        self.selection = (new_edge, high_val)
                    else:
                        self.selection = (high_val, new_edge)
                        self._drag_edge = "high"

                self.component.request_redraw()
        return

class SingleLineTool(AbstractOverlay,DragTool):
    component = Instance(Component)
    line = Instance(Line, args=())
    points = List #data space points
    startpoint = Array
    endpoint = Array
    currpoint = Array
    color = ColorTrait
    
    finished = Event
    
    event_state=Enum('normal','dragging')
    
    def __init__(self,component,**traits):
        DragTool.__init__(self,component)
        AbstractOverlay.__init__(self,component)
        
    
    def overlay(self, component, gc, view_bounds=None, mode="normal"):
        drawfunc = getattr(self,self.event_state+'_draw')
           
        gc.save_state()
        gc.clip_to_rect(component.x, component.y, component.width-1, component.height-1)
        drawfunc(gc)
        gc.restore_state()
        
    def drag_start(self, event):
        self.startpoint = self.component.map_data(np.array((event.x,event.y)))
        self.currpoint = self.startpoint 
        self.event_state = 'dragging'
        self.request_redraw()
        
    def drag_end(self, event):
        self.endpoint = self.component.map_data(np.array((event.x,event.y)))
        self.currpoint = tuple()
        self.finished = np.array((self.startpoint,self.endpoint)) 
        self.event_state = 'normal'
        self.request_redraw()
    
    def drag_cancel(self, event):
        self.endpoint = tuple() #empty
        self.event_state = 'normal'
        self.request_redraw()
        
    def dragging(self, event):
        self.currpoint = self.component.map_data(np.array((event.x,event.y))) 
        self.request_redraw()
    
    def normal_draw(self,gc):
        return
        
    def dragging_draw(self,gc):
        self.line.points = list(self.component.map_screen((self.startpoint,self.currpoint)))
        self.line._draw(gc)
        
    def _color_changed(self):
        self.line.line_color = self.color
        self.line.vertex_color = self.color
        
class LineHighlighterOverlay(AbstractOverlay):
    
    extents = Array(shape=(None,2),dtype=float)
    continuua = List
    selected = Trait((None,Int))
    
    highlightcolor = ColorTrait((1.0,0,0,0.4))
    edgecolor = ColorTrait('black')
    selcolor = ColorTrait((1.0,0,0,0.7))
    selectedline = Trait((None,Int))
    
    def overlay(self, component, gc, view_bounds=None, mode="normal"):
        from operator import isSequenceType
        
        gc.save_state()
        gc.clip_to_rect(component.x, component.y, component.width-1, component.height-1)
        
        for i,((ex1,ex2),cont) in enumerate(zip(self.extents,self.continuua)):
            
            if self.selected is not None and i== self.selected:
                gc.set_fill_color(self.selcolor_)
            else:
                gc.set_fill_color(self.highlightcolor_)
            
            x1 = self.component.index_mapper.map_screen(ex1) 
            x2 = self.component.index_mapper.map_screen(ex2) 
            y1 = self.component.y
            y2 = self.component.y + component.height - 1
            
            gc.rect(x1,y1,x2-x1,y2-y1)
            gc.fill_path()
            
            
            gc.begin_path()
            gc.move_to(x1,y1)
            gc.line_to(x1,y2)
            gc.move_to(x2,y2)
            gc.line_to(x2,y1)
            gc.set_stroke_color(self.edgecolor_)
            gc.set_line_width(2 if i == self.selectedline else 1)
            gc.draw_path()
            
            if isSequenceType(cont) and len(cont)==2:
                yc1 = self.component.value_mapper.map_screen(cont[0])
                yc2 = self.component.value_mapper.map_screen(cont[1])
                gc.begin_path()
                gc.move_to(x1,yc1)
                gc.line_to(x2,yc2)
                r,g,b,a = self.selcolor_
                gc.set_stroke_color((r,g,b,1.0))
                gc.set_line_width(2)
                gc.draw_path()
            
        gc.restore_state()

class FeaturesAdapter(TabularAdapter):
    columns = [('Extent','extent'),
               ('Line Center','center'),
               ('Flux','flux'),
               ('EW','ew'),
               ('IDed Name','idname')]
               
    can_edit = Bool(False) #TODO:why doesn't this work?
               
    extent_text = Property
    extent_can_edit = Bool(True)
    
    def _get_extent_text(self):
        return '%.6g,%.6g'%self.content
    
    def _set_extent_text(self,val):
        vals = tuple([float(v) for v in val.split(',')])
        setattr(self.item,self.column_id,vals if vals[0] <= vals[1] else vals[::-1])
        self.object._update_featurelist()
        


class SpylotHandler(Handler):
    def nextSpec(self,info):
        info.object.specright = True
    
    def prevSpec(self,info):
        info.object.specleft = True
        
    def close(self,info, is_ok):
        info.object.lastspec = None
        return True
    
    def object_editfeatures_changed(self,info):
        #event handler placed here to pass in the correct parent
        info.object.edit_traits(view='features_view',parent=info.ui.control)
    
spylotkeybindings = KeyBindings(
                       KeyBinding(binding1='<',binding2='Ctrl-,',
                                  description='Previous spectrum',      
                                  method_name='prevSpec'),
                       KeyBinding(binding1='>',binding2='Ctrl-.',
                                  description='Next spectrum',      
                                  method_name='nextSpec')
                    )

class Spylot(HasTraits):
    defaultlines = 'galaxy' #can be 'galaxy' or 'stellar' or None
    
    plot = Instance(Plot)
    histplot = Bool(True)
    
    majorlineeditor = LineListEditor
    minorlineeditor = LineListEditor
    linenamelist = Property(Tuple)
    showmajorlines = Bool(True)
    showminorlines = Bool(False)
    
    labels = List(DataLabel)
    showlabels = Bool(False)
    editmajor = Button('Edit Major')
    editminor = Button('Edit Minor')
    
    specs = List(Instance(spec.Spectrum))
    currspeci = Int
    currspecip1 = Property(depends_on='currspeci')
    lowerspecip1 = Property(depends_on='currspeci')
    upperspecip1 = Property(depends_on='currspeci')
    currspec = Property
    lastspec = Instance(spec.Spectrum)
    z = Float
    lowerz = Float(0.0)
    upperz = Float(1.0)
    coarserz = Button('Coarser')
    finerz = Button('Finer')
    _zql,_zqh = min(spec.Spectrum._zqmap.keys()),max(spec.Spectrum._zqmap.keys())
    zqual = Range(_zql,_zqh,-1)
    
    spechanged = Event
    specleft = Button('<')
    specright = Button('>')
    
    scaleerr = Bool(False)
    scaleerrfraclow = Range(0.0,1.0,1.0)
    scaleerrfrachigh = Float(1.0)
    fluxformat = Button('Flux Line Format')
    errformat = Button('Error Line Format')
    showcoords = Bool(False)
    showgrid = Bool(True)
    
    dosmoothing = Bool(False)
    smoothing = Float(3)
    
    contsub = Button('Fit Continuum...')
    contclear = Button('Clear Continuum')
    showcont = Bool(False)
    contformat = Button('Continuum Line Format')
    
    _titlestr = Str('Spectrum 0/0')
    _oldunit = Str('')
    maintool = Tuple(Instance(Interactor),Instance(AbstractOverlay))
    
    featureselmode = Enum(['No Selection','Click Select','Range Select','Base Select','Click Delete'])
    editfeatures = Button('Features...')
    showfeatures = Bool(True)
    featurelocsmooth = Float(None)
    featurelocsize = Int(200)
    featurelist = List(Instance(spec.SpectralFeature))
    
    selectedfeatureindex = Int
    deletefeature = Button('Delete')
    idfeature = Button('Identify')
    recalcfeature = Button('Recalculate')
    clearfeatures = Button('Clear')
    
    delcurrspec = Button('Delete Current')
    saveloadfile = File(filter=['*.specs'])
    savespeclist = Button('Save Spectra')
    loadspeclist = Button('Load Spectra')
    loadaddfile = File(filter=['*.fits'])
    loadaddspec = Button('Add Spectrum')
    loadaddspectype = Enum('wcs','deimos','astropysics')
    
    titlegroup = HGroup(Item('specleft',show_label=False,enabled_when='currspeci>0'),
                         spring,
                         Label('Spectrum #',height=0.5),
                         Item('currspecip1',show_label=False,editor=RangeEditor(low_name='lowerspecip1',high_name='upperspecip1',mode='spinner')),
                         Item('_titlestr',style='readonly',show_label=False),
                         spring,
                         Item('specright',show_label=False,enabled_when='currspeci<(len(specs)-1)'))
                         
    speclistgroup = HGroup(Label('Spectrum List:'),spring,
                           Item('delcurrspec',show_label=False,enabled_when='len(specs)>1'),
                           Item('saveloadfile',show_label=False), 
                           Item('savespeclist',show_label=False),
                           Item('loadspeclist',show_label=False,enabled_when='os.path.exists(saveloadfile)'),
                           Item('loadaddfile',show_label=False), 
                           Item('loadaddspec',show_label=False,enabled_when='os.path.exists(saveloadfile)'),
                           Item('loadaddspectype',show_label=False),
                           spring)
                    
    plotformatgroup = HGroup(spring,
                             Item('fluxformat',show_label=False),
                             Item('errformat',show_label=False),
                             Item('scaleerr',label='Scale Error?'),
                             Item('scaleerrfraclow',label='Lower',enabled_when='scaleerr',editor=TextEditor(evaluate=float)),
                             Item('scaleerrfrachigh',label='Upper',enabled_when='scaleerr'),
                             Item('showgrid',label='Grid?'),
                             Item('showcoords',label='Coords?'),
                             spring)
                    
    featuregroup = HGroup(spring,
                          Item('showmajorlines',label='Show major?'),
                          Item('editmajor',show_label=False),
                          Item('showlabels',label='Labels?'),
                          Item('showminorlines',label='Show minor?'),
                          Item('editminor',show_label=False),
                          spring,
                          Item('editfeatures',show_label=False),
                          Item('featureselmode',show_label=False),
                          spring)
                          
    continuumgroup = HGroup(spring,
                            Item('contsub',show_label=False),
                            Item('contclear',show_label=False),
                            Item('showcont',label='Continuum line?'),
                            Item('contformat',show_label=False),
                            Item('dosmoothing',label='Smooth?'),
                            Item('smoothing',show_label=False,enabled_when='dosmoothing'),
                            spring)
                    
    zgroup = VGroup(Item('z',editor=RangeEditor(low_name='lowerz',high_name='upperz',format='%.4f')),
                    HGroup(Item('lowerz',show_label=False),
                           Item('coarserz',show_label=False),
                           spring,
                           Item('zqual',style='custom',label='Z quality',editor=RangeEditor(cols=_zqh-_zql+1,low=_zql,high=_zqh)),
                           spring,
                           Item('finerz',show_label=False),
                           Item('upperz',show_label=False))
                   )
                   
    features_view = View(VGroup(HGroup(Item('showfeatures',label='Show'),
                                       Item('featurelocsmooth',label='Locator Smoothing'),
                                       Item('featurelocsize',label='Locator Window size')
                                       ),
                                Item('featurelist',
                                     editor=TabularEditor(adapter=FeaturesAdapter(),
                                     selected_row='selectedfeatureindex'),
                                     show_label=False),
                                HGroup(Item('deletefeature',show_label=False,enabled_when='len(featurelist)>0 and selectedfeatureindex>=0'),
                                       Item('idfeature',show_label=False,enabled_when='len(featurelist)>0 and selectedfeatureindex>=0'),
                                       Item('recalcfeature',show_label=False,enabled_when='len(featurelist)>0 and selectedfeatureindex>=0'),
                                       Item('clearfeatures',show_label=False,enabled_when='len(featurelist)>0')
                                )),
                    resizable=True, 
                    title='Spylot Features')
                   
    traits_view = View(VGroup(Include('titlegroup'),
                              Include('speclistgroup'),
                              Include('plotformatgroup'),
                              Include('featuregroup'),     
                              Include('continuumgroup'),  
                              Item('plot',editor=ComponentEditor(),show_label=False,width=768),
                              Include('zgroup')
                           ),
                           resizable=True, 
                           title='Spectrum Plotter',
                           handler = SpylotHandler(),
                           key_bindings = spylotkeybindings
                       )
                      
    def __init__(self,specs,**traits):
        #pd = ArrayPlotData(x=[1],x0=[1],flux=[1],err=[1]) #reset by spechanged event
        pd = ArrayPlotData(x=[1],flux=[1],err=[1]) #reset by spechanged event
        pd.set_data('majorx',[1,1])#reset by majorlines change
        pd.set_data('majory',[0,1])#reset by majorlines change
        pd.set_data('minorx',[1,1])#reset by minorlines change
        pd.set_data('minory',[0,1])#reset by minorlines change
        pd.set_data('continuum',[0,0])#reset
        
        self.plot = plot = Plot(pd,resizeable='hv')
        ploi = plot.draw_order.index('plot')
        plot.draw_order[ploi:ploi] = ['continuum','err','flux','annotations']
        plot.plot(('x','flux'),name='flux',type='line',line_style='solid',color='blue',draw_layer='flux',unified_draw=True)
        plot.plot(('x','err'),name='err',type='line',line_style='dash',color='green',draw_layer='err',unified_draw=True)
        
        topmapper = LinearMapper(range=DataRange1D())
        plot.x_mapper.range.on_trait_change(self._update_upperaxis_range,'updated')
        plot.x_mapper.on_trait_change(self._update_upperaxis_screen,'updated')
        self.upperaxis = PlotAxis(plot,orientation='top',mapper=topmapper)
        plot.overlays.append(self.upperaxis)
        
        self.errmapperfixed = plot.plots['err'][0].value_mapper
        self.errmapperscaled = LinearMapper(range=DataRange1D(high=1.0,low=0))
        plot.x_mapper.on_trait_change(self._update_errmapper_screen,'updated')
        
        plot.padding_top = 30 #default is a bit much
        plot.padding_left = 70 #default is a bit too little
        
        majorlineplot = plot.plot(('majorx','majory'),name='majorlineplot',type='line',line_style='dash',color='red')[0]
        majorlineplot.set(draw_layer='annotations',unified_draw=True)
        majorlineplot.value_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
        majorlineplot.visible = self.showmajorlines
        del plot.x_mapper.range.sources[-1]  #remove the line plot from the x_mapper sources so scaling is only on the spectrum
        self.majorlineeditor = LineListEditor(lineplot=majorlineplot)
        
        minorlineplot = plot.plot(('minorx','minory'),name='minorlineplot',type='line',line_style='dot',color='red')[0]
        minorlineplot.set(draw_layer='annotations',unified_draw=True)
        minorlineplot.value_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
        minorlineplot.visible = self.showminorlines
        del plot.x_mapper.range.sources[-1]  #remove the line plot from the x_mapper sources so scaling is only on the spectrum
        self.minorlineeditor = LineListEditor(lineplot=minorlineplot)
        
        self.contline = plot.plot(('x','continuum'),name='continuum',type='line',line_style='solid',color='black')[0]
        self.contline.set(draw_layer='continuum',unified_draw=True)
        self.contline.visible = self.showcont
#        idat = ArrayDataSource((0.0,1.0))
#        vdat = ArrayDataSource((0.0,0.0))
#        self.zeroline = LinePlot(index=idat,value=vdat,line_style='solid',color='black')
#        self.zeroline.index_mapper = LinearMapper(range=DataRange1D(high=0.9,low=0.1))
#        self.zeroline.value_mapper = self.plot.y_mapper
#        self.zeroline.visible = self.showcont
#        self.plot.add(self.zeroline)
        
        if Spylot.defaultlines:
            defaultlines = _get_default_lines(Spylot.defaultlines)
            self.majorlineeditor.candidates = defaultlines[0]
            self.majorlineeditor.selectednames = defaultlines[1]
            self.minorlineeditor.candidates = defaultlines[0]
            self.minorlineeditor.selectednames = defaultlines[2]
        
        plot.tools.append(PanTool(plot))
        plot.tools.append(ZoomTool(plot))
        plot.overlays.append(plot.tools[-1])
        
        self.coordtext = TextBoxOverlay(component=plot,align='ul')
        plot.overlays.append(self.coordtext)
        plot.tools.append(MouseMoveReporter(overlay=self.coordtext,plot=plot))
        self.coordtext.visible = self.showcoords
        
        self.linehighlighter = lho = LineHighlighterOverlay(component=plot)
        lho.visible = self.showfeatures
        plot.overlays.append(lho)
         
        if specs is None:
            specs = []
        elif isinstance(specs,spec.Spectrum):
            specs = [specs]
        self.specs = specs
        
        self.spechanged = True
        
        self.on_trait_change(self._majorlines_changed,'majorlineeditor.selectedobjs')
        self.on_trait_change(self._minorlines_changed,'minorlineeditor.selectedobjs')
        
        super(Spylot,self).__init__(**traits)
        
    def __del__(self):
        try:
            self.currspec._features = list(self.currspec._features)
        except (AttributeError,IndexError),e:
            pass
        #no __del__ present in parents?
        #finally:
        #    super(Spylot,self).__del__()
        
    def _update_upperaxis_range(self):
        newlow = self.plot.x_mapper.range.low/(self.z+1)
        newhigh = self.plot.x_mapper.range.high/(self.z+1)
        self.upperaxis.mapper.range.set_bounds(newlow,newhigh)
        
    def _update_upperaxis_screen(self):
        self.upperaxis.mapper.screen_bounds = self.plot.x_mapper.screen_bounds
        
    def _update_errmapper_screen(self):
        self.errmapperscaled.screen_bounds = self.plot.y_mapper.screen_bounds
        
    def _get_currspec(self):
        return self.specs[self.currspeci]
    def _set_currspec(self,val):
        self.currspeci = self.specs.index(val)
        
    def _get_lowerspecip1(self):
        return 1
    def _get_upperspecip1(self):
        return len(self.specs)
    
    def _get_currspecip1(self):
        return self.currspeci+1
    def _set_currspecip1(self,val):
        self.currspeci = val - 1
    
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
        
    def _coarserz_changed(self):
        rng = self.upperz - self.lowerz
        newrng = rng*2
        lowerz = self.z-newrng/2
        if lowerz > 0:
            self.lowerz = lowerz
            self.upperz = self.z+newrng/2
        else:
            self.lowerz = 0
            self.upperz = newrng
    
    def _finerz_changed(self):
        rng = self.upperz - self.lowerz
        newrng = rng/2
        lowerz = self.z-newrng/2
        if lowerz > 0:
            self.lowerz = lowerz
            self.upperz = self.z+newrng/2
        else:
            self.lowerz = 0
            self.upperz = newrng
    
    def _z_changed(self):
        self.currspec.z = self.z
        #self.plot.data.set_data('x0',_hist_sample_x(self.currspec.x0) if self.histplot else self.currspec.x0)
        self._majorlines_changed()
        self._minorlines_changed()
        self._update_upperaxis_range()
        
    def _zqual_changed(self):
        self.currspec.zqual = self.zqual
        
    def _scaleerr_changed(self,new):
        if new:
            self._scaleerrfrac_changed()
            self._update_errmapper_screen() #not sure why this is necessary, but apparently it is
        self.spechanged = True #TODO: control more finely?
        
    @on_trait_change('scaleerrfraclow,scaleerrfrachigh')
    def _scaleerrfrac_changed(self):
        if self.scaleerrfraclow==1:
            self.errmapperscaled.range.low = self.plot.data.get_data('err').min()
        else:
            self.errmapperscaled.range.low = self.scaleerrfraclow*self.plot.data.get_data('err').max()
        self.errmapperscaled.range.high = self.scaleerrfrachigh*self.plot.data.get_data('err').max()
        if self.scaleerr:
            self._scaleerr_changed(False)
    
    @on_trait_change('smoothing,dosmoothing')    
    def _update_smoothing(self):
        self.spechanged = True #TODO: control more finely?
        
    def _spechanged_fired(self):
        p = self.plot
        s = self.currspec
        self.lastspec = s #performs any necessary cleanup
        
        if s.unit != self._oldunit:
            self._minorlines_changed()
            self._majorlines_changed()
            self._oldunit = s.unit
        
        x = s.x
        #x0 = s.x0
        if self.dosmoothing:
            smoothing = self.smoothing
            if smoothing < 0:
                smoothing *=-1
                smtype = 'boxcar'
            else:
                smtype = 'gaussian'
            flux,err = s.smooth(smoothing,smtype,replace=False)
        else:
            flux = s.flux.copy()
            err = s.err.copy()
            
        if self.histplot:
            x = _hist_sample_x(x)
            #x0 = _hist_sample_x(x0)
            flux = _hist_sample_y(flux)
            err = _hist_sample_y(err)
        
        if self.scaleerr:
            p.plots['err'][0].value_mapper = self.errmapperscaled
        else:
            p.plots['err'][0].value_mapper = self.errmapperfixed
           
        flux[~np.isfinite(flux)] = 0
        err[~np.isfinite(err)] = 0
            
        p.data.set_data('x',x)
        #p.data.set_data('x0',x0)
        p.data.set_data('flux',flux)
        p.data.set_data('err',err)
        cmodel = self.currspec.continuum
        p.data.set_data('continuum',np.zeros_like(x) if cmodel is None else cmodel(x))
        
        rng = self.upperz - self.lowerz
        self.z = z = s.z
        if (z - rng/2.0) >=0:
            self.upperz,self.lowerz = z+rng/2.0,z-rng/2.0
        else:
            self.lowerz = 0
            self.upperz = rng
        self.zqual = s.zqual
        
        #must have the SAME feature list instance in the spylot gui and the Spectrum
        self.featurelist = s._features
        s._features = self.featurelist 
        
        units = tuple(s.unit.split('-'))
        p.x_axis.title = '%s/%s'%units
        self.upperaxis.title = 'rest %s/%s'%units
        p.y_axis.title ='Flux/(erg s^-1 cm^-2 %s^-1)'%units[1]
        p.request_redraw()
        
        #self._titlestr = 'Spectrum %i/%i: %s'%(self.currspeci+1,len(self.specs),s.name)
        self._titlestr = '/%i: %s'%(len(self.specs),s.name)
        
    def _lastspec_changed(self,old,new):
        if old is not None:
            #necessary to not carry arround traited feature lists
            old._features = list(old._features)
        
    def _specleft_fired(self):
        if self.currspeci > 0:
            self.currspeci -= 1
    
    def _specright_fired(self):
        if self.currspeci < len(self.specs)-1:
            self.currspeci += 1
            
    def _showlabels_changed(self):
        if len(self.labels) == 0:
            self._rebuild_labels()
        else:
            for l in self.labels:
                l.visible = self.showlabels
        self.plot.request_redraw()
        
    def _rebuild_labels(self):
        for dl in self.labels:
                if dl in self.plot.overlays:
                    self.plot.overlays.remove(dl)
        self.labels = []
        
        if self.showlabels:
            majorlines = self.majorlineeditor.selectedobjs
            linestrs = [str(l) for l in majorlines]
            linelocs = [l.loc for l in majorlines]
            lfact = self.z+1
            xd,yd = self.plot.data.get_data('x'),self.plot.data.get_data('flux')
            for s,l in zip(linestrs,linelocs):
                x = l*lfact
                y = np.interp(x,xd,yd)
#                xs,ys = self.plot.map_screen((x,y))
#                x,y = self.plot.map_data((xs,ys+50))
                dl = DataLabel(self.plot,data_point=(x,y),lines=[s],
                               label_position="top", padding_bottom=10,
                               arrow_visible=False)
                self.labels.append(dl)
                self.plot.overlays.append(dl)
                
        

    def _majorlines_changed(self):
        majorlines = self.majorlineeditor.selectedobjs
        if len(majorlines)>0:
            u = self.currspec.unit
            linelocs = [l.getUnitLoc(u) for l in majorlines]
            linex = np.repeat(linelocs,2)
            linex.sort()
            liney = np.tile((0,1,1,0),np.ceil(len(majorlines)/2))[:len(majorlines)*2]
            self.plot.data.set_data('majorx',linex*(self.z+1))
            self._rebuild_labels()
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
        
        
    #TODO:possibly cache
    def _get_linenamelist(self):
        mjos = self.majorlineeditor.selectedobjs
        mnos = self.minorlineeditor.selectedobjs
        lst = []
        for o in self.majorlineeditor.selectedobjs:
            lst.append(o.name)
        sepi = len(lst)
        for o in self.minorlineeditor.selectedobjs:
            lst.append(o.name)
        lst.insert(sepi,'-'*max([len(n) for n in lst]))
        return tuple(lst)
        
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
        
    def _contformat_fired(self):
        from copy import copy
        plot = self.plot.plots['continuum'][0]
        v = copy(plot.trait_view())
        v.title = 'Continuum Line Format'
        plot.edit_traits(view=v)
        
    def _showgrid_changed(self,new):
        p = self.plot
        p.x_grid.visible = new
        p.y_grid.visible = new
        p.request_redraw()
    
    def _showcoords_changed(self):    
        self.coordtext.visible = self.showcoords
        self.plot.request_redraw()
        
    def _showcont_changed(self,new):
        self.contline.visible = new
        self.plot.request_redraw()
        
    def _contsub_fired(self):
        self.currspec.fitContinuum(interactive='reuse')
        pd = self.plot.data
        cmodel = self.currspec.continuum
        x = pd.get_data('x')
        pd.set_data('continuum',np.zeros_like(x) if cmodel is None else cmodel(x))

    def _contclear_fired(self):
        self.currspec.continuum = None
        x = self.plot.data.get_data('x')
        self.plot.data.set_data('continuum',np.zeros_like(x))
        
    def _maintool_changed(self,old,new):
        otool,oover = old
        ntool,nover = new
        
        if otool is not None:
            self.plot.tools.remove(otool)
        if oover is not None:
            self.plot.overlays.remove(oover)
        if ntool is not None:
            self.plot.tools.append(ntool)
        if nover is not None:
            self.plot.overlays.append(nover)
        
    def _featureselmode_changed(self,old,new):
        pl = self.plot
            
        if new == 'No Selection':
            pl.tools[0].drag_button = 'left'
            self.maintool = (None,None)

        elif new == 'Click Select':
            pl.tools[0].drag_button = 'right'
            fctool = FeatureClickTool(plot=self.plot)
            self.maintool = (fctool,None)
            fctool.on_trait_change(self._add_feature_point,'mouse_clicked')
            
        elif new == 'Range Select':
            pl.tools[0].drag_button = 'right'
            fplot = pl.plots['flux'][0]
            
            #TODO:clean up these hacky bugfix techniques if possible
            rstool = RangeSelectionBugfix(component=fplot,
                                          left_button_selects=True,
                                          enable_resize=False,
                                          mapper=self.plot.x_mapper)
            rstool.plot = fplot
            rstool.component = pl
            
            rsoverlay = RangeSelectionOverlay(component=fplot,
                                              mapper=self.plot.x_mapper)
            rsoverlay.plot = fplot
            rsoverlay.component = pl
            self.maintool = (rstool,rsoverlay)
            rstool.on_trait_change(self._add_feature_range,'selection_completed')
            
        elif new == 'Base Select':
            pl.tools[0].drag_button = 'right'
            slt = SingleLineTool(component=self.plot)
            slt.component = self.plot
            self.maintool = (slt,slt)
            slt.on_trait_change(self._add_feature_base,'finished')
        elif new == 'Click Delete':
            pl.tools[0].drag_button = 'right'
            fctool = FeatureClickTool(plot=self.plot)
            self.maintool = (fctool,None)
            fctool.on_trait_change(self._del_feature_point,'mouse_clicked')    
        else:
            assert True,'enum invalid!'
        
        
            
    def _add_feature_range(self,selection):
        rstool = self.maintool[0]
        x1,x2 = selection
        rstool.deselect(Event(window=rstool.downwindow))
        self.currspec.addFeatureRange(x1,x2)
        
    def _add_feature_point(self,dataxy):
        datax,datay = dataxy
        smooth = self.featurelocsmooth if self.featurelocsmooth != 0 else None
        window = self.featurelocsize
        self.currspec.addFeatureLocation(datax,window=window,smooth=smooth)
            
    def _del_feature_point(self,dataxy):
        datax,datay = dataxy
        self.currspec.removeFeatureLocation(datax)
        
    def _add_feature_base(self,dataarr):
        (x1,y1),(x2,y2) = dataarr
        self.currspec.addFeatureRange(x1,x2,continuum=(y1,y2)) 
        
            
##    Moved to handler
#    def _editfeatures_fired(self):
#        self.edit_traits(view='features_view')

    def _deletefeature_fired(self):
        del self.featurelist[self.selectedfeatureindex]
        
    def _idfeature_fired(self):
        raise NotImplementedError
    
    def _recalcfeature_fired(self):
        sf = self.featurelist[self.selectedfeatureindex]
        sf.computeFeatureData(self.currspec,edges='interp')
        #TODO:figure out how to make the tabulareditor update
    
    def _clearfeatures_fired(self):
        for i in range(len(self.featurelist)):
            del self.featurelist[0]
        
    def _selectedfeatureindex_changed(self):
        i = self.selectedfeatureindex
        if i >=0:
            self.linehighlighter.selectedline = i
        else:
            self.linehighlighter.selectedline = None
        if self.showfeatures:
            self.plot.request_redraw()
        
    def _showfeatures_changed(self):
        self.linehighlighter.visible = self.showfeatures
        self.linehighlighter.request_redraw()
    
    @on_trait_change('featurelist_items,featurelist')
    def _update_featurelist(self):
        self.linehighlighter.extents = [f.extent for f in self.featurelist] if  len(self.featurelist)>0 else np.ndarray((0,2))
        self.linehighlighter.continuua = [f.continuum for f in self.featurelist]
        if self.showfeatures:
            self.plot.request_redraw()
            
    def _saveloadfile_changed(self,new):
        if not os.path.exists(new):
            if not new.endswith('.specs'):
                self.saveloadfile = new+'.specs'
                
    def _delcurrspec_fired(self):
        del self.specs[self.currspeci]
        
        self.spechanged = True
        
        if self.currspeci == len(self.specs) and self.currspeci > 0:
            self.currspeci -= 1
    
    def _savespeclist_fired(self):
        from cPickle import dump
        
        fn = self.saveloadfile
        with open(fn,'w') as f:
            dump(list(self.specs),f)
        
    
    def _loadspeclist_fired(self):
        from cPickle import load
        from operator import isSequenceType
        
        fn = self.saveloadfile
        with open(fn) as f:
            obj = load(f)
            
        if isSequenceType(obj):
            self.specs = obj
        else:
            self.specs = [obj]
        self.spechanged = True #doesn't always fire with the assignment above
        
    def _loadaddspec_fired(self):
        if loadaddspec == 'wcs':
            spec = spec.load_wcs_spectrum(self.loadaddfile)
        elif loadaddspec == 'deimos':
            spec = spec.load_deimos_spectrum(self.loadaddfile)
        elif loadaddspec == 'astropysics':
            spec = spec.Spectrum.load(self.loadaddfile)
            
        self.specs.append(spec)
        self.loadaddfile = ''
            
def _get_default_lines(linetypes):
    candidates = spec.load_line_list(linetypes)
    if linetypes == 'galaxy':
        major=['Ly_alpha','H_alpha','H_beta','H_gamma','[OIII]_5007','[OIII]_4959','[OII]_3727','Balmer_limit','[SII]_6717','[SII]_6731']
        return candidates,[str(c) for c in candidates if c.name in major],[str(c) for c in candidates if c.name not in major]
    elif linetypes == 'stellar':
        major = ['Halpha']
        return candidates,[str(c) for c in candidates if c.name in major],[str(c) for c in candidates if c not in major]
    else:
        return candidates,[],[]
    
def spylot_specs(specs):
    """
    Generates a Spylot instance containing the supplied sequence of spectra and 
    displays it.  A GUI application instance must already exist (e.g. 
    interactive mode of ipython)
    
    returns the Spylot instance
    """
    sp = Spylot(specs)
    sp.edit_traits()
    return sp
