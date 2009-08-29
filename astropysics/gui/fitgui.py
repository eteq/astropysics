#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the FitGui gui.
"""
from __future__ import division,with_statement
import numpy as np

from enthought.traits.api import HasTraits,Instance,Int,Float,Bool,Button, \
                                 Event,Property,on_trait_change,Array,List, \
                                 Tuple,Str,Dict,cached_property,Color,Enum, \
                                 TraitError
from enthought.traits.ui.api import View,Handler,Item,Label,Group,VGroup, \
                                    HGroup, InstanceEditor,EnumEditor, \
                                    ListEditor, TupleEditor
from enthought.traits.ui.menu import ModalButtons
from enthought.chaco.api import Plot,ArrayPlotData,jet,ColorBar,HPlotContainer,\
                                ColorMapper,LinearMapper,ScatterInspectorOverlay,\
                                LassoOverlay,AbstractOverlay
from enthought.chaco.tools.api import PanTool, ZoomTool,SelectTool,LassoSelection,\
                                      ScatterInspector
from enthought.enable.api import ColorTrait,ComponentEditor

#TODO:move to remove deps?
from enthought.tvtk.pyface.scene_editor import SceneEditor 
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene


from ..models import list_models,get_model,Model,binned_weights

class ColorMapperFixSingleVal(ColorMapper):
    coloratval = ColorTrait('black')
    val = 0
    
    def map_screen(self, data_array):
        res = super(ColorMapperFixSingleVal,self).map_screen(data_array)
        res[data_array==self.val] = self.coloratval_
        return res

#_cmap = jet
def _cmapblack(range, **traits):
    _data =   {'red':   ((0,1,1), (0.3, .8, .8), (0.5, 0, 0), (0.75,0.75, 0.75),(.875,.2,.2),
                         (1, 0, 0)),
               'blue': ((0., 0, 0), (0.3,0, 0), (0.5,0, 0), (0.75,.75, .75),
                         (0.875,0.75,0.75), (1, 1, 1)),
               'green':  ((0.,0, 0),(0.3,.8,.8), (0.4, 0.4, 0.4),(0.5,1,1), (0.65,.75, .75), (0.75,0.1, 0.1),
                         (1, 0, 0))}

    return ColorMapperFixSingleVal.from_segment_map(_data, range=range, **traits)

def _cmap(range, **traits):
    _data =   {'red':   ((0,1,1), (0.3, .8, .8), (0.5, 0, 0), (0.75,0.75, 0.75),(.875,.2,.2),
                         (1, 0, 0)),
               'blue': ((0., 0, 0), (0.3,0, 0), (0.5,0, 0), (0.75,.75, .75),
                         (0.875,0.75,0.75), (1, 1, 1)),
               'green':  ((0.,0, 0),(0.3,.8,.8), (0.4, 0.4, 0.4),(0.5,1,1), (0.65,.75, .75), (0.75,0.1, 0.1),
                         (1, 0, 0))}

#    """ inverted version of 'jet' colormap"""

#    _data =   {'red':   ((0., 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89,1, 1),
#                         (1, 0.5, 0.5)),
#               'green': ((0., 0, 0), (0.125,0, 0), (0.375,1, 1), (0.64,1, 1),
#                         (0.91,0,0), (1, 0, 0)),
#               'blue':  ((0., 0.5, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.65,0, 0),
#                         (1, 0, 0))}
#    for k,v in _data.items():
#        _data[k] = tuple(reversed([(v[-1-i][0],t[1],t[2]) for i,t in enumerate(v)]))
#    
    return ColorMapper.from_segment_map(_data, range=range, **traits)

class _TraitedModel(HasTraits):
    from inspect import isclass
    
    model = Instance(Model,allow_none=True)
    modelname = Property(Str)
    updatetraitparams = Event
    paramchange = Event
    fitdata = Event
    
    def __init__(self,model,**traits):
        super(_TraitedModel,self).__init__(**traits)
        
        from inspect import isclass
        
        if isinstance(model,basestring):
            model = get_model(model)()
        elif isclass(model):
            model = model()
            
        self.model = model
        
    def default_traits_view(self):
        if self.model is None:
            g = Group()
            g.content.append(Label('No Model Selected'))
        else:
            #g = Group(label=self.modelname,show_border=False,orientation='horizontal',layout='flow')
            g = Group(label=self.modelname,show_border=True,orientation='horizontal',scrollable=True)
            for p in self.model.params:
                gi = Group(orientation='horizontal',label=p)
                self.add_trait(p,Float)
                setattr(self,p,getattr(self.model,p))
                self.on_trait_change(self._param_change_handler,p)
                gi.content.append(Item(p,show_label=False))
                
                ffp = 'fixfit_'+p
                self.add_trait(ffp,Bool)
                setattr(self,ffp,False)
                self.on_trait_change(self._param_change_handler,ffp)
                gi.content.append(Item(ffp,label='Fix?'))
                
                g.content.append(gi)
            
        return View(g,buttons=['Apply','Revert','OK','Cancel'])
    
    def _param_change_handler(self,name,new):
        setattr(self.model,name,new)
        self.paramchange = name
    
    def _updatetraitparams_fired(self):
        m = self.model
        for p in m.params:
            setattr(self,p,getattr(m,p))
        self.paramchange = True
        
    def _fitdata_fired(self,new):
        from operator import isSequenceType,isMappingType
        
        if self.model is not None:
            if isSequenceType(new) and len(new) == 2:
                kw={'x':new[0],'y':new[1]}
            elif isMappingType(new):
                kw = dict(new)
                for xy in ('x','y'):
                    if xy not in new:
                        if self.model.fiteddata:
                            new[xy] = self.model.fiteddata[0]
                        else:
                            raise ValueError('not pre-fitted data available')
            elif new is True:
                if self.model.fiteddata:
                    fd = self.model.fiteddata
                    kw= {'x':fd[0],'y':fd[1]}
                else:
                    raise ValueError('No data to fit')
            else:
                raise ValueError('unusable fitdata event input')
            
            if 'fixedpars' not in kw:
                 kw['fixedpars'] = [tn.replace('fixfit_','') for tn in self.traits() if tn.startswith('fixfit_') if getattr(self,tn)]
            
            self.model.fitData(**kw)
            self.updatetraitparams = True
    
    def _get_modelname(self):
        if self.model is None:
            return 'None'
        else:
            return self.model.__class__.__name__
        
class NewModelSelector(HasTraits):
    modelnames = List
    selectedname = Str('No Model')
    modelargnum = Int(2)
    selectedmodelclass = Property
    isvarargmodel = Property(depends_on='modelnames')
    
    traits_view = View(Item('selectedname',label='Model Name:',editor=EnumEditor(name='modelnames')),
                       Item('modelargnum',label='Extra Parameters:',enabled_when='isvarargmodel'),
                       buttons=['OK','Cancel'])
    
    def __init__(self,include_models=None,exclude_models=None,**traits):
        super(NewModelSelector,self).__init__(**traits)
        
        self.modelnames = list_models(1,include_models,exclude_models)
        self.modelnames.insert(0,'No Model')
        self.modelnames.sort()
        
    def _get_selectedmodelclass(self):
        n = self.selectedname
        if n == 'No Model':
            return None
        else:
            return get_model(n)
        
    def _get_isvarargmodel(self):
        cls = self.selectedmodelclass
        
        if cls is None:
            return False
        else:
            return cls._args is None
        
#class WeightFillOverlay(AbstractOverlay):
#    weightval = Float(0)
#    color = ColorTrait('black')
#    plot = Instance(Plot)
    
#    def overlay(self, component, gc, view_bounds=None, mode="normal"):
#        from enthought.chaco.scatterplot import render_markers
        
#        plot = self.component
#        scatter = plot.plots['data'][0]
#        if not plot or not scatter or not scatter.index or not scatter.value:
#            return
        
#        w = plot.data.get_data('weights')
#        inds = w==self.weightval
        
#        index_data = scatter.index.get_data()
#        value_data = scatter.value.get_data()
#        screen_pts = scatter.map_screen(np.array([index_data[inds],value_data[inds]]).T)
#        screen_pts = screen_pts+[plot.x,plot.y]
        
#        props = ('line_width','marker_size','marker')
#        markerprops = dict([(prop,getattr(scatter,prop)) for prop in props])
        
#        markerprops['color']=self.color_
#        markerprops['outline_color']=self.color_
        
#        if markerprops.get('marker', None) == 'custom':
#            markerprops['custom_symbol'] = scatter.custom_symbol
        
#        gc.save_state()
#        gc.clip_to_rect(scatter.x+plot.x, scatter.y+plot.y, scatter.width, scatter.height)
#        render_markers(gc, screen_pts, **markerprops)
#        gc.restore_state()

class FGHandler(Handler):
    def object_selbutton_changed(self,info):
        info.object.edit_traits(parent=info.ui.control,view='selection_view')
        
    def object_datasymb_changed(self,info):
        kind = info.ui.rebuild.__name__.replace('ui_','') #TODO:not hack!
        info.object.plot.plots['data'][0].edit_traits(parent=info.ui.control,
                                                      kind=kind)
    
    def object_modline_changed(self,info):
        kind = info.ui.rebuild.__name__.replace('ui_','') #TODO:not hack!
        info.object.plot.plots['model'][0].edit_traits(parent=info.ui.control,
                                                       kind=kind)
    
class FitGui(HasTraits):
    plot = Instance(Plot)
    colorbar = Instance(ColorBar)
    plotcontainer = Instance(HPlotContainer)
    tmodel = Instance(_TraitedModel,allow_none = False)
    nomodel = Property
    newmodel = Button('New Model...')
    fitmodel = Button('Fit Model')
    updatemodelplot = Button('Update Model Plot')
    autoupdate = Bool(True)
    data = Array(dtype=float,shape=(2,None))
    weights = Array
    weighttype = Enum(('custom','equal','lin bins','log bins'))
    weightsvary = Property(Bool)
    weights0rem = Bool(True)
    modelselector = NewModelSelector
    
    selbutton = Button('Selection...')    
    scattertool = Enum(None,'clicktoggle','clicksingle','clickimmediate','lassoadd','lassoremove','lassoinvert')
    selectedi = Property #indecies of the selected objects
    weightchangesel = Button('Set Selection To')
    weightchangeto = Float(1.0)
    delsel = Button('Delete Selected')
    unselectonaction = Bool(True)
    clearsel = Button('Clear Selections')
    lastselaction = Str('None')
    
    datasymb = Button('Data Symbol...')
    modline = Button('Model Line...')
    
    savews = Button('Save Weights')
    loadws = Button('Load Weights')
    _savedws = Array
    
    plotname = Property
    
    nmod = Int(1024)
    #modelpanel = View(Label('empty'),kind='subpanel',title='model editor') 
    modelpanel = View
    
    traits_view = View(VGroup(
                       Item('plotcontainer', editor=ComponentEditor(),show_label=False),
                       HGroup(VGroup(HGroup(Item('weighttype',label='Weights:'),
                                            Item('savews',show_label=False),
                                            Item('loadws',enabled_when='_savedws',show_label=False)),
                                Item('weights0rem',label='Remove 0-weight points for fit?'),
                                HGroup(Item('newmodel',show_label=False),
                                       Item('fitmodel',show_label=False),
                                       Item('selbutton',show_label=False))  
                              ),
                              VGroup(HGroup(Item('autoupdate',label='Auto?'),
                              Item('updatemodelplot',show_label=False,enabled_when='not autoupdate')),
                              Item('nmod',label='Nmodel'),
                              HGroup(Item('datasymb',show_label=False),Item('modline',show_label=False)))),
                       Item('tmodel',show_label=False,style='custom',editor=InstanceEditor(kind='subpanel'))
                      ),
                    handler=FGHandler(),
                    resizable=True, 
                    title='Data Fitting',
                    buttons=['OK','Cancel'],
                    width=700,
                    height=700
                    )
                    
    panel_view = View(VGroup(
                       Item('plot', editor=ComponentEditor(),show_label=False),
                       HGroup(Item('tmodel.modelname',show_label=False,style='readonly'),
                              Item('nmod',label='Number of model points'),
                              Item('updatemodelplot',show_label=False,enabled_when='not autoupdate'),
                              Item('autoupdate',label='Auto?'))
                      ),
                    title='Model Data Fitter'
                    )
                    
                    
    selection_view = View(Group(
                           Item('scattertool',label='Selection Mode',
                                 editor=EnumEditor(values={None:'1:No Selection',
                                                           'clicktoggle':'3:Toggle Select',
                                                           'clicksingle':'2:Single Select',
                                                           'clickimmediate':'7:Immediate',
                                                           'lassoadd':'4:Add with Lasso',
                                                           'lassoremove':'5:Remove with Lasso',
                                                           'lassoinvert':'6:Invert with Lasso'})),
                           Item('unselectonaction',label='Clear Selection on Action?'), 
                           Item('clearsel',show_label=False),
                           Item('weightchangesel',show_label=False),
                           Item('weightchangeto',label='Weight'),
                           Item('delsel',show_label=False)
                         ),title='Selection Options')
    
    def __init__(self,xdata,ydata,mod=None,weights=None,include_models=None,exclude_models=None,**kwargs):
        self.modelpanel = View(Label('empty'),kind='subpanel',title='model editor')
        
        self.tmodel = _TraitedModel(mod)

        self.on_trait_change(self._paramsChanged,'tmodel.paramchange')
        
        self.modelselector = NewModelSelector(include_models,exclude_models)
        
        self.data = [xdata,ydata]
        
        if weights is None:
            self.weights = np.ones_like(xdata)
            self.weighttype = 'equal'
        else:
            self.weights = np.array(weights,copy=True)
            self.savews = True
        
        pd = ArrayPlotData(xdata=self.data[0],ydata=self.data[1],weights=self.weights)
        self.plot = plot = Plot(pd,resizable='hv')
        
        self.scatter = plot.plot(('xdata','ydata','weights'),name='data',
                         color_mapper=_cmapblack if self.weights0rem else _cmap,
                         type='cmap_scatter', marker='circle')[0]
                        
        if not isinstance(mod,Model):
            self.fitmodel = True
            
        self.updatemodelplot = False #force plot update - generates xmod and ymod
        plot.plot(('xmod','ymod'),name='model',type='line',line_style='dash',color='black',line_width=2)
        
        self.on_trait_change(self._rangeChanged,'plot.index_mapper.range.updated')
        
        plot.tools.append(PanTool(plot,drag_button='right'))
        plot.overlays.append(ZoomTool(plot))
#        self.filloverlay = WeightFillOverlay(plot)
#        if self.weights0rem:
#            self.plot.overlays.append(self.filloverlay)
        
        
        self.scattertool = None
        
        #scatter.tools.append(ScatterInspector(scatter))
        self.scatter.overlays.append(ScatterInspectorOverlay(self.scatter, 
                        hover_color = "black",
                        selection_color="black",
                        selection_outline_color="red",
                        selection_line_width=2))
                        
#        self.ls = lasso_selection = LassoSelection(component=scatter,selection_datasource=scatter.index)
#        scatter.active_tool = lasso_selection
#        lasso_overlay = LassoOverlay(lasso_selection=lasso_selection,
#                                     component=scatter)
#        scatter.overlays.append(lasso_overlay)
        
        self.colorbar = colorbar = ColorBar(index_mapper=LinearMapper(range=plot.color_mapper.range),
                                            color_mapper=plot.color_mapper.range,
                                            plot=plot,
                                            orientation='v',
                                            resizable='v',
                                            width = 30,
                                            padding = 5)
        colorbar.padding_top = plot.padding_top
        colorbar.padding_bottom = plot.padding_bottom
        colorbar._axis.title = 'Weights'
        
        self.plotcontainer = container = HPlotContainer(use_backbuffer=True)
        container.add(plot)
        container.add(colorbar)
        
        super(FitGui,self).__init__(**kwargs)
        
    def _weights0rem_changed(self,old,new):
        if new:
            self.plot.color_mapper = _cmapblack(self.plot.color_mapper.range)
        else:
            self.plot.color_mapper = _cmap(self.plot.color_mapper.range)
        self.plot.request_redraw()
#        if old and self.filloverlay in self.plot.overlays:
#            self.plot.overlays.remove(self.filloverlay)
#        if new:
#            self.plot.overlays.append(self.filloverlay)
#        self.plot.request_redraw()
        
    def _paramsChanged(self):
        self.updatemodelplot = True
            
    def _nmod_changed(self):
        self.updatemodelplot = True
    
    def _rangeChanged(self):
        self.updatemodelplot = True
    
    def _updatemodelplot_fired(self,new):
        #if False (e.g. button click), update regardless, otherwise check for autoupdate
        if new and not self.autoupdate:
            return

        mod = self.tmodel.model
        if mod:
            #xd = self.data[0]
            #xmod = np.linspace(np.min(xd),np.max(xd),self.nmod)
            xl = self.plot.index_range.low
            xh = self.plot.index_range.high
            xmod = np.linspace(xl,xh,self.nmod)
            ymod = self.tmodel.model(xmod)
            
            self.plot.data.set_data('xmod',xmod)
            self.plot.data.set_data('ymod',ymod)

        else:
            self.plot.data.set_data('xmod',[])
            self.plot.data.set_data('ymod',[])
    
    def _fitmodel_fired(self):
        preaup = self.autoupdate
        try:
            self.autoupdate = False
            xd,yd = self.data
            kwd = {'x':xd,'y':yd}
            if xd.shape == self.weights.shape:
                w = self.weights
                if self.weights0rem:
                    m = w!=0
                    w = w[m]
                    kwd['x'] = kwd['x'][m]
                    kwd['y'] = kwd['y'][m]
                kwd['weights'] = w
            self.tmodel.fitdata = kwd
        finally:
            self.autoupdate = preaup
            
        self.updatemodelplot = True
        
    def _newmodel_fired(self,newval):
        if isinstance(newval,basestring):
            self.tmodel = _TraitedModel(newval)
        else:
            if self.modelselector.edit_traits(kind='modal').result:
                cls = self.modelselector.selectedmodelclass
                if cls is None:
                    self.tmodel = _TraitedModel(None)
                elif cls._args is None:
                    self.tmodel = _TraitedModel(cls(self.modelselector.modelargnum))
                else:
                    self.tmodel = _TraitedModel(cls())
            else: #cancelled
                return      
            
        #TODO:generalize this to a model auto-setting method?
        from ..models import SmoothSplineModel
        if isinstance(self.tmodel.model,SmoothSplineModel):
                self.tmodel.model.s = len(self.data[0])
        self.fitmodel = True
                
    def _get_nomodel(self):
        return self.tmodel.model is None
    
    def _get_weightsvary(self):
        w = self.weights
        return np.any(w!=w[0])if len(w)>0 else False
    
    def _get_plotname(self):
        xlabel = self.plot.x_axis.title
        ylabel = self.plot.y_axis.title
        if xlabel == '' and ylabel == '':
            return ''
        else:
            return xlabel+' vs '+ylabel
    def _set_plotname(self,val):
        if isinstance(val,basestring):
            val = val.split('vs')
            if len(val) ==1:
                val = val.split('-')
            val = [v.strip() for v in val]
        self.x_axis.title = val[0]
        self.y_axis.title = val[1]
    
    
    #selection-related
    def _scattertool_changed(self,old,new):
        if old is not None and 'lasso' in old:
            if new is not None and 'lasso' in new:
                #connect correct callbacks
                self.lassomode = new.replace('lasso','')
                return
            else:
                #TODO:test
                self.scatter.tools[-1].on_trait_change(self._lasso_handler,
                                            'selection_changed',remove=True) 
                del self.scatter.overlays[-1]
                del self.lassomode
        elif old == 'clickimmediate':
            self.scatter.index.on_trait_change(self._immediate_handler,
                                            'metadata_changed',remove=True)        
                
        self.scatter.tools = []    
        if new is None:
            pass
        elif 'click' in new:
            smodemap = {'clickimmediate':'single','clicksingle':'single',
                        'clicktoggle':'toggle'}
            self.scatter.tools.append(ScatterInspector(self.scatter,
                                      selection_mode=smodemap[new]))
            if new == 'clickimmediate':
                self.clearsel = True
                self.scatter.index.on_trait_change(self._immediate_handler,
                                                    'metadata_changed')
        elif 'lasso' in new:
            lasso_selection = LassoSelection(component=self.scatter,
                                    selection_datasource=self.scatter.index)
            self.scatter.tools.append(lasso_selection)
            lasso_overlay = LassoOverlay(lasso_selection=lasso_selection,
                                         component=self.scatter)
            self.scatter.overlays.append(lasso_overlay)
            self.lassomode = new.replace('lasso','')
            lasso_selection.on_trait_change(self._lasso_handler,
                                            'selection_changed')
        else:
            raise TraitsError('invalid scattertool value')
        
    def _weightchangesel_fired(self):
        self.weights[self.selectedi] = self.weightchangeto
        if self.unselectonaction:
            self.clearsel = True
            
        self._sel_alter_weights()
        self.lastselaction = 'weightchangesel'
    
    def _delsel_fired(self):
        self.weights[self.selectedi] = 0
        if self.unselectonaction:
            self.clearsel = True
        
        self._sel_alter_weights()
        self.lastselaction = 'delsel'
        
    def _sel_alter_weights(self):
        if self.weighttype != 'custom':
            self._customweights = self.weights
            self.weighttype = 'custom'
        self.weightsChanged()
            
    def _clearsel_fired(self,event):
        if isinstance(event,list):
            self.scatter.index.metadata['selections'] = event
        else:
            self.scatter.index.metadata['selections'] = list()
        
    def _lasso_handler(self):
        lassomask = self.scatter.index.metadata['selection'].astype(int)
        clickmask = np.zeros_like(lassomask)
        clickmask[self.scatter.index.metadata['selections']] = 1
        
        if self.lassomode == 'add':
            mask = clickmask | lassomask
        elif self.lassomode == 'remove':
            mask = clickmask & ~lassomask
        elif self.lassomode == 'invert':
            mask = np.logical_xor(clickmask,lassomask)
        else:
            raise TraitsError('lassomode is in invalid state')
        
        self.scatter.index.metadata['selections'] = list(np.where(mask)[0])
        
    def _immediate_handler(self):
        sel = self.selectedi
        if len(sel) > 1:
            self.clearsel = True
            raise TraitsError('selection error in immediate mode - more than 1 selection')
        elif len(sel)==1:
            if self.lastselaction != 'None':
                setattr(self,self.lastselaction,True)
            del sel[0]
            
    def _savews_fired(self):
        self._savedws = self.weights.copy()
    
    def _loadws_fired(self):
        self.weights = self._savedws
            
    def _get_selectedi(self):
        return self.scatter.index.metadata['selections']
    
    @on_trait_change('data',post_init=True)
    def dataChanged(self):        
        pd = self.plot.data
        #TODO:make set_data apply to both simultaneously?
        pd.set_data('xdata',self.data[0])
        pd.set_data('ydata',self.data[1])
        
        self.updatemodelplot = True
        
    @on_trait_change('weights',post_init=True)    
    def weightsChanged(self):
        self.plot.data.set_data('weights',self.weights)
        self.plot.plots['data'][0].color_mapper.range.refresh()
        
        if self.weightsvary:
            if self.colorbar not in self.plotcontainer.components:
                self.plotcontainer.add(self.colorbar)
                self.plotcontainer.request_redraw()
        elif self.colorbar in self.plotcontainer.components:
                self.plotcontainer.remove(self.colorbar)
                self.plotcontainer.request_redraw()
            
    
    def _weighttype_changed(self, name, old, new):
        if old == 'custom':
            self._customweights = self.weights
        
        if new == 'custom':
            self.weights = self._customweights #if hasattr(self,'_customweights') else np.ones_like(self.data[0])
        elif new == 'equal':
            self.weights = np.ones_like(self.data[0])
        elif new == 'lin bins':
            self.weights = binned_weights(self.data[0],10,False)
        elif new == 'log bins':
            self.weights = binned_weights(self.data[0],10,True)
        else:
            raise TraitError('Invalid Enum value on weighttype')
        
    def getModelInitStr(self):
        """
        returns a python code string that can be used to generate a 
        model with parameters matching that of this FitGui
        """
        mod = self.tmodel.model
        if mod is None:
            return 'None'
        else:
            parstrs = []
            for p,v in mod.pardict.iteritems():
                parstrs.append(p+'='+str(v))
            if mod.__class__._args is None: #varargs need to have the first argument give the right number
                varcount = len(mod.params)-len(mod.__class__._statargs)
                parstrs.insert(0,str(varcount))
            return '%s(%s)'%(mod.__class__.__name__,','.join(parstrs))
    
    def getModelObject(self):
        return self.tmodel.model
            
            
def fit_data(xdata,ydata,model=None,**kwargs):
    """
    fit a 2d data set using the FitGui interface.  Returns the model or None
    if fitting is cancelled.
    """
    
    fg = FitGui(xdata,ydata,model,**kwargs)
    if model is not None and not isinstance(model,Model):
        fg.fitmodel = True
    res = fg.configure_traits(kind='livemodal')
    
    if res:
        return fg.getModelObject()
    else:
        return None
    
class MultiFitGui(HasTraits):
    """
    data should be c x N where c is the number of data columns/axes and N is 
    the number of points
    """
    doplot3d = Bool(False)
    replot3d = Button('Replot 3D')
    scalefactor3d = Float(0)
    do3dscale = Bool(False)
    nmodel3d = Int(1024)
    usecolor3d = Bool(False)
    color3d = Color((0,0,0))
    scene3d = Instance(MlabSceneModel,())
    plot3daxes = Tuple(('x','y','z'))
    data = Array(shape=(None,None))
    weights = Array(shape=(None,))
    curveaxes = List(Tuple(Int,Int))
    axisnames = Dict(Int,Str)
    invaxisnames = Property(Dict,depends_on='axisnames')
    
    fgs = List(Instance(FitGui))
    
    
    traits_view = View(HGroup(VGroup(Item('doplot3d',label='3D Plot?'),
                              Item('scene3d',editor=SceneEditor(scene_class=MayaviScene),show_label=False,resizable=True,visible_when='doplot3d'),
                              Item('plot3daxes',editor=TupleEditor(cols=3,labels=['x','y','z']),label='Axes',visible_when='doplot3d'),
                              HGroup(Item('do3dscale',label='Scale by weight?',visible_when='doplot3d'),
                              Item('scalefactor3d',label='Point scale',visible_when='doplot3d'),
                              Item('nmodel3d',label='Nmodel',visible_when='doplot3d')),
                              HGroup(Item('usecolor3d',label='Use color?',visible_when='doplot3d'),Item('color3d',label='Relation Color',visible_when='doplot3d',enabled_when='usecolor3d')),
                              Item('replot3d',show_label=False,visible_when='doplot3d'),
                              ),
                              Item('fgs',editor=ListEditor(use_notebook=True,page_name='.plotname'),style='custom',show_label=False)),
                              resizable=True,width=1280,height=800,buttons=['OK','Cancel'],title='Multiple Model Data Fitters')
    
    def __init__(self,data,names=None,models=None,weights=None,**traits):
        super(MultiFitGui,self).__init__(**traits)
        self._lastcurveaxes = None
        #self.scene3d = None
        
        data = np.array(data,copy=False)
        if weights is None:
            self.weights = np.ones(data.shape[1])
        else:
            self.weights = np.array(weights)
        
        self.data = data
        if data.shape[0] < 2:
            raise ValueError('Must have at least 2 columns')
        
        
        if isinstance(names,basestring):
            names = names.split(',')
        if names is None:
            if len(data) == 2:
                self.axisnames = {0:'x',1:'y'}
            elif len(data) == 3:
                self.axisnames = {0:'x',1:'y',2:'z'}
            else:
                self.axisnames = dict((i,str(i)) for i in data)
        elif len(names) == len(data):
            self.axisnames = dict([t for t in enumerate(names)])
        else:
            raise ValueError("names don't match data")
        
        #default to using 0th axis as parametric
        self.curveaxes = [(0,i) for i in range(len(data))[1:]]
        
        if models is not None:
            if len(models) != len(data)-1:
                raise ValueError("models don't match data")
            for i,m in enumerate(models):
                fg = self.fgs[i]
                fg.tmodel = _TraitedModel(m)
                if not isinstance(m,Model):
                    fg.fitmodel = True
        
        
        
    def _data_changed(self):
        self.curveaxes = [(0,i) for i in range(len(self.data))[1:]]

    def _axisnames_changed(self):
        for ax,fg in zip(self.curveaxes,self.fgs):
            fg.plot.x_axis.title = self.axisnames[ax[0]] if ax[0] in self.axisnames else ''
            fg.plot.y_axis.title = self.axisnames[ax[1]] if ax[1] in self.axisnames else ''
        self.plot3daxes = (self.axisnames[0],self.axisnames[1],self.axisnames[2] if len(self.axisnames) > 2 else self.axisnames[1])
        
    @on_trait_change('curveaxes[]')
    def _curveaxes_update(self,names,old,new):
        ax=[]
        for t in self.curveaxes:
            ax.append(t[0])
            ax.append(t[1])
        if set(ax) != set(range(len(self.data))):
            self.curveaxes = self._lastcurveaxes
            return #TOOD:check for recursion
            
        if self._lastcurveaxes is None:
            self.fgs = [FitGui(self.data[t[0]],self.data[t[1]],weights=self.weights) for t in self.curveaxes]
            for ax,fg in zip(self.curveaxes,self.fgs):
                fg.plot.x_axis.title = self.axisnames[ax[0]] if ax[0] in self.axisnames else ''
                fg.plot.y_axis.title = self.axisnames[ax[1]] if ax[1] in self.axisnames else ''
        else:
            for i,t in enumerate(self.curveaxes):
                if  self._lastcurveaxes[i] != t:
                    self.fgs[i] = fg = FitGui(self.data[t[0]],self.data[t[1]],weights=self.weights)
                    ax = self.curveaxes[i]
                    fg.plot.x_axis.title = self.axisnames[ax[0]] if ax[0] in self.axisnames else ''
                    fg.plot.y_axis.title = self.axisnames[ax[1]] if ax[1] in self.axisnames else ''
        
        self._lastcurveaxes = self.curveaxes
        
    def _doplot3d_changed(self,new):
        if new:
            self.replot3d = True
    def _plot3daxes_changed(self):
        self.replot3d = True
        
    @on_trait_change('weights',post_init=True)
    def weightsChanged(self):
        for fg in self.fgs:
            if fg.weighttype != 'custom':
                fg.weighttype = 'custom'
            fg.weights = self.weights
            
        
    @on_trait_change('data','fgs','replot3d','weights')
    def _do_3d(self):
        if self.doplot3d:
            M = self.scene3d.mlab
            try:
                xi = self.invaxisnames[self.plot3daxes[0]]
                yi = self.invaxisnames[self.plot3daxes[1]]
                zi = self.invaxisnames[self.plot3daxes[2]]
                
                x,y,z = self.data[xi],self.data[yi],self.data[zi]
                w = self.weights

                M.clf()
                if self.scalefactor3d == 0:
                    sf = x.max()-x.min()
                    sf *= y.max()-y.min()
                    sf *= z.max()-z.min()
                    sf = sf/len(x)/5
                    self.scalefactor3d = sf
                else:
                    sf = self.scalefactor3d
                glyph = M.points3d(x,y,z,w,scale_factor=sf)
                glyph.glyph.scale_mode = 0 if self.do3dscale else 1
                M.axes(xlabel=self.plot3daxes[0],ylabel=self.plot3daxes[1],zlabel=self.plot3daxes[2])
                
                try:
                    xs = np.linspace(np.min(x),np.max(x),self.nmodel3d)
                    
                    #find sequence of models to go from x to y and z
                    ymods,zmods = [],[]
                    for curri,mods in zip((yi,zi),(ymods,zmods)):
                        while curri != xi:
                            for i,(i1,i2) in enumerate(self.curveaxes):
                                if curri==i2:
                                    curri = i1
                                    mods.insert(0,self.fgs[i].tmodel.model)
                                    break
                            else:
                                raise KeyError
                    
                    ys = xs
                    for m in ymods:
                        ys = m(ys)
                    zs = xs
                    for m in zmods:
                        zs = m(zs)
                    
                    if self.usecolor3d:
                        c = (self.color3d[0]/255,self.color3d[1]/255,self.color3d[2]/255)
                        M.plot3d(xs,ys,zs,color=c)
                    else:
                        M.plot3d(xs,ys,zs,np.arange(len(xs)))
                except (KeyError,TypeError):
                    M.text(0.5,0.75,'Underivable relation')
            except KeyError:
                M.clf()
                M.text(0.25,0.25,'Data problem')
                
            
            
    @cached_property        
    def _get_invaxisnames(self):
        d={}
        for k,v in self.axisnames.iteritems():
            d[v] = k
        return d
    
def fit_data_multi(data,names=None,weights=None,models=None):
    """
    fit a data set consisting of a variety of curves simultaneously
    
    returns a tuple of models e.g. [xvsy,xvsz]
    """
    
    if len(data.shape) !=2 or data.shape[0]<2:
        raise ValueError('data must be 2D with first dimension >=2')
    
    if models is not None and len(models) != data.shape[0]:
        raise ValueError('Number of models does not match number of data sets')
    
    mfg = MultiFitGui(data,names,models,weights=weights)

    res = mfg.configure_traits(kind='livemodal')

    if res:
        return tuple([fg.tmodel.model for fg in mfg.fgs])
    else:
        return None
