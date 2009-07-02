#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the FitGui gui.
"""
from __future__ import division,with_statement
import numpy as np

from enthought.traits.api import HasTraits,Instance,Int,Float,Bool,Button,Event, \
                                 Property,on_trait_change,Array,List,Str
from enthought.traits.ui.api import View,Item,Label,Group,VGroup,HGroup, \
                                    InstanceEditor,EnumEditor
from enthought.traits.ui.menu import ModalButtons
from enthought.chaco.api import Plot,ArrayPlotData,jet
from enthought.chaco.tools.api import PanTool, ZoomTool,SelectTool
from enthought.enable.component_editor import ComponentEditor


from ..models import list_models,get_model,Model

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
        
        
        if model is None:
            g = Group()
            g.content.append(Label('No Model Selected'))
        else:
            g = Group(columns=8,label=self.modelname,show_border=True)
            for p in model.params:
                self.add_trait(p,Float)
                setattr(self,p,getattr(model,p))
                self.on_trait_change(self._param_change_handler,p)
                g.content.append(Item(p))
                
                ffp = 'fixfit_'+p
                self.add_trait(ffp,Bool)
                setattr(self,ffp,False)
                self.on_trait_change(self._param_change_handler,ffp)
                g.content.append(Item(ffp,label='Fix?'))
        
        view = View(g,buttons=['Apply','Revert','OK','Cancel'])
        self.trait_view('traits_view',view)
    
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
        
class _NewModelSelector(HasTraits):
    modelnames = List
    selectedname = Str('No Model')
    selectedmodelclass = Property
    
    traits_view = View(Item('selectedname',editor=EnumEditor(name='modelnames')),
                       buttons=['OK','Cancel'])
    
    def __init__(self,**traits):
        super(_NewModelSelector,self).__init__(**traits)
        
        self.modelnames = list_models()
        self.modelnames.insert(0,'No Model')
        
    def _get_selectedmodelclass(self):
        n = self.selectedname
        if n == 'No Model':
            return None
        else:
            return get_model(n)
    
class FitGui(HasTraits):
    plot = Plot
    tmodel = Instance(_TraitedModel,allow_none = False)
    nomodel = Property
    newmodel = Button('New Model...')
    fitmodel = Button('Fit Model')
    updatemodelplot = Button('Update Model Plot')
    autoupdate = Bool(True)
    data = Array(dtype=float,shape=(2,None))
    weights = Array
    
    nmod = Int(100)
    #modelpanel = View(Label('empty'),kind='subpanel',title='model editor')
    modelpanel = View
    
    traits_view = View(VGroup(
                       Item('plot', editor=ComponentEditor(),show_label=False),
                       HGroup(Item('newmodel',show_label=False),
                              Item('fitmodel',show_label=False),
                              Item('nmod',label='Nmodel'),
                              Item('updatemodelplot',show_label=False,enabled_when='not autoupdate'),
                              Item('autoupdate',label='Auto?')),
                       Item('tmodel',show_label=False,style='custom',editor=InstanceEditor(kind='subpanel'))
                      ),
                    resizable=True, title='Data Fitting',buttons=['OK','Cancel']
                    )
                    
    panel_view = View(VGroup(
                       Item('plot', editor=ComponentEditor(),show_label=False),
                       HGroup(Item('tmodel.modelname',show_label=False,style='readonly'),
                              Item('nmod',label='Number of model points'),
                              Item('updatemodelplot',show_label=False,enabled_when='not autoupdate'),
                              Item('autoupdate',label='Auto?'))
                      ),
                    title='Data Fitting'
                    )
    
    def __init__(self,xdata,ydata,mod=None,weights=None,**kwargs):
        self.modelpanel = View(Label('empty'),kind='subpanel',title='model editor')

        self.tmodel = _TraitedModel(mod)

        self.on_trait_change(self._paramchanged,'tmodel.paramchange')
        
        self.data = [xdata,ydata]
        self.weights = np.ones_like(xdata) if weights is None else weights
        
        pd = ArrayPlotData(xdata=self.data[0],ydata=self.data[1],weights=self.weights)
        self.plot = plot = Plot(pd,resizable='hv')
        
        plot.plot(('xdata','ydata','weights'),name='data',type='cmap_scatter',color_mapper=jet,marker='circle')
        if not isinstance(mod,Model):
            self.fitmodel = True
        self.updatemodelplot = False #force plot update - generates xmod and ymod
        plot.plot(('xmod','ymod'),name='model',type='line',line_style='dash',color='black',line_width=2)
        
        plot.tools.append(PanTool(plot))
        plot.tools.append(ZoomTool(plot))
        
        super(FitGui,self).__init__(**kwargs)
        
    def _paramchanged(self,new):
        self.updatemodelplot = True
            
    def _nmod_changed(self):
        self.updatemodelplot = True
        
    def _updatemodelplot_fired(self,new):
        #if False (e.g. button click), update regardless, otherwise check for autoupdate
        if new and not self.autoupdate:
            return

        mod = self.tmodel.model
        if mod:
            xd = self.data[0]
            xmod = np.linspace(np.min(xd),np.max(xd))
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
                kwd['weights'] = self.weights
            self.tmodel.fitdata = kwd
        finally:
            self.autoupdate = preaup
            
        self.updatemodelplot = True
        
    def _newmodel_fired(self):
        msel = _NewModelSelector()
        if msel.edit_traits(kind='modal').result:
            cls = msel.selectedmodelclass
            if cls is None:
                self.tmodel = _TraitedModel(None)
            else:
                self.tmodel = _TraitedModel(cls())
            self.fitmodel = True
                
                
    def _get_nomodel(self):
        return self.tmodel.model is None
    
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
        
        
        
    def getModelInitStr(self):
        """
        returns a string that can be used to generate a model like the one in 
        this FitGui
        """
        mod = self.tmodel.model
        if mod is None:
            return 'None'
        else:
            parstrs = []
            for p,v in mod.pardict:
                parstrs.append(p+'='+str(v))
            return '%s(%s)'%(mod.__class__.__name__,','.join(parstrs))
    
    def getModelObject(self):
        return self.tmodel.model
            
def fit1d(xdata,ydata,model=None,**kwargs):
    """
    fit a 1d data set using the FitGui interface.  Returns the model or None
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