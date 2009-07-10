#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the FitGui gui.
"""
from __future__ import division,with_statement
import numpy as np

from enthought.traits.api import HasTraits,Instance,Int,Float,Bool,Button, \
                                 Event,Property,on_trait_change,Array,List, \
                                 Tuple,Str,Dict,cached_property,Color
from enthought.traits.ui.api import View,Item,Label,Group,VGroup,HGroup, \
                                    InstanceEditor,EnumEditor,ListEditor, \
                                    TupleEditor
from enthought.traits.ui.menu import ModalButtons
from enthought.chaco.api import Plot,ArrayPlotData,jet
from enthought.chaco.tools.api import PanTool, ZoomTool,SelectTool
from enthought.enable.component_editor import ComponentEditor

#TODO:move to remove deps?
from enthought.tvtk.pyface.scene_editor import SceneEditor 
from enthought.mayavi.tools.mlab_scene_model import MlabSceneModel
from enthought.mayavi.core.ui.mayavi_scene import MayaviScene


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
        
        self.modelnames = list_models(1)
        self.modelnames.insert(0,'No Model')
        self.modelnames.remove('polynomial')
        self.modelnames.sort()
        
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
    plotname = Property
    
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
        plot.overlays.append(ZoomTool(plot))
        
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
            xmod = np.linspace(np.min(xd),np.max(xd),self.nmod)
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
        
    def _newmodel_fired(self,newval):
        if isinstance(newval,basestring):
            self.tmodel = _TraitedModel(newval)
        else:
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
            
            
def fit_curve(xdata,ydata,model=None,**kwargs):
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
    nmodel3d = Int(100)
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
                              resizable=True,width=1240,height=650,buttons=['OK','Cancel'])
    
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
    
def fit_curve_multi(data,names=None,weights=None,models=None):
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
