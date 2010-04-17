#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This module contains GUI applications, elements, and adapters for various 
parts of astropysics.

Note that this module makes heavy use of Enthought Traits and TraitsGUI 
(http://code.enthought.com/projects/traits/) -- if it is not installed, 
nearly everything in this package will raise an exception.

GUI tools (each in their own modules):

* spectarget: a tool for making slitmasks from fiducial CMDs and photometric
  catalogs.  Currently only has direct support for Keck/DEIMOS
* fitgui: tools for interactively fitting ``astropysics.models`` models 
  to data.
* spylot: spectrum plotter

Convinence functions for launching the apps and the primary application classes
are also aliased in this module for ease of use.
 
"""
try:
    from spectarget import spec_target,SpecTarget
    from fitgui import fit_data,FitGui
    from spylot import spylot_specs,Spylot
except ImportError,e:
    if 'traits' in e.message:
        from warnings import warn
        warn('Traits or Traits GUI not installed! most of gui module will be unuseable')
    else:
        raise
        
    
#Use lines below as imports for new gui apps
#from __future__ import division
#from enthought.traits.api import HasTraits
#from enthought.traits.ui.api import View    

