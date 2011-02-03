#Copyright 2009 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
This module contains GUI applications, elements, and adapters for various 
parts of astropysics.

Note that this module makes heavy use of Enthought Traits and TraitsGUI 
(http://code.enthought.com/projects/traits/) -- if it is not installed, 
nearly everything in this package will raise an exception.

GUI tools (each in their own modules):

* spectarget: a tool for making slitmasks from fiducial CMDs and photometric
  catalogs.  Currently only has direct support for Keck/DEIMOS
* spylot: spectrum plotter

Convinence functions for launching the apps and the primary application classes
are also aliased in this module for ease of use. This also includes the
:func:`pymodelfit.fitgui.fit_data` function, used for launching the FitGUI
interactive curve-fitter.
 
"""
try:
    from spectarget import spec_target,SpecTarget
    from pymodelfit.fitgui import fit_data
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

