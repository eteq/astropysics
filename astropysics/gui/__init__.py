#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains GUI applications, elements, and adapters for various 
parts of astropysics.

Note that this module makes heavy use of Enthought Traits and TraitsGUI 
(http://code.enthought.com/projects/traits/) -- if it is not installed, 
nearly everything in this package will raise an exception.
"""
from __future__ import division

try:
    import enthought.traits.api as _tapi
    import enthought.traits.ui.api as _uapi
except ImportError:
    from warnings import warn
    warn('Traits or Traits GUI not installed! most of gui module will be unuseable')
    
#Use lines below as imports for new gui apps
#from __future__ import division
#from enthought.traits.api import HasTraits
#from enthought.traits.ui.api import View    