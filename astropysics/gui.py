#©2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains GUI elements and adapters for various parts of astropysics.

Note that this module makes heavy use of Enthought Traits and TraitsGUI 
(http://code.enthought.com/projects/traits/) -- if it is not installed, nearly 
everything in this package will raise an exception
"""

from __future__ import division
import numpy as np

try:
    import enthought.traits.api as tapi
    import enthought.traits.ui.api as uapi
except ImportError:
    from warnings import warn
    warn('Traits or Traits GUI not installed! most of gui module will be worthless')