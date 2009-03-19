#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains functions and classes for specialized plotting functions

generally, 2D plots depend on matplotlib, 3D plots depend on MayaVI
"""

from __future__ import division
import numpy as np

try:
    import enthought.mayavi
except ImportError:
    from warnings import warn
    warn('MayaVI 2 not found -3D plotting probably will not work')
    
try:
    import matplotlib
except ImportError:
    from warnings import warn
    warn('Matplotlib not found -most 2D plotting probably will not work')
    