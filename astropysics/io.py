#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains functions and classes for various forms of data/file
input and output As well as library retrieval for built-in data.
"""

from __future__ import division
import numpy as np

try:
    import pyfits
except ImportError:
    from warnings import warn
    warn('pyfits not found - all FITS-related IO will not work')