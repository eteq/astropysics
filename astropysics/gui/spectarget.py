#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 
"""
This package contains the internals for the SpecTarget gui.
"""
from __future__ import division
from enthought.traits.api import HasTraits
from enthought.traits.ui.api import View

class SpecTarget(HasTraits):
    def __init__(self):
        self.view = View()