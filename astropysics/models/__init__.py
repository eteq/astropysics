"""
Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

This module contains objects and functions for fitting data to models as well as
calculations and estimates from these models.

The aim of these classes are mostly for easily and quickly generating a range of
models - subclasses with just a function "f" will do all the expected things 
right out of the box.

Currently, the main fitting algorithms are those from scipy.optimize and the 
PyMC package (http://code.google.com/p/pymc/)
"""
from modcore import *
from modbuiltins import *
from pca import Pca

del ABCMeta,abstractmethod,abstractproperty,np,pi #clean up namespace