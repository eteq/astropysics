#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""

======
models
======

The :mod:`models` module contains objects and functions for fitting data to
models as well as calculations and estimates from these models. The available
fitting algorithms are those from :mod:`scipy.optimize` or the `PyMC
<http://code.google.com/p/pymc/>`_ .

The aim of these classes are mostly for easily and quickly generating a range of
models. There are classes in :mod:`core` for various dimensionalities, and all
that is necessary is to subclass one of the base classes that end in "Auto" and
implement just a function :meth:`f` describing a parameterized model. The
parameters will be automatically inferred from the :meth:`f` function and
fitting, plotting, and the :mod:`astropysics.gui.fitgui` gui will be available.

The module includes a model registry to easily access the builtin models and any
others the user wishes to implement. Given a subclass `MyModel`, the following
will register it::

    from astropysics.model import register_model
    register_model(MyModel,[name])
    
and it will subsequently be available as a model under whatever name you
provided (or a default name based on the class name with "Model" removed) in any
place where a model needs to be specified in :mod:`astropysics`.

.. todo:: more specific examples




API Overview
------------

The :mod:`astropysics.models` module is composed of two submodules that are
both imported into the main module. The first, :mod:`astropysics.models.core`
contains the classes and functions that structure and do most of the work of the
models. The second, :mod:`astropysics.models.builtins` contains a default
set of models.  There is also a :mod:`astropysics.models.pca` module for 
performing n-dimensional Principal Component Analysis and plotting the results.

models.core
-----------

.. automodule:: astropysics.models.core
   :members:
   :undoc-members:
   :show-inheritance:


models.builtins
---------------

.. automodule:: astropysics.models.builtins
   :members:
   :undoc-members:
   :show-inheritance:

models.pca
----------

.. automodule:: astropysics.models.pca
   :members:
   :undoc-members:
   :show-inheritance:

"""
from core import *
from builtins import *
from pca import Pca

del ABCMeta,abstractmethod,abstractproperty,np,pi #clean up namespace
