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

==================================================
models -- model and data-fitting classes and tools
==================================================

The :mod:`~astropysics.models` module contains objects and functions for fitting
data to models as well as calculations and estimates from these models. The
available fitting algorithms are those from :mod:`scipy.optimize` or the `PyMC
<http://code.google.com/p/pymc/>`_ .

The aim of these classes are mostly for easily and quickly generating a range of
models. There are classes in :mod:`~astropysics.models.core` for various
dimensionalities, and helper classes to make it trivial to implement one of
these classes. Further, the module includes a model registry to easily access
the builtin models and any others the user wishes to implement. See the examples
below for details.


all that is necessary is to subclass one of the base
classes that end in "Auto" and implement just a function :meth:`f` describing a
parameterized model. The parameters will be automatically inferred from the
:meth:`f` function and fitting, plotting, and the :mod:`astropysics.gui.fitgui`
gui will be available.

Examples
--------

Accessing Models
^^^^^^^^^^^^^^^^

To access a builtin model and any custom models that have been registered, the
:func:`get_model_instance` function is the way to get a new instance of a model
using the model's name::

    >>> linmod = get_model_instance('linear')
    >>> linmod.m = 1.25
    >>> linmod.b = 3
    >>> linmod([0,1,2])
    array([ 3.  ,  4.25,  5.5 ])
    >>> linmod.__class__
    <class 'astropysics.models.builtins.LinearModel'>
    
Alternatively, a new instance can be directly instantianted from the class
itself::
    
    linmod = LinearModel(m=1.25,b=3)
    
A list of available models can be obtained with the :func:`list_models`, and
this list can be used to get only models of a specific dimensionality. The two
examples below are for scalar->scalar models, and 2D vector->scalar models)::

    list_models(baseclass=FunctionModel1D)
    list_models(baseclass=FunctionModel2DScalar)


Creating a Custom Model
^^^^^^^^^^^^^^^^^^^^^^^

While a custom model can be written by directly subclassing
:class:`FunctionModel` subclasses of various dimensionality, the much easier way
is to inherit from one of the models that ends in 'Auto'::

    from astropysics.model import FunctionModel1DAuto
    class CuberootLinearModel(FunctionModel1DAuto):
        def f(self,x,a=1,b=5,c=0):
            return a*x**(1/3) + b*x + c
        
The 'Auto' classes take care of assigning the parameter names and defaults based
on the signature of the `f` function. The resulting class thus represents a
mixed cubed-root/linear model, and all of the model-fitting, plotting, and
related tools are immediately available simply by doing::

    model = CuberootLinearModel()

Registering a Custom Model
^^^^^^^^^^^^^^^^^^^^^^^^^^

While the the `CuberootLinearModel` model from above can be used as-is with
standard python class syntax, for it to be visible in astropysics for tools like
the :mod:`fitgui` GUI interface, it must be registered with astropysics using
the :func:`register_model` function::

    from astropysics.model import register_model
    register_model(CuberootLinearModel)
    
The model will now be available in the :func:`get_model_class` and
:func:`get_model_instance` functions under the name 'cuberootlinear' (the
default name can be changed using the `name` parameter of
:func:`register_model`). Thus, the new model will be visible to anything in
astropysics that uses functional models (e.g. spectral line fitting forms,
fitgui, etc.).



Module Organization
-------------------

The :mod:`~astropysics.models` module is composed of two submodules that are
both imported into the main module. The first, :mod:`~astropysics.models.core`
contains the classes and functions that structure and do most of the work of the
models. The second, :mod:`~astropysics.models.builtins` contains a default set
of models. There is also a :mod:`~astropysics.models.pca` module for performing
n-dimensional Principal Component Analysis and plotting the results.

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
