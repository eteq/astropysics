#Copyright 2010 Erik Tollerud
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
======================================
utils -- utility classes and functions
======================================

Overview
--------

The :mod:`utils` module contains classes and funtions of general utility used in
multiple places throughout `astropysics`. Some of these are
astropyhysics-specific algorithms while others are more python tricks.

The :mod:`~astropysics.utils` module is composed of three submodules to make
organization clearer. The submodules are fairly different from each other, but
the main uniting theme is that all of these submodules are not specific to a
particularly astronomical application. Hence, they are available for re-use in
astropysics or anywhere else they are deemed useful. The detailed documentation
for each of the sub-modules is described below.

.. note::

    The :mod:`~astropysics.utils` module is composed of submodules that can be
    accessed seperately. However, they are all also included in the base module.
    Thus, as an example, ``astropysics.utils.gen.DataObjectRegistry`` and
    ``astropysics.utils.DataObjectRegistry`` are different names for the same
    class (:class:`~astropysics.utils.gen.DataObjectRegistry`). The
    ``astropysics.utils.DataObjectRegistry`` usage is preferred as this allows
    the internal organization to be changed if it is deemed necessary.

    
    
utils.gen -- general-purpose (python-specific) utilities
--------------------------------------------------------

.. automodule:: astropysics.utils.gen
   :members:
   :undoc-members:
   :show-inheritance:


utils.alg -- common/repeated algorithms
---------------------------------------

.. automodule:: astropysics.utils.alg
   :members:
   :undoc-members:
   :show-inheritance:
   
utils.stats -- statistcs-related tools
--------------------------------------

.. automodule:: astropysics.utils.stats
   :members:
   :undoc-members:
   :show-inheritance:
   
utils.io -- input/output and various file format tools
------------------------------------------------------

.. automodule:: astropysics.utils.io
   :members:
   :undoc-members:
   :show-inheritance:

"""


from gen import *
from alg import *
from io import *
from stats import *

del np #clean up namespace
