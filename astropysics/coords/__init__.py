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

=======================================================
coords -- coordinate classes and coordinate conversions
=======================================================


Overview
--------

The :mod:`coords` module contains classes and functions specifying locations of
objects on the sky as well as coordinate converstions and distance computations,
including cosmological distance and redshift calculations.

Some of the calculations involved make use of the currently selected cosmology
(see :mod:`astropysics.constants`) and hence may not function as expected if a
particularly strange cosmology is in use.

.. note::
    Timekeeping/conversion functions are kept in :mod:`astropysics.obstools`,
    although some are crucial to coordinate systems.

.. seealso::

    `Kapteyn libraries <http://www.astro.rug.nl/software/kapteyn/index.html>`_
        A set of python libraries with excellent coordinate transform and other
        capabilities.
        
    `Meeus, Jean H. "Astronomical Algorithms" ISBN 0943396352 <http://www.willbell.com/MATH/mc1.htm>`_ 
        An authoritative reference on coordinates, ephemerides, and related
        transforms in astronomy.
   
    `Standards Of Fundamental Astronomy (SOFA) <http://www.iausofa.org/>`_ 
        The IAU reference implementations for coordinates and earth rotation.  
   
    `USNO Circular 179 <http://aa.usno.navy.mil/publications/docs/Circular_179.pdf>`_ 
        An excellent description of the IAU 2000 resolutions and related
        background for defining ICRS, CIO, and related standards.
    
.. note::
    Some functionality in this module makes use of derived versions of the `SOFA
    <http://www.iausofa.org/>`_ routines. Use of these derived works requires
    inclusion of the SOFA license (included in the licenses/SOFA_LICENSE file of
    the source distribution) and notation indicating how these routines have
    been modified. In *all* cases where SOFA-derived routines are used, the
    routine's functionality is replicated exactly as the SOFA implementation
    (barring the possibility of typos), but adapted to the python language.
   
   
.. todo:: Tutorials that unify the sub-modules?



    
The :mod:`~astropysics.coords` module is composed of three submodules to make
organization clearer. The :mod:`~astropysics.coords.coordsys` module implements
classes representing useful coordinate systems in astronomy and a framework to
add additional coordinate systems as desired. It also implements standard
transformations between the various celestial and terrestrial coordinates
(although transformation to local horizontal coordinates is done with methods of
:class:`astropysics.obstools.Site`). :mod:`~astropysics.coords.ephems`
implements ephemerides for solar system objects and proper motions. Finally,
:mod:`~astropysics.coords.funcs` contains a range of utility functions including
cartesian<->spherical and other canonical transforms, as well as cosmological
distance calculations. The documentation for each of the sub-modules is
described below.

.. note::

    The :mod:`~astropysics.coords` module is composed of submodules that can be
    accessed seperately. However, they are all also included in the base module.
    Thus, as an example, ``astropysics.coords.coordsys.ICRSCoordinates`` and
    ``astropysics.coords.ICRSCoordinates`` are different names for the same
    class (:class:`~astropysics.coords.coordsys.ICRSCoordinates`). The
    ``astropysics.coords.ICRSCoordinates`` usage is preferred as this allows the
    internal organization to be changed if it is deemed necessary.

coords.coordsys -- coordinate systems and transforms
----------------------------------------------------

.. automodule:: astropysics.coords.coordsys
   :members:
   :undoc-members:
   :show-inheritance:


coords.ephems -- ephemerides and proper motions
-----------------------------------------------

.. automodule:: astropysics.coords.ephems
   :members:
   :undoc-members:
   :show-inheritance:

coords.funcs -- coordinate and distance utility functions
---------------------------------------------------------

.. automodule:: astropysics.coords.funcs
   :members:
   :undoc-members:
   :show-inheritance:

"""

#TODO: WCSlib or similar support a la Kapteyn?

from coordsys import *
from ephems import *
from funcs import *

del ABCMeta,abstractmethod,abstractproperty,np,pi #clean up namespace
