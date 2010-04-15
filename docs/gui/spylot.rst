
Spylot -- Spectrum Plotter
==========================

This application is a spectrum plotting/interactive analysis tool based on `Traits <http://code.enthought.com/projects/traits/>`_.  It is essentially a GUI wrapper and interface around a collection of :class:`astropysics.spec.Spectrum` objects.

It can be run as part of a python script or interactively in ipython by generating :class:`astropysics.gui.spylot.Spylot` objects or calling the :func:`astropysics.gui.spylot.spylot_specs` function as detailed below.  A command-line script 'spylot' is also installed when you install astropysics.  The command-line tool takes the name of the spectrum file to display and the following options:  

-h, --help                             show help message and exit
-t TYPE, --type=TYPE                   file format of spectra: "wcs", "deimos", or  "astropysics"(default)
-e EXT, --extension=EXT                Fits extension number
-c CONFIGFILE, --config=CONFIGFILE     File to save Spylot configuration to


.. todo:: Write a Tutorial/examples for both code and command-line script


.. autoclass:: astropysics.gui.spylot.Spylot
   :members:
   :undoc-members: 
   :exclude-members: majorlineeditor,minorlineeditor
   
.. autofunction:: astropysics.gui.spylot.spylot_specs
