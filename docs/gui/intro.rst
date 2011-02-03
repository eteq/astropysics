GUI applications
================

Astropysics includes a variety of graphical applications for various data
analysis tasks. Their full power is acheived when used interactively or as part
of scripts, but some operate as stand-alone command-line tools when the
situation warrants.

Note that these applications make heavy use of `Enthought Traits
<http://code.enthought.com/projects/traits/>`_ and the associated `Enthought
Tools Suite <http://code.enthought.com/projects/index.php>`_ if they are not
installed, most of these will not function.
 
Another important related GUI used in astropysics is the Fit GUI from 
:mod:`pymodelfit`.  This GUI is used wherever interactive curve-fitting is
needed (and was, in fact, originally made for astropysics).

.. toctree::
   :maxdepth: 2
   
   spylot
   spectarget
   
   
the :mod:`gui` module also imports the convinience functions and primary
application classes for these GUIs. These include:

* :func:`fit_data <astropysics.gui.fitgui.fit_data>` (imported from :mod:`pymodelfit`)
* :class:`Spylot <astropysics.gui.spylot.Spylot>`
* :func:`spylot_specs <astropysics.gui.spylot.spylot_specs>`
* :class:`SpecTarget <astropysics.gui.spectarget.SpecTarget>`
* :func:`spec_target <astropysics.gui.spectarget.spec_target>`
