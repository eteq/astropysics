
Fitgui -- Interactive Curve Fitting
===================================

This application is an interactive 1D curve-fitting tool based on `Traits <http://code.enthought.com/projects/traits/>`_.  The models for this application are based on the model framework of :mod:`astropysics.models`, using subclasses of :class:`astropysics.models.core.FunctionModel1D`.

.. todo:: Write a Tutorial/examples


.. autoclass:: astropysics.gui.fitgui.FitGui
   :members:
   :undoc-members:
   :exclude-members: modelpanel,modelselector
   
.. autofunction:: astropysics.gui.fitgui.fit_data
