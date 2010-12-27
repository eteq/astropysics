Getting Started/Tutorial
========================

If you have not done so, install astropysics as described in :doc:`install`, and
be sure to run the setup command ``astpys-setup`` as described in
:ref:`inst-setup`. Be sure you have IPython installed for the rest of this
section to function correctly.

Interactive Environment
-----------------------

Astropysics uses `IPython <http://ipython.scipy.org/>`_ to provide an
interactive environment to run python. To start using astropysics in ipython,
just run the helper script ``ipyastpys`` - that will run ipython with a
customized profile that automatically imports commonly-used parts of of
astropysics (and :mod:`numpy`).  


Projects
--------

The ``ipyastpys`` environment also supports "projects", allowing the interactive
environment to be started such that it automatically moves to a given directory
and runs a given script.  This directory can then hold all necessary data files
and the script can load data and store functions to generate plots for the 
project/paper.  A project can be created using the command::

    ipyastpys -c projectname 
    
This will create a project named "projectname", with a directory
"/path/to/currentdir/projectname" and script "projectname/projectname.py". If
the directory or script don't exist, they will be created. A different directory or
script name can be used by calling the script as::

    ipyastpys -c projectname /path/to/projectdir projectscript.py
    
Note that the script must be inside the project directory (in the above example,
the script's full path is "/path/to/projectdir/projectscript.py").

The interactive environment should then be started using::

    ipyastpys -p projectname
    
And it will start ipython in /path/to/projectdir, with the projectscript.py
script automatically run inside the interactive environment (the ``-s`` option
can be used to give projectscript.py any necessary command line arguments).


If a project script does not already exist, it will be generated from
:download:`a template </../astropysics/data/project_template.py>` that allows
for easy interactive work to create plots. In the script, add as many functions
as necessary with the function decorator :func:`plotfunc` (see the template
itself for detailed syntax), and then load any necessary data variables at
the end of the script. The :func:`make_plots` function can then be used to 
generate the plots by passing it the name of the plot function.  The script can
also be used at the command line to auto-generate all plots.


Module Tutorials
----------------

Most of the modlues in astropysics have their own tutorials and examples that
are relevant to their specific functionality. This is the best place to go to
find examples of how to actually use astropysics. See :doc:`coremods/intro` and
:doc:`gui/intro` for a list of these modules.

    .. todo:: Add links to the relevant tutorials in some automated way
