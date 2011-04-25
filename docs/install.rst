Installing Astropysics
======================


Requirements
------------
Before you do anything with astropysics, you'll need:

    * `Python <http://www.python.org/>`_ 2.5 or higher (2.6 highly recommended), although 3.x is not yet supported.
    * `numpy <http://numpy.scipy.org>`_ 
    * `scipy <http://www.scipy.org/>`_
    
Follow the instructions on those sites, or, far easier, install them as packages from your operating system (e.g. apt-get or the synaptic GUI on Ubuntu, `Macports <http://www.macports.org/>`_ on OS X, etc.).  


Install
-------

Once you have the requirements satisfied, you have a few options for installing astropysics.  

.. note::
    On most unix-like systems, you will need to either execute these commands as the root user, or preface them with ``sudo``.

If you have `pip <http://pypi.python.org/pypi/pip>`_ (the new, better easy installer) or `easy_install/setuptools <http://pypi.python.org/pypi/setuptools>`_ (you should probably install pip...), just run either::

    pip install astropysics

or::

    easy_install astropysics

If you are installing from source code, instead, just do::

    python setup.py install

If you plan on using the most up-to-date version of astropysics or if you wish
to alter the source code (see :doc:`develop`), a useful way to immediately
see changes without having to re-install every time is to use the command::

    python setup.py develop
    
    
.. _inst-setup:

Setup
-----
    
After the package has been installed, at the command line, run::

    astpys-setup
    
This script does two things:

* Prompts you to select which optional packages you want to download and install.
* Configures IPython to support the ``ipyastpys`` script (described in :doc:`getstarted`).

The first of these involves an interactive process that downloads and installed
requested packages. To install all of them, type ``a`` (and hit enter),
otherwise choose a number to install that package. If you want to quit before
installing all of the package (for example, if some don't install correctly),
choose ``q``. For information on a package, type ``i#`` (where # is the number
for the package).

Note that if you can't get any packages to install, you might try running the
script as::

    astpys-setup -s
    
Depending on your operating system, you may want to use your package management
system to install the recommended packages, before running the setup (although 
you may need the more up-to-date versions given here).
    

.. _inst-rec:

Recommended Packages
--------------------

A number of other packages are necessary for added functionality in astropysics
or to provide functionality that has no need to be duplicated. These packages
can be installed with the ``astpys-setup`` script as described in
:ref:`inst-setup`, but if available from your system's package management
system, it may be better to try installing that way, instead.

    * `Matplotlib <http://matplotlib.sourceforge.net/index.html>`_
        *highly recommended*, as it is necessary for all plotting (aside from 
        the GUI applications).
            
    * `IPython <http://ipython.scipy.org/>`_
        *highly recommended*, as it is a far better interactive shell than the
        python default and has many wonderful features. Necessary for the
        ``ipyastpys`` script.
        
    * `NetworkX <http://networkx.lanl.gov/>`_
        *highly recommended*, as it is used for a variety of internal purposes
        as well as any place where a network/graph is plotted.

    * `PyGraphviz <http://networkx.lanl.gov/pygraphviz/>`_
        It might also be useful to have a closely related package for generating
        `graphviz <http://www.graphviz.org/>`_ graphs from networkx.
    
    * `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
        *highly recommended*, necessary for reading FITS files (the most common 
        astronomy data format).
        
    * `asciitable <http://cxc.cfa.harvard.edu/contrib/asciitable/>`
        A valuable tool for loading and writing ASCII tables.
        
    * `ATpy <http://atpy.github.com/>`
        Astronomical Tables in Python - a general tool for dealing with tabular
        data, both ASCII (uses asciitable) and other formats.
        
    * `pidly <http://astronomy.sussex.ac.uk/~anthonys/pidly/>`
        IDL within Python.  For those times when someone sends you an IDL code 
        that you don't have the time to convert to python, but want to be able
        to call from inside python.  Requires IDL to be installed.
        
            
    * `Traits <http://code.enthought.com/projects/traits/>`_, `TraitsGUI <http://code.enthought.com/projects/traits_gui/>`_, `Chaco <http://code.enthought.com/projects/chaco/>`_, and `Mayavi <http://code.enthought.com/projects/mayavi/>`_.  Alternatively, `ETS <http://code.enthought.com/projects/index.php>`_ is all bundled in one.
        Necessary for the interfaces in the :mod:`gui` modules::
        
            pip install ETS
            
        or::
        
            pip install traits
            pip install traitsGUI
            pip install chaco
            pip install mayavi
        
    
Astropysics also includes pythonic wrappers around some astronomy-related tools that need to be installed seperately if their functionality is desired:

    * `SExtractor <http://www.astromatic.net/software/sextractor>`_
    * `Kcorrect <http://howdy.physics.nyu.edu/index.php/Kcorrect>`_