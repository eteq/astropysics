.. Astropysics documentation master file, created by
   sphinx-quickstart on Wed Feb 24 00:07:20 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Astropysics: astrophysics utilities for python
==============================================
:Author: `Erik Tollerud <http://www.physics.uci.edu/~etolleru/>`_

Astropysics is a library containing a variety of utilities and algorithms for reducing, analyzing, and visualizing astronomical data.  Best of all, it encourages the user to leverage the existing capabilities of Python to make this quick, easy, and as painless as cutting-edge science can even actually be.  There do exist other Python packages with some of the capabilities of this project, but the goal of this project is to integrate all these tools together and make them interact in the most straightforward ways possible.  

(And to that end, if you are running one of those other projects, I'd love to help integrate our projects into a common framework!)


Installation
============

Requirements
------------
Before you touch astropysics, you'll need:

    * `Python <http://www.python.org/>`_ 2.5 or higher (2.6 highly recommended), although 3.x is not yet compatible with numpy or scipy.
    * `numpy <http://numpy.scipy.org>`_ 
    * `scipy <http://www.scipy.org/>`_
    
Follow the relevant instructions on those sites, or, far easier, install them as packages from your operating system (e.g. apt-get or synaptic on Ubuntu).  


Install
-------

Once you have the requirements satisfied, you have a few options for installing astropysics.

If you have `pip <http://pypi.python.org/pypi/pip>`_ (the new, better easy installer) or `easy_install/setuptools <http://pypi.python.org/pypi/setuptools>`_ (you should probably install pip...), just run either::

    pip install astropysics

or::

    easy_install astropysics

If you are on Ubuntu or a similar linux distribution, you may need to prefix those commands with ``sudo``.

If you are installing from source code, instead, just do::

    python setup.py install
    
where again, you may have to prefix the command with ``sudo`` if it doesn't work the first time.


If you want to be sure you have the most recent version (or want to help improve the code), you can pull the current development version from launchpad.  You must have `bazaar <http://bazaar.canonical.com/>`_ installed, and execute::

    bzr branch lp:astropysics <dirname>
    
and this will create a branch with the requested directory name containing the latest and greatest version of astropysics.  If at any time you want to update this version, go into the directory and do::

    bzr update
    
and re-install following the directions above.


Extra packages
--------------

A number of other packages are necessary for added functionality in astropysics (and all are recommended):
    * `Matplotlib <http://matplotlib.sourceforge.net/index.html>`_
        *highly recommended*, as it is necessary for all plotting (aside from the GUI applications).
    * `Traits <http://code.enthought.com/projects/traits/>`_, `TraitsGUI <http://code.enthought.com/projects/traits_gui/>`_, `Chaco <http://code.enthought.com/projects/chaco/>`_, and `Mayavi <http://code.enthought.com/projects/mayavi/>`_.  Alternatively, `ETS <http://code.enthought.com/projects/index.php>`_ is all bundled in one.
        Necessary for the interfaces in the :mod:`gui` modules.
    * `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
        Necessary for reading FITS files.
    * `vo.table <http://stsdas.stsci.edu/astrolib/vo/html/>`_
        Necessary to write VOTable files, and makes reading them much better, as well.

Documentation
=============

Astropysics is divided into two major subcomponents - the core modules that perform the calculation and contain a variety of useful classes for everyone's use, and the gui module that contains a number of useful applications.

.. toctree::
   :maxdepth: 2
   
   coremods/intro
   gui/intro
   

Bug Reports
===========

The best place to report bugs is via the `launchpad bug tracker <https://bugs.launchpad.net/astropysics>`_.  That way they won't be forgotten unless an asteroid impact destroys all of humanity (or all the launchpad servers...).


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

