.. Astropysics documentation master file, created by
   sphinx-quickstart on Wed Feb 24 00:07:20 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Astropysics: astrophysics utilities for python
==============================================
:Author: `Erik Tollerud <http://www.physics.uci.edu/~etolleru/>`_

Astropysics is a library containing a variety of utilities and algorithms for reducing, analyzing, and visualizing astronomical data.  Best of all, it encourages the user to leverage the existing capabilities of Python to make this quick, easy, and as painless as cutting-edge science can even actually be.  There do exist other Python packages with some of the capabilities of this project, but the goal of this project is to integrate all these tools together and make them interact in the most straightforward ways possible.  

(And to that end, if you are running one of those other projects, I'd love to help integrate our projects into a common framework!)

Documentation
=============

Astropysics is divided into two major subcomponents - the core modules that perform the calculation and contain a variety of useful classes for everyone's use, and the gui module that contains a number of useful applications.

.. toctree::
   :maxdepth: 2
   
   coremods/intro
   gui/intro


Installation
============

Requirements
------------
Before you touch astropysics, you'll need:

    * `Python <http://www.python.org/>`_ 2.5 or higher (2.6 highly recommended), although 3.x is not yet compatible with numpy or scipy.
    * `numpy <http://numpy.scipy.org>`_ 
    * `scipy <http://www.scipy.org/>`_
    
Follow the instructions on those sites, or, far easier, install them as packages from your operating system (e.g. apt-get or the synaptic GUI on Ubuntu, `Macports <http://www.macports.org/>`_ on OS X, etc.).  


Install
-------

Once you have the requirements satisfied, you have a few options for installing astropysics.

If you have `pip <http://pypi.python.org/pypi/pip>`_ (the new, better easy installer) or `easy_install/setuptools <http://pypi.python.org/pypi/setuptools>`_ (you should probably install pip...), just run either::

    pip install astropysics

or::

    easy_install astropysics

If you are on Ubuntu or a similar linux distribution, you may need to prefix those commands with ``sudo`` e.g. ``sudo pip install astropysics``.

If you are installing from source code, instead, just do::

    python setup.py install
    
where again, you may have to prefix the command with ``sudo`` if it doesn't work the first time.



Recommended packages
--------------------

A number of other packages are necessary for added functionality in astropysics (and all are recommended):
    * `Matplotlib <http://matplotlib.sourceforge.net/index.html>`_
        *highly recommended*, as it is necessary for all plotting (aside from the GUI applications). Install with::
        
            pip install matplotlib
            
    * `Ipython <http://ipython.scipy.org/>`_
        *highly recommended*, as it is a far better interactive shell than the python default and has many wonderful features. Install with::
        
            pip install ipython
            
        and if you have matplotlib, run as::
            
            ipython -pylab
            
    * `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
        Necessary for reading FITS files::
        
            pip install pyfits
        
            
    * `Traits <http://code.enthought.com/projects/traits/>`_, `TraitsGUI <http://code.enthought.com/projects/traits_gui/>`_, `Chaco <http://code.enthought.com/projects/chaco/>`_, and `Mayavi <http://code.enthought.com/projects/mayavi/>`_.  Alternatively, `ETS <http://code.enthought.com/projects/index.php>`_ is all bundled in one.
        Necessary for the interfaces in the :mod:`gui` modules::
        
            pip install ETS
            
        or::
        
            pip install traits
            pip install traitsGUI
            pip install chaco
            pip install mayavi
    
    * `vo.table <http://stsdas.stsci.edu/astrolib/vo/html/>`_
        Necessary to write VOTable files, and makes reading them much better, as well. Download from `<http://www.stsci.edu/trac/ssb/astrolib>`_ at the bottom of the page.
        
        
    Note that you can install all of these at once if you install astropysics using the following command::
    
        easy_install "astropysics[all]"
        
    or if there are problems installing the ETS tools, use::
    
        easy_install "astropysics[allnogui]"
        
    and as with all of these, you may need to prefix with ``sudo``.

Development/Advanced Install
----------------------------

If you want to be sure you have the most recent version (or want to help improve the code), you can pull the current development version from launchpad.  You must have `bazaar <http://bazaar.canonical.com/>`_ installed, and execute::

    bzr branch lp:astropysics astropysics-dev
    
and this will create a branch with the name 'astropysics-dev' containing the latest and greatest version of astropysics.  If at any time you want to update this version, go into the directory and do::

    bzr update
    
and re-install following the directions above.  If you plan on editing the astropysics source code (please do so, and submit patches/new features!), a useful way to immediately see changes without having to re-install every time is to use the command::

    python setup.py develop

possibly prefixed with ``sudo`` depending on your OS.  This will install links to the source code directory instead of copying over the source code, so any changes you make to a module can be seen just be doing ``reload(module)``.

   

Bug Reports
===========

The best place to report bugs is via the `launchpad bug tracker <https://bugs.launchpad.net/astropysics>`_.  That way they won't be forgotten unless an asteroid impact destroys all of humanity (or all the launchpad servers...).

Logo Image Credit
=================

The multiwavelength image of M81 was put together by the folks at the Chandra X-Ray Observatory (http://chandra.harvard.edu/photo/2008/m81/), and they credit: "X-ray: NASA/CXC/Wisconsin/D.Pooley & CfA/A.Zezas; Optical: NASA/ESA/CfA/A.Zezas; UV: NASA/JPL-Caltech/CfA/J.Huchra et al.; IR: NASA/JPL-Caltech/CfA".  The Python logo can be found at `<http://www.python.org/community/logos/>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

