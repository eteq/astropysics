.. Astropysics documentation master file, created by
   sphinx-quickstart on Wed Feb 24 00:07:20 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Astropysics: astrophysics utilities for python
==============================================
:Author: `Erik Tollerud <http://www.physics.uci.edu/~etolleru/>`_

Astropysics is a library containing a variety of utilities and algorithms for
reducing, analyzing, and visualizing astronomical data. Best of all, it
encourages the user to leverage the existing capabilities of Python to make this
quick, easy, and as painless as cutting-edge science can even actually be. There
do exist other Python packages with some of the capabilities of this project,
but the goal of this project is to integrate all these tools together and make
them interact in the most straightforward ways possible.

(And to that end, if you are running one of those other projects, I'd love to 
help integrate our projects into a common framework!)

.. _contents:

Contents
========

Astropysics is divided into two major subcomponents - the core modules that
contain functions and classes to the calculations and organize data, and the gui
module that contains a number of useful small-scale astronomy applications.

.. toctree::
   :maxdepth: 2
   
   install
   getstarted
   coremods/intro
   gui/intro
   develop


Quick Install
=============

See :doc:`install` for full install instructions, including prerequisite
packages.

To install a current release of astropysics, the simplest approach is::

    pip install astropysics

(on unix-like systems or OS X, add "sudo " before this command)

If you want the must up-to-date (possible unstable) version, do::

    hg clone https://astropysics.googlecode.com/hg/ astropysics-dev
    cd astropysics-dev
    python setup.py develop

(note that `mercurial <http://mercurial.selenic.com/>`_ must be installed, and
on some systems the last command may need to have "sudo " at the beginning) 

You can also alter the source code if you use this approach (see :doc:`develop`
for guidelines of working contributing source code).

In either case, afterwords you can run::

    astpys-setup
    
to install optional packages and setup the environment.

       
Bug Reports
===========

The best place to report bugs is via the `google code bug tracker <http://code.google.com/p/astropysics/issues>`_.  That way they won't be forgotten unless an asteroid impact destroys all of google's servers.

People
============

Lead Developer
--------------

* Erik Tollerud

Contributors
------------

* Sebastian Gurevich
* Alex Hagen
* Frederic Grollier


Logo Image Credit
=================

The multiwavelength image of M81 was put together by the folks at the Chandra X-Ray Observatory (http://chandra.harvard.edu/photo/2008/m81/), and they credit: "X-ray: NASA/CXC/Wisconsin/D.Pooley & CfA/A.Zezas; Optical: NASA/ESA/CfA/A.Zezas; UV: NASA/JPL-Caltech/CfA/J.Huchra et al.; IR: NASA/JPL-Caltech/CfA".  The Python logo can be found at `<http://www.python.org/community/logos/>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

