Astropysics
===========

A python library for astronomy/astrophysics calculations and data analysis.

Web Page: http://packages.python.org/Astropysics/

License: Apache License v2.0

Note that Astropysics is now in "maintainence mode" only.  That is,
further development is not expected, as my efforts have shifted to the
`Astropy project <http://www.astropy.org>`_.  It is a much larger effort
that contains most of the functionality in Astropysics and much more.
Astropysics will remain so that code written against it will continue to
function, though, and I will continue to accept bug fixes as necessary.

Installation
------------

If you have numpy (http://numpy.scipy.org/) and scipy (http://www.scipy.org/) installed, just do::

  python setup.py install

On some linux distributions, this may need to be::

  sudo python setup.py install
  
Note also, if you are using `pip` to install, the latest version might not work do to some issues with SSL. In that case you may have luck using the latest dev version::

  pip install git+https://github.com/eteq/astropysics.git
  

Documentation
-------------

Documentation (in the "docs" directory) is meant to be used with sphinx (http://sphinx.pocoo.org).  Typical usage is::

  python setup.py build_sphinx

to build html documentation and place it in docs/_build/html.  Note that this requires Sphinx >=1.0

Source Distribution Directory Structure
---------------------------------------

* astropysics/ - source code for all astropysics modules
* docs/ - Sphinx documentation and example code
* logo/ - astropysics logo in various forms
* licenses/ - copyright/license notices for all open source code
* scripts/ - command-line scripts installed with astropysics
* tests/ - unit tests (for use with Nose) and algorithmic experiments (Note unit tests require `networkx <http://pypi.python.org/pypi/networkx>`_ and `pymodelfit <http://pypi.python.org/pypi/PyModelFit>`_ to run).

