Developer Guidelines for Astropysics
====================================

Astropysics welcomes contributions from anyone interested, from simple bugfixes to contributing new functionality.  
If you want to check out the most recent version to contribute (or just want the latest and greatest), you can pull the current development version from  the `github page <http://github.com/eteq/astropysics>`_.  


Getting the Code
----------------

You will need to have `git <http://git-scm.com/>`_ installed - it is available for all major OSes.  Once you have it working, simply execute::

    git clone git://github.com/eteq/astropysics.git astropysics-dev
    
This will create a directory with the name ``astropysics-dev`` containing the latest and greatest version of astropysics.  
If at any time you want to update this version, go into the directory and do::

    git pull
    
then re-install following the directions at :doc:`install`.  Usage of git to edit source code is well-documented elsewhere - github's `help page <http://help.github.com/>`_ covers the important basics. 

If you plan on editing the astropysics source code (please do so, and submit patches/new features!), a useful way to immediately see changes without having to re-install every time is to use the command::

    python setup.py develop

possibly prefixed with ``sudo`` depending on your OS.  This will install links to the source code directory instead of copying over the source code, so any changes you make to a module can be seen just be doing ``reload(module)``.

Cloning the Repository to Submit Code
-------------------------------------

If you intend to regularly contribute changes or patches to astropysics, a great way to submit changes is to use github's "fork" feature.  Just go to the  `astropysics page <http://github.com/eteq/astropysics>`_
and click the "fork" button in the upper-right.  Make a new branch in your fork with whatever change you want to make, and when you've fixed the bug of added whatever new fancy
thing you want to add, just submit a pull request (button in the upper-right) to the main astropysics project.
   
Coding Style Guidelines
-----------------------

Naming Conventions
^^^^^^^^^^^^^^^^^^

For general coding style, `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ provides the main coding style for astropysics, with an important exception: PEP8 calls for methods (i.e. functions that are in the scope of a class) to use the ``lower_case_with_underscores`` convention. 
Astropysics instead uses ``camelCase`` (e.g. first letter of each word upper case, first letter of the method lower case) for method names.  Functions that are not methods (i.e. directly in a module scope) remain ``lowercase_with_underscores``.
This allows methods and functions to be easily distinguished (but keeps method names distinct from class names, for which the first letter is always upper case).  This could change in the future to be fully PEP8 compliant, but for now, given the already existing codebase, ``camelCase`` it is.  

To summarize, the naming conventions are:

* Module names are always ``lowercase``.
* Class names are always ``CamelCaseFirstLetterUppercase``.
* Methods (including static methods and class methods) are always ``camelCaseFirstLetterLowercase``.
* Functions are always ``lowercase_with_underscore``.
* Variables and class attributes should be ``lowercase`` and kept to a single word if possible.
* Private/internal functions, methods, or variables should be ``_leadingUnderscore`` with the appropriate convention for the type.  Python and sphinx both know to hide these unless you specifically request them.  Python also supports ``__doubleLeadingUnderscore`` for private class methods (the double-underscore is `mangled <http://docs.python.org/tutorial/classes.html#private-variables>`_), but this generally just leads to confusion if you're not careful, so it should be avoided unless there's some very good reason.


Documentation Standards
^^^^^^^^^^^^^^^^^^^^^^^

Documentation should be provided with every object, using `sphinx <http://sphinx.pocoo.org/>`_ REStructuredText markup. 
Functions and methods should use `info field lists <http://sphinx.pocoo.org/domains.html#info-field-lists>`_ To specify input parameters, return values, and exceptions (where relevant). Below is an example of the standard format::

    def a_function(arg1,arg2=None,flag=False):
        """
        A short description of the function.
        
        Multiple paragraphs if necessary, i.e. background is needed.
        
        :param arg1: 
            This argument is important input data, although I'm not specifying
            exactly what it's type is (maybe it'd duck-typed?)  Also, the 
            description is more than one line, so it has to start underneath
            and indented.
        :param arg2: This argument is an optional input.
        :type arg2: You can specify a type here if you want.
        :param bool flag: You can also give the type in param if it fits.
        
        :except ValueError: 
            If you raise an exception, specify here what type it is and why.
            
        :returns: A description of the return value, if there is one.
    
        **Examples**
        
        If an examples are needed, they should go here, ideally in doctest 
        format so that they can be used as tests:
    
        >>> inpt = 'something for a_function'
        >>> a_function(inpt,flag=True)
        'whatever a_function should output'
        
        """


Classes with public attributes can document using the sphinx construct for documenting class attributes::

    class MyClass(object):
        
        #: Documentation for :attr:`value` attribute.
        value = None
    
        def __init__(self,value):
            self.value = value
        
Testing Astropysics
-------------------

There is a test suite that should be periodically run to ensure everything that has tests is still working correctly.  It requires `nose <http://pypi.python.org/pypi/nose>`_.
It can be run from the astropysics source directory (where setup.cfg lives) with the command::

    nosetests

Note that this is also set up to easily debug in the event that some of the tests fail.  Simply do::

    nosetest --failed

And nose will only run those tests that failed the last time around.  If you want to run a particular test, do::

    nostest --with-id 3

Where the '3' can be replaced by whatever number test you want.

When writing  functionality in astropysics, it's a good idea to add tests.  These should go in the 'tests' directory, and should have module names with the word 'test' in them, along with the function names themselves.
This naming is necessary to allow nose to find all the tests.  Alternatively, snippets of code as they would appear on the python interpreter (*with* output) can be placed directly in the docstrings, and they will be automatically included in the tests.
        
