Developer Guidelines for Astropysics
====================================

Astropysics welcomes contributions from anyone interested, from simple bugfixes to contributing new functionality.  
If you want to check out the most recent version to contribute (or just want the latest and greatest), you can pull the current development version from  the `google code page <http://code.google.com/p/astropysics/>`_.  


Getting the Code
----------------

You will need to have `mercurial <http://mercurial.selenic.com/>`_ installed.  Once you do, simply execute::

    hg clone https://astropysics.googlecode.com/hg/ astropysics-dev
    
This will create a directory with the name ``astropysics-dev`` containing the latest and greatest version of astropysics.  
If at any time you want to update this version, go into the directory and do::

    hg pull
    hg update
    
then re-install following the directions above.  If you plan on editing the astropysics source code (please do so, and submit patches/new features!), a useful way to immediately see changes without having to re-install every time is to use the command::

    python setup.py develop

possibly prefixed with ``sudo`` depending on your OS.  This will install links to the source code directory instead of copying over the source code, so any changes you make to a module can be seen just be doing ``reload(module)``.

Cloning the Repository to Submit Code
-------------------------------------

If you intend to regularly contribute changes or patches to astropysics, a more convinient way to submit changes is with a public clone of the main astropysics repository.
Go to the `source tab  <http://code.google.com/p/astropysics/source/checkout>`_ of the `google code project <http://code.google.com/p/astropysics>`_, and click on the ``create clone`` button.  
Fill in the necessary information, and clone *your* repository to your computer instead of the main astropysics repository.  
Make your changes, using ``hg commit -m "a message"`` to describe changes as you make them.
You can then use ``hg push`` to send changes back to your repository on google code, and those can easily be merged with the main repository with a pull request.
   
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
        """


Classes with public attributes can document using the sphinx construct for documenting class attributes::

    class MyClass(object):
        
        #: Documentation for :attr:`value` attribute.
        value = None
    
        def __init__(self,value):
            self.value = value
        
        
