#Copyright 2010 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
This package contains small external packages that are included with astropysics 
to simplify installation.  Currently, this includes:

    * ConfigObj v4.7.2  (http://www.voidspace.org.uk/python/configobj.html)

"""

def setup_module(module):
    """
    Used by doctests
    """
    from nose.plugins.skip import SkipTest
    raise SkipTest('Skipping tests for external modules/packages')