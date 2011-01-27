#Copyright 2008 Erik Tollerud
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
This package contains a variety of utilities and algorithms for processing
and visualizing astronomical data.

**Modules**

* constants: physical constants and conversions, Cosmology objects for choosing 
  related sets of constants
* coords: Astronomical coordinate systems, distance measurements,
  and related objects
* obstools: Tools and corrections for observations (mostly optical)
* models: Fitting functions/models and related calculations
* objcat: Object Catalog objects and functions
* phot: Photometry objects and related functions
* spec: Spectrum objects and related functions
* ccd: Functions and tools for processing CCD images.

"""

#If True, the add_docs and similar functions will only replace with empty 
#strings - this is used by sphinx to 
_ignore_add_docs = False

#Large-scale TODO list:

#Add Obsplots 
#ccd tweaks/todos 
#instrument registry along with site
#objcat testing w/ M31 catalog
#ZODB/Web objcat integration
#pipeline gui and mid-saving
#Phot reworking/clean up docs