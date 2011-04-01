""" 
The astropysics interactive IPython configuration file for ipython versions < .11
"""

import IPython.ipapi
#this is only necessary if ipythonrc isn't in use, so if it is, this import can be commented out
import ipy_defaults

ip = IPython.ipapi.get()

ip.ex("from __future__ import division")
ip.ex("import numpy")
ip.ex("import numpy as np")
ip.ex("from numpy import *")
ip.ex("from numpy.random import rand,randn,randint")

ip.ex("import scipy")
ip.ex("from scipy import stats,optimize,ndimage,integrate,interpolate,special")

#import pyfits and asciitable if they are present
try:
    ip.ex("import pyfits")
except ImportError:
	pass
	try:
ip.ex("import asciitable")
	except ImportError:
		pass
		
try:
    ip.ex("import astropysics")
    #import typical modules
    ip.ex("from astropysics import phot,spec,coords,models,constants,objcat,obstools,plotting,utils")
    try:
        ip.ex("from astropysics import gui")
    except ImportError:
        pass #this just means traits isn't installed
except ImportError:
    print "Unable to start astropysics profile, try re-running astpys-setup (or re-installing astropysics)"
    



