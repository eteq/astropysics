"""
The astropysics interactive IPython configuration file for ipython versions >=0.11
"""

load_subconfig('ipython_config.py')
c = get_config()

lines = """
import numpy
import numpy as np
from numpy import *
from numpy.random import rand,randn,randint

import scipy
from scipy import stats,optimize,ndimage,integrate,interpolate,special

try:
    import astropysics
    from astropysics import phot,spec,coords,models,constants,objcat,obstools,plotting,utils
except ImportError:
    print "Unable to start astropysics profile, try re-running astpys-setup (or re-installing astropysics)"

#silently ignore pyfits and asciitable if they are not present,as they are optional
try:
	import pyfits
except ImportError:
	pass
try:
	import asciitable
except ImportError:
	pass

"""

mpllines = """
import matplotlib
matplotlib.interactive(True)
matplotlib.use('{MPLBACK}')
guiapp = %gui {GUITK}
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
"""

c.Global.exec_lines.append(mpllines)
c.Global.exec_lines.append(lines)


