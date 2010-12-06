"""
The astropysics interactive IPython configuration file for ipython versions >=0.11
"""

load_subconfig('ipython_config.py')
c = get_config()

lines = """
import numpy
import numpy as np
from numpy import *

import scipy
from scipy import stats,optimize,ndimage,integrate,interpolate,special

import pyfits
try:
    import astropysics
    from astropysics import phot,spec,coords,models,constants,objcat,obstools,io,plotting
except ImportError:
    print "Unable to start astropysics profile, is astropysics installed?"

"""

mpllines = """
import matplotlib
matplotlib.interactive(True)
matplotlib.use('{MPLBACK}')
wxapp = %gui {GUITK}
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
"""

c.Global.exec_lines.append(mpllines)
c.Global.exec_lines.append(lines)


