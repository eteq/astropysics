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

import pyfits
try:
    import astropysics
    from astropysics import phot,spec,coords,models,constants,objcat,obstools,plotting,utils
except ImportError:
    print "Unable to start astropysics profile, try re-running astpys-setup (or re-installing astropysics)"

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


