from __future__ import division
from math import pi
import numpy as np

try:
    import spylot as sp
except ImportError:
    print 'Spylot not found - many spec package functions will not work correctly'