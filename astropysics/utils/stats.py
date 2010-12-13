#Copyright 2009 Erik Tollerud
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

The :mod:`stats` module contains classes and functions for statistics and
statistical analysis. These tools in this module are mostly inspecific to
astrophysics - the applications are all in the relevant other modules.

.. seealso:: 
    :mod:`scipy.stats` - :mod:`astropysics,utils.stats` is intended only to
    provide utilites and interfaces that are not present in :mod:`scipy.stats` -
    when possible, :mod:`scipy.stats` should be used.

.. todo:: examples?

.. add this back in if classes are implemented::

    Classes and Inheritance Structure
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. inheritance-diagram:: astropysics.utils.stats
       :parts: 1

Module API
^^^^^^^^^^

"""



from __future__ import division,with_statement
import numpy as np
from scipy import stats as spystats
import re as _re

try:
    #requires Python 2.6
    from abc import ABCMeta
    from abc import abstractmethod
    from abc import abstractproperty
    from collections import MutableMapping
except ImportError: #support for earlier versions
    abstractmethod = lambda x:x
    abstractproperty = property
    ABCMeta = type
    MutableMapping = type
    


def moments(arr,ms,axes=None,bgmethod=None,norm=True,std=False):
    """
    Compute the moments of the provided n-d array.  That is 
    
    :param arr: The input array for which the moments are desired.
    :param ms: 
        The desired order of the moment to be returned.  Can be:
            
            * A scalar
                The same order will be used for each dimension, and the return 
                type will be an array of the moment along each dimension.
            * A sequence 
                The sequence must match the number of dimensions in `arr`, and
                specifies the order along each dimension.  e.g. for a 3D 
                cartesian array, ms=[1,3,2] means the moemnt comupted will be
                :math:`x^1 y^3 z^2`.
    :param axes: 
        If None, the spatial axes are taken to be 0,1,2...,nd-1 for each of the
        dimenstions in the input. Otherwise, `axes` must be a seqence of arrays
        specifying the spatial locations along each axis. This sequence must be
        of length matching the number of dimensions in the input, and Each
        element must be of the same length as the corresponding input array
        dimension.
    :param bgmethod: 
        A background-estimation method (see :func:`estimate_background` for
        options), a scalar to subtract from the array, or None to do no
        subtraction.
    :param bool norm: 
        If True, the moment will be normalized - i.e. ``sum(x^m*arr)/sum(arr)``
        instead of just ``sum(x^m*arr)``
    :param bool std: 
        If True, the output will be standardized (mth moment divided by standard
        deviation to the mth power)
    
    :returns: 
        Either the computed moment if `ms` is a sequence, or a 1D array of
        moments for each dimension if `ms` is a scalar.
    
    """
    from .alg import estimate_background
   
    arr = np.array(arr,copy=False)
    if bgmethod:
        arr = arr - estimate_background(arr,bgmethod)
    shp = arr.shape
    
    #setup/check axes
    if axes is None:
        axes = [np.arange(s) for s in shp]
    elif len(axes) != len(shp):
        raise ValueError('incorrect number of axes provided')
    else:
        axmatch = np.array([len(ax)==s for ax,s in zip(axes,shp)])
        if np.any(~axmatch):
            raise ValueError('axes %s do not match input array'%np.where(~axmatch))
        
    newax = []
    for i,ax in enumerate(axes):
        shi = np.ones(len(shp))
        shi[i] = len(ax)
        newax.append(ax.reshape(shi))
    axes = newax
    
    if np.isscalar(ms):
        res = np.array([np.sum(arr*ax**ms) for ax in axes])
    else:
        if len(ms) != len(shp):
            raise ValueError('moment sequence does not match data')
        #bcast = np.broadcast_arrays(*[ax**m for m,ax in zip(ms,axes)])
        #res = np.sum(arr*np.prod(bcast,axis=0))
        axprod = reduce(np.multiply,[ax**m for m,ax in zip(ms,axes)])
        res = np.sum(arr*axprod)
    
    if std:
        if np.isscalar(ms):
            res = res/np.array([np.std(ax)**ms for ax in axes])
        else:
            res = res/np.prod([np.std(ax)**m for m,ax in zip(ms,axes)])
    
    if norm:
        res/=np.sum(arr)
    return res    
    

#<--------------------------Robust statistics---------------------------------->    
def interquartile_range(values,scaletonormal=False):
    """
    Computes the interquartile range for the provided sequence of values, a
    more robust estimator than the variance.
    
    :param values: the values for which to compute the interquartile range
    :type values: array-like, will be treated as 1D
    :param scaletonormal: Rescale so that a normal distribution returns 1
    :type scaletonormal: bool
    
    :returns: the interquartile range as a float
    
    """
    from scipy.stats import scoreatpercentile
    from scipy.special import erfinv
    
    x = np.array(values,copy=False).ravel()
    res = scoreatpercentile(x,75) - scoreatpercentile(x,25)
    if scaletonormal:
        nrm = 8**0.5*erfinv(.5)
        return res/nrm
    else:
        return res
    
def median_absolute_deviation(values,scaletonormal=False,cenestimator=np.median):
    """
    Computes the median_absolute_deviation for the provided sequence of values, 
    a more robust estimator than the variance.
    
    :param values: the values for which to compute the MAD
    :type values: array-like, will be treated as 1D
    :param scaletonormal: Rescale the MAD so that a normal distribution is 1
    :type scaletonormal: bool
    :param cenestimator: 
        A function to estimate the center of the values from a 1D array of
        values. To actually be the "median" absolute deviation, this must be
        left as the default (median).    
    :type cenestimator: callable
    
    :returns: the MAD as a float
    
    """
    from scipy.special import erfinv
    
    x = np.array(values,copy=False).ravel()
    res = np.median(np.abs(x-np.median(x)))
    
    if scaletonormal:
        nrm = (2**0.5*erfinv(.5))
        return res/nrm
    else:
        return res
    
def biweight_midvariance(values,influencescale=9,cenestimator=np.median):
    """
    Computes the biweight midvariance of a sequence of data points, a robust 
    statistic of scale.  
    
    For normal and uniform distributions, it is typically close to, but a bit
    above the variance.
    
    :param values: the values for which to compute the biweight
    :type values: array-like, will be treated as 1D
    :param influencescale: 
        The number of MAD units away at which a data point has no weight
    :type influencescale: int
    :param cenestimator: 
        A function to estimate the center of the values from a 1D array of
        values. To be a true standard biweight midvariance, this must be left as
        the default (median).
    :type cenestimator: callable
    
    :returns: biweight,median tuple (both floats)
    
    """
       
    x = np.array(values,copy=False).ravel()
    
    Q = cenestimator(x)
    MAD = median_absolute_deviation(x,cenestimator=cenestimator)
       
    ui = (x-Q)/(influencescale*MAD)
    uisq=ui**2
    
    I = uisq <= 1
    
    numer = np.sum(I*(x-Q)**2*(1-uisq)**4)
    denom = np.sum(I*(1-uisq)*(1-5*uisq))**2
    
    
    return x.size*numer/denom,Q

del ABCMeta,abstractmethod,abstractproperty #clean up namespace
