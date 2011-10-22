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
The :mod:`alg` module contains basic algorithms, numerical tricks, and generic
data processing tasks used in more than one place in `astropysics`.

.. seealso: :mod:`astropysics.stats`

.. todo:: examples?


.. add this back in if classes are implemented::

    Classes and Inheritance Structure
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. inheritance-diagram:: astropysics.utils.alg
       :parts: 1

Module API
^^^^^^^^^^

"""

from __future__ import division,with_statement
import numpy as np

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
    
def nearestsorted(a,val):
    """
    Search a sorted sequence for the nearest value.
    
    :param a: A sorted sequence to be searched.
    :type a: array-like
    :param val: The value(s) at which to find the nearest match in `a`.
    
    :returns: nearest,inds
        Where `nearest` are the elements of `a` that are nearest to the `val`
        values, and `inds` are the indecies into `a` that give those values.
    
    .. seealso::
        :func:`numpy.sort`, :func:`numpy.searchsorted`
    
    """
    a = np.asarray(a)
    val = np.asarray(val)
    scalar = val.shape == ()
    val = np.atleast_1d(val)
    
    i = np.searchsorted(a,val)
    
    iabove = i>=len(a)
    if np.any(iabove):
        i[iabove] -= 1
    ai = a[i]
    am = a[i-1]
    mbetter = np.abs(am-val)<np.abs(ai-val)
    i[mbetter] = i[mbetter]-1
    
    if scalar:
        return a[i][0],i[0]
    else:
        return a[i],i

def lin_to_log_rescale(val,lower=1,upper=3,base=10):
    """
    Linearly rescales input values onto the range [base^lower,base^upper] and
    then applies the requested base of logarithm to the rescaled values.
    
    :param val: input linear values
    :type val: array-like
    :param lower: lower value for output values
    :type lower: scalar
    :param upper: upper value for output values
    :type upper: scalar
    :param base: base of logarithm
    :type base: base of logarithm
    
    :returns: logarithm of rescaled input
    
    """
    if lower > upper:
        raise ValueError('lower must be less than upper')
    
    lower = base**lower
    upper = base**upper
    
    val = np.array(val,copy=False)
    #offset to [0,something]
    val = val - val.min()
    #reacale to [0,range]
    val *= ((upper-lower)/val.max())
    val += lower

    if base is None:
        return np.log(val)
    elif base==10:
        return np.log10(val)
    else:
        return np.log(val)/np.log(base)

def crossmask(x,threshold=0,belowtoabove=True):
    """
    Generates a boolean mask for the point where the input crosses or has
    complted crossing (if between) from below (or above) a threshold.
    
    If `belowtoabove` is True, the returned masks is for where the input
    transitions from below to above.  Otherwise, from above to below.
    
    :param x: input (will be flattened to 1D)
    :type x: array-like
    :param threshold: the transition value for where the crossings occur
    :type threshold: scalar
    :param belowtoabove: 
        If True, the returned mask is for where the input transitions from below
        to above. Otherwise, from above to below.
    :type belowtoabove: bool
    
    :returns: A mask that is True where crossing occurs, False everywhere else.
    :rtype: bool :class:`~numpy.ndarray`
    
    **Examples**
    
    >>> from numpy import array,where
    
    >>> xup = [-2,-1,0,1,2]
    >>> xdn = [2,1,0,-1,-2]
    
    >>> print crossmask(xup,0,True)
    [False False  True False False]
    >>> print crossmask(xup,-0.5,True)
    [False False  True False False]
    >>> print crossmask(xup,0.5,True)
    [False False False  True False]
    
    >>> print crossmask(xdn,0,True)
    [False False False False False]
    
    >>> print crossmask(xdn,0,False)
    [False False  True False False]
    >>> print crossmask(xdn,0.5,False)
    [False False  True False False]
    >>> print crossmask(xdn,-0.5,False)
    [False False False  True False]
    
    >>> xupdnup = [-2,-1,0,1,2,1,0,-1,-2,-1,0,1,2]
    >>> where(crossmask(xupdnup,0.5,True))
    (array([ 3, 11]),)
    >>> print array(xupdnup)[crossmask(xupdnup,0,True)]
    [0 0]


        
    """
    x = np.array(x,copy=False).ravel()
    if belowtoabove:
        a = x <  threshold
        b = x >= threshold
    else:
        a = x >  threshold
        b = x <= threshold
        
    mask = np.roll(a,1)&b
    mask[0] = False
    return mask
    
def nd_grid(*vecs):
    """
    Generates a grid of values given a sequence of 1D arrays. The inputs 
    will be converted to 1-d vectors and the output is an array with
    dimensions (nvecs,n1,n2,n3,...) varying only on the dimension 
    corresponding to its input order
    
    
    **Examples**::
    
        x = linspace(-1,1,10)
        y = x**2+3
        z = randn(13)
        result = nd_grid(x,y,z)
        xg,yg,zg = result
    
    result will be a (3,10,10,13) array, and each of xg,yg, and zg are 
    (10,10,13)
    """
    vecs = [np.array(v,copy=False).ravel() for v in vecs]
    shape = tuple([v.size for v in vecs])
    sz = np.prod(shape)
    
    arrs = []
    for i,(v,n) in enumerate(zip(vecs,shape)):
        newshape = list(shape)
        newshape.insert(0,newshape.pop(i))
        
        arr = np.repeat(v,sz/n).reshape(newshape)
        
        transorder = range(len(arr.shape))[1:]
        transorder.insert(i,0)
        
        arrs.append(arr.transpose(transorder))
    return np.array(arrs)
    
    

def centroid(val,axes=None,offset=None):
    """
    Convinience function calling :func:`moments` to get the first normalized
    moment (i.e. the centroid).
    
    :param val: n-d array for which to compute the centroid.
    :type val: array-like
    :param axes: 
        None to use default 0-based location scheme (see
        :func:`astropysics.utils.stats.moments`) or an array of poisition values
        for the centriod (must have as many elements as the dimension of `val`)
    :type axes: array-like or None
    :param offset: 
        A fixed offset to subtract from the data before centroiding, or a
        `bgmethod` like those accepted by :func:`estimate_background`.
        
    .. seealso:: :func:`astropysics.utils.stats.moments`
        
    """
    from .stats import moments
    
    if axes is not None and np.isscalar(axes[0]):
        axes = (axes,)
    return moments(val,1,axes,offset,True,False)



def sigma_clip(data,sig=3,iters=1,cenfunc='median',varfunc=np.var,maout=False):
    """
    This performs the sigma clipping algorithm - i.e. the data will be iterated
    over, each time rejecting points that are more than a specified number of
    standard deviations discrepant.
    
    :param data: input data (will be flattened to 1D)
    :type data: array-like
    :param sig: 
        The number of standard deviations to use as the clipping limit, or 
        the square root of the variance limit.
    :type sig: scalar
    :param iters: 
        The number of iterations to perform clipping for, or None to clip until
        convergence is achieved
    :param cenfunc: 
        The technique to compute the center for the clipping - may be any valid
        input to :func:`estimate_background`
    :param varfunc: 
        The method to compute the variance about the. Should take a 1D array as
        input and output a scalar. This will be compared to the square of the
        data as if it is the variance.
    :type varfunc: a callable
    :param maout: If True, return a masked array (see return value for details).
    :type maout: bool
    
    :returns: 
        A :class:`numpy.ma.Maskedarray` with the rejected points masked, if
        `maout` is True. If maout is False, a tuple (filtereddata,mask) is
        returned where the mask is False for rejected points (and matches the
        shape of the input).
    
    """
    data = np.array(data,copy=False)
    oldshape = data.shape
    data = data.ravel()
    
    mask = np.ones(data.size,bool)
    if iters is None:
        lastrej = sum(mask)+1
        while(sum(mask)!=lastrej):
            lastrej = sum(mask)
            do = data-estimate_background(data[mask],cenfunc)
            mask = do*do <= varfunc(data[mask])*sig**2
    else:
        for i in range(iters):
            do = data-estimate_background(data[mask],cenfunc)
            mask = do*do <= varfunc(data[mask])*sig**2
        
    if maout:
        return np.ma.MaskedArray(data,~mask,copy='maout'=='copy')
    else:
        return data[mask],mask.reshape(oldshape)
    
    
def estimate_background(arr,method='median'):
    """
    Estimates the background of the provided array following a technique
    specified by the `method` keyword:
    
    * 'median' : median of all the data
    * 'mean' : mean of all the data
    * '32' : 3*median - 2*mean
    * '2515' : 2.5*median - 1.5*mean
    * '21' : 2*median - 1*mean
    * a callable that takes a 1D array and outputs a scalar
    * a scalar : always returns that value
    * None : always returns a 0 background
    
    outputs a scalar background estimate
    """
    arr = np.array(arr,copy=False).ravel()
    
    if method is None:
        res = 0
    elif isinstance(method,basestring):
        if method == 'median':
            res = np.median(arr)
        elif method == 'mean':
            res = np.mean(arr)
        elif method == '32':
            res = (3*np.median(arr) - 2*np.mean(arr))
        elif method == '2515':
            res = (2.5*np.median(arr) - 1.5*np.mean(arr))
        elif method == '21':
            res = (2*np.median(arr) - np.mean(arr))
        else:
            raise ValueError('unrecognized estimate_background method string')
    elif np.isscalar(method) or (hasattr(method,'shape') and method.shape is tuple()):
        res = method
    elif callable(method):
        res = method(arr)
    else:
        raise ValueError('unrecognized estimate_background method type')
    
    return res

    
#<-----------------------Rotations--------------------------------------------->
    
def rotation_matrix(angle,axis='z',degrees=True):
    """
    Generate a 3x3 rotation matrix in cartesian coordinates for rotation about
    the requested axis.
    
    :param axis:
        Either 'x','y', 'z', or a (x,y,z) specifying an axis to rotate about. If
        'x','y', or 'z', the rotation sense is counterclockwise looking down the
        + axis (e.g. positive rotations obey left-hand-rule).
    :type axis: string or 3-sequence
    :param degrees: If True the input angle is degrees, otherwise radians.
    :type degrees: boolean
    
    :returns: A :class:`numpy.matrix` unitary rotation matrix.
    """
    from math import sin,cos,radians,sqrt
    if degrees:
        angle = radians(angle)
        
    
    
    if axis == 'z':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, s, 0),
                          (-s, c, 0),
                          ( 0, 0, 1)))
    elif axis == 'y':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( c, 0,-s),
                          ( 0, 1, 0),
                          ( s, 0, c)))
    elif axis == 'x':
        s = sin(angle)
        c = cos(angle)
        return np.matrix((( 1, 0, 0),
                          ( 0, c, s),
                          ( 0,-s, c)))
    else:
        x,y,z = axis
        w = cos(angle/2)
        
        #normalize
        if w == 1:
            x=y=z=0
        else:
            l = sqrt((x*x + y*y + z*z)/(1-w*w))
            x /= l
            y /= l
            z /= l
        
        wsq = w*w
        xsq = x*x
        ysq = y*y
        zsq = z*z
        return np.matrix((( wsq+xsq-ysq-zsq, 2*x*y-2*w*z, 2*x*z+2*w*y),
                          ( 2*x*y+2*w*z, wsq-xsq+ysq-zsq,2*y*z-2*w*x),
                          ( 2*x*z-2*w*y, 2*y*z+2*w*x, wsq-xsq-ysq+zsq)))
                          
                          
def angle_axis(matrix,degrees=True):
    """
    Computes the angle of rotation and the rotation axis for a given rotation
    matrix.
    
    :param matrix: the rotation matrix
    :type matrix: a 3x3 :class:`numpy.ndarray`
    :param degrees: if True, output is in degrees
    :type degrees: boolean
    
    :returns:
        an (angle,axis) tuple where the angle is in degrees if `degrees` is
        True, otherwise in radians
        
    """
    from math import sin,cos,acos,degrees,sqrt
    
    m = np.asmatrix(matrix)
    if m.shape != (3,3):
        raise ValueError('matrix is not 3x3')
    
    
    
    angle = acos((m[0,0] + m[1,1] + m[2,2] - 1)/2)
    denom = sqrt(2*((m[2,1]-m[1,2])+(m[0,2]-m[2,0])+(m[1,0]-m[0,1])))
    axis = np.array((m[2,1]-m[1,2],m[0,2]-m[2,0],m[1,0]-m[0,1]))/denom
    axis /= sqrt(np.sum(axis**2)) 
    
    if degrees:
        return degrees(angle),axis
    else:
        return angle,axis
    
def rotation_matrix_xy(x,y):
    """
    Computes the rotation matrix that moves the +z-axis (pole) to a new vector
    (x,y,z') where x and y are specified and z' is constrained by requiring
    the vector to be length-1.
    
    :param x: The x-value of the new pole in the starting coordinate system
    :type x: float
    :param y: The y-value of the new pole in the starting coordinate system
    :type y: float
    
    :returns: A :class:`numpy.matrix` unitary rotation matrix.
    
    """
    from math import sqrt
    
    xsq = x*x
    ysq = y*y
    xy = x*y
    
    z = sqrt(1-xsq-ysq)
    b = 1/(1+z)
    
    return np.matrix([[1-b*xsq,-b*xy,-x],
                      [-b*xy,1-n*ysq,-y],
                      [x,y,1-b*(xsq+ysq)]])
                      

def intrinsic_to_observed_ellipticity(ei,i,degrees=True):
    """
    Converts intrinsic ellipticity to observed where :math:`e^2 = 1-(b/a)^2`
    
    :param ei: intrinsic ellipticity
    :type ei: float or array-like
    :param i: inclination angle
    :type i: float or array-like
    :param degrees: if True, the inclination is assumed to be in degrees
    
    :returns: observed ellipticity
    """
    if degrees:
        i = np.radians(i)
        
    return 1 - ((1-ei)**2*np.sin(i)**2+np.cos(i)**2)**0.5
    
def observed_to_intrinsic_ellipticity(eo,i,degrees=True):
    """
    Converts observed ellipticity to intrinsic where :math:`e^2 = 1-(b/a)^2`
    
    :param eo: observed ellipticity
    :type eo: float or array-like
    :param i: inclination angle
    :type i: float or array-like
    :param degrees: if True, the inclination is assumed to be in degrees
    
    :returns: intrinsic ellipticity
    """
    if degrees:
        i = np.radians(i)
        
    return 1 - ((1-eo)-np.cos(i)**2)**0.5/np.sin(i)

del ABCMeta,abstractmethod,abstractproperty #clean up namespace
