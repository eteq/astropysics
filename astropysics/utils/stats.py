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
    

def binned_weights(values,n,log=False):
    """
    Produces an array of weights that are generated by subdividing the values
    into n bins such that each bin has an equal share of the total number of
    values.
    
    :param values: the input values
    :type values: array-like
    :param n: number of bins
    :type n: int
    :param log: 
        If True, the values are evenly-spaced on logarithmic intervals,
        otherwise, linear.
    :type log: bool
    
    :returns: An array of weights on [0,1] with shape matching `values`
    
    """
    
    if log:
        values = np.log(values).ravel()
    else:
        values = np.array(values,copy=False).ravel()
        
    mx,mi = np.max(values),np.min(values)
    
    n,edges = np.histogram(values)
    ws = np.zeros_like(values)
    
    wsr = ws.ravel()
    for i,w in enumerate(1.0/n):
        m = (edges[i]<=values) & (values<edges[i+1])
        wsr[m] = w 
    wsr[edges[-1]==values] = w
    
    return ws


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

class Pca(object):
    """
    A class for Principal Component Analysis (PCA).
    
    When mentioned in docstrings, `p` is the number of dimensions, and `N` is
    the number of data points.
    """
    _colors=('r','g','b','c','y','m','k') #defaults
    
    def __calc(self):
        A = self.A
        M = A-np.mean(A,axis=0)
        N = M/np.std(M,axis=0)
        
        self.M = M
        self.N = N
        self._eig = None
    
    def __init__(self,data,names=None):
        """
        
        :param data: 
            A `p` X `N` array-like object that is the data upon which PCA is to
            be performed. 
        :param names: 
            A sequence of names for each of the `p` axes (used in plots).
            
        """
        from warnings import warn
        self.A = A = np.array(data).T
        n,p = A.shape
        self.n,self.p = n,p
        if p > n:
            warn('p > n - intentional?')
        self._origA = A.copy()
        self.__calc()
        
        self._colors= np.tile(self._colors,int((p-1)/len(self._colors))+1)[:p]
        if names is not None and len(names) != p:
            raise ValueError('names must match data dimension')
        self.names = None if names is None else tuple([str(n) for n in names])
        
        
    def getCovarianceMatrix(self):
        """
        Computes  the covariance matrix for the dataset.
        
        :returns: A `p` x `p` covariance matrix.
        
        """
        return np.cov(self.N.T)
        
    def getEigensystem(self):
        """
        Computes the eigenvalues and eigenvectors of the data set.
        
        :returns: A 2-tuple (eigenvalues,eigenvectors) 
        """
        if self._eig is None:
            res = np.linalg.eig(self.getCovarianceMatrix())
            sorti=np.argsort(res[0])[::-1]
            res=(res[0][sorti],res[1][:,sorti])
            self._eig=res
        return self._eig
        
    def getEigenvalues(self):
        """
        Computes and returns a length `p` array with the eigenvalues.
        """
        return self.getEigensystem()[0]
        
    def getEigenvectors(self):
        """
        Computes and returns a `p` x `p` array with the eigenvalues.
        """
        return self.getEigensystem()[1]
    
    def getEnergies(self):
        """
        Computes and returns a length `p` array with the eigenvalues normalized
        so that they sum to 1.
        """
        v=self.getEigenvalues()
        return v/np.sum(v)
        
    def plot2d(self,ix=0,iy=1,clf=True):
        """
        Generates a 2-dimensional plot of the data set and principle components 
        using matplotlib.
        
        ix specifies which p-dimension to put on the x-axis of the plot
        and iy specifies which to put on the y-axis (0-indexed)
        """
        import matplotlib.pyplot as plt
        x,y=self.N[:,ix],self.N[:,iy]
        if clf:
            plt.clf()
        plt.scatter(x,y)
        vals,evs=self.getEigensystem()
        #evx,evy=evs[:,ix],evs[:,iy]
        xl,xu=plt.xlim()
        yl,yu=plt.ylim()
        dx,dy=(xu-xl),(yu-yl)
        for val,vec,c in zip(vals,evs.T,self._colors):
            plt.arrow(0,0,val*vec[ix],val*vec[iy],head_width=0.05*(dx*dy/4)**0.5,fc=c,ec=c)
        #plt.arrow(0,0,vals[ix]*evs[ix,ix],vals[ix]*evs[iy,ix],head_width=0.05*(dx*dy/4)**0.5,fc='g',ec='g')
        #plt.arrow(0,0,vals[iy]*evs[ix,iy],vals[iy]*evs[iy,iy],head_width=0.05*(dx*dy/4)**0.5,fc='r',ec='r')
        if self.names is not None:
            plt.xlabel('$'+self.names[ix]+'/\\sigma$')
            plt.ylabel('$'+self.names[iy]+'/\\sigma$')
        
    def plot3d(self,ix=0,iy=1,iz=2,clf=True):
        """
        Generates a 3-dimensional plot of the data set and principle components
        using mayavi.  
        
        ix, iy, and iz specify which of the input p-dimensions to place on each of
        the x,y,z axes, respectively (0-indexed).
        """
        import enthought.mayavi.mlab as M
        if clf:
            M.clf()
        z3=np.zeros(3)
        v=(self.getEigenvectors()*self.getEigenvalues())
        M.quiver3d(z3,z3,z3,v[ix],v[iy],v[iz],scale_factor=5)
        M.points3d(self.N[:,ix],self.N[:,iy],self.N[:,iz],scale_factor=0.3)
        if self.names:
            M.axes(xlabel=self.names[ix]+'/sigma',ylabel=self.names[iy]+'/sigma',zlabel=self.names[iz]+'/sigma')
        else:
            M.axes()
        
    def sigclip(self,sigs):
        """
        clips out all data points that are more than a certain number
        of standard deviations from the mean.  
        
        sigs can be either a single value or a length-p sequence that
        specifies the number of standard deviations along each of the 
        p dimensions.
        """
        if np.isscalar(sigs):
            sigs=sigs*np.ones(self.N.shape[1])
        sigs = sigs*np.std(self.N,axis=1)
        n = self.N.shape[0]
        m = np.all(np.abs(self.N) < sigs,axis=1)
        self.A=self.A[m]
        self.__calc()
        return n-sum(m)
        
    def reset(self):
        self.A = self._origA.copy()
        self.__calc()
        
        
    def project(self,vals=None,enthresh=None,nPCs=None,cumen=None):
        """
        projects the normalized values onto the components
        
        enthresh, nPCs, and cumen determine how many PCs to use
        
        if vals is None, the normalized data vectors are the values to project.
        Otherwise, it should be convertable to a p x N array
        
        returns n,p(>threshold) dimension array
        """
        nonnones = sum([e != None for e in (enthresh,nPCs,cumen)])
        if nonnones == 0:
            m = slice(None)
        elif nonnones > 1:
            raise ValueError("can't specify more than one threshold")
        else:
            if enthresh is not None:
                m = self.energies() > enthresh
            elif nPCs is not None:
                m = slice(None,nPCs)
            elif cumen is not None:
                m = np.cumsum(self.energies()) <  cumen
            else:
                raise RuntimeError('Should be unreachable')
        
        if vals is None:
            vals = self.N.T
        else:
            vals = np.array(vals,copy=False)
            if self.N.T.shape[0] != vals.shape[0]:
                raise ValueError("shape for vals doesn't match")
        proj = np.matrix(self.getEigenvectors()).T*vals
        return proj[m].T
            
    def deproject(self,A,normed=True):
        """
        input is an n X q array, where q <= p
        
        output is p X n
        """
        A=np.array(A,ndmin=2)
        n,q = A.shape
        p = self.A.shape[1]
        if q > p :
            raise ValueError("q > p")
        
        evinv=np.linalg.inv(np.matrix(self.getEigenvectors()).T)
        
        zs = np.zeros((n,p))
        zs[:,:q]=A
        
        proj = evinv*zs.T
        
        if normed:
            return np.array(proj.T).T
        else:
            mns=np.mean(self.A,axis=0)
            sds=np.std(self.M,axis=0)
            return (np.array(proj.T)*sds+mns).T
    
    def subtractPC(self,pc,vals=None):
        """
        pc can be a scalar or any sequence of pc indecies
        
        if vals is None, the source data is self.A, else whatever is in vals
        (which must be p x m)
        """
        if vals is None:
            vals = self.A
        else:
            vals = vals.T
            if vals.shape[1]!= self.A.shape[1]:
                raise ValueError("vals don't have the correct number of components")
        
        pcs=self.project()
        zpcs=np.zeros_like(pcs)
        zpcs[:,pc]=pcs[:,pc]
        upc=self.deproject(zpcs,False)
        
        A = vals.T-upc
        B = A.T*np.std(self.M,axis=0)
        return B+np.mean(self.A,axis=0)



del ABCMeta,abstractmethod,abstractproperty #clean up namespace
