"""
Tests exploring how well centroiding performs in various conditions
"""

from __future__ import division
import numpy as np
from numpy.random import poisson,randn,rand

from astropysics import models,utils

def do_centroid(model,size,sampling):
    
    F = poisson(m.pixelize(-size/2,size/2,-size/2,size/2,
                           size,size,sampling=sampling))
    x = np.linspace(-size/2,size/2,size)
    y = np.linspace(-size/2,size/2,size)

    I = np.sum(F,axis=1)
    J = np.sum(F,axis=0)

    Isub = I-np.sum(I)/size
    Isub[Isub<0] = 0
    Jsub = J-np.sum(J)/size
    Jsub[Jsub<0] = 0
    
    truex,truey = model.mux,model.muy
    xc_howell = np.sum(Isub*x)/np.sum(Isub)
    yc_howell = np.sum(Jsub*y)/np.sum(Jsub)
    res_howell = (xc_howell-truex,yc_howell-truey)
    
    resx,resy = utils.centroid(F,axes=(x,y))
    cenres = (resx-truex,resy-truey)
    
    Fsub = F-np.mean(F)
    Fsub[Fsub<0] = 0
    meancenx,meanceny = utils.centroid(Fsub,axes=(x,y))
    meancenres = (meancenx-truex,meanceny-truey)

    return res_howell,meancenres,cenres

n = 100
size = 50
sampling = 10
xs = rand(n)*20-10 #[-10,10]
ys = rand(n)*20-10 #[-10,10]
sigs = rand(n)*8+2 #[2,10]
A = 1000

print 'Gaussian:'
res = []
for i,(x,y,s) in enumerate(zip(xs,ys,sigs)):
    print 'doing',i+1,'of',n
    m = models.Gaussian2DModel(A=A,mux=x,muy=y)
    m.sigx = m.sigy = s
    res.append(do_centroid(m,size,sampling))
res = np.array(res)
ressep = np.sum(res**2,axis=-1)**0.5
hm = np.mean(ressep[:,1]-ressep[:,0])
hc = np.mean(ressep[:,2]-ressep[:,0])
mc = np.mean(ressep[:,2]-ressep[:,1])
print 'howell','better than' if hm>0 else 'worse than','meancen',hm
print 'howell','better than' if hc>0 else 'worse than','bgcentroid',hc
print 'meancen','better than' if mc>0 else 'worse than','bgcentroid',mc

