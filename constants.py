from __future__ import division
from math import pi
import numpy as np

#all in cgs including gaussian esu for charge
G=6.673e-8
mp=1.67262171e-24
me=9.1093897e-28
e=4.8032068e-10
Ms=1.9891e33
Mj=1.8986e30
Me=5.9742e27
Rs=6.96e10
Rj=7.1492e9
Re=6.371e8
De=1.49597887e13
kb=1.3807e-16 
#sigma=5.6704e-5
c=2.99792458e10

#unit conversion factors
evtoerg=1.60217646e-12
secperyr=3.15576926e7
secpergyr=secperyr*1e9
pcpercm=1/3.08568025e18
pcperly=1/3.26

def flambda_to_fnu_l(flambda,lamb):
    return flambda*lamb*lamb/c

def fnu_to_flambda_l(fnu,lamb):
    return fnu*c/lamb/lamb

def flambda_to_fnu_n(flambda,nu):
    return flambda*c/nu/nu

def fnu_to_flambda_n(fnu,nu):
    return fnu*nu*nu/c

