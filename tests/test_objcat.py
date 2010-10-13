#!/usr/bin/env python
from __future__ import division,with_statement
from astropysics.constants import pi
import numpy as np
from astropysics.objcat import StructuredFieldNode,Catalog,Field

class Test1(StructuredFieldNode):
    num = Field('num',(float,int),(4.2,1,2))
    num2 = Field('num2',(float,int),(5.6,1.5,0.3))
    
    @StructuredFieldNode.derivedFieldFunc
    def f(num='num',num2='num2'):
        return 2*num+1+num2
    
    @StructuredFieldNode.derivedFieldFunc(ferr=True)
    def f2(num='num',num2='num2'):
        num,unum,lnum = num
        val = 2/num
        uerr = unum*2*num**-2
        lerr = lnum*2*num**-2
        return val,uerr,lerr
    
    @StructuredFieldNode.derivedFieldFunc
    def f3(num='num',num2='num2'):
        return num2*num
    
class Test2(StructuredFieldNode):
    val = Field('val',float,4.2)
    
    @StructuredFieldNode.derivedFieldFunc(num='num')
    def d1(val='val',d2='d2'):
        return val+np.exp(d2)
    
    @StructuredFieldNode.derivedFieldFunc(num='num')
    def d2(d1='d1'):
        if d1 is not None:
            return np.log(d1)
    

def test_cat():
    c = Catalog()
    t1 = Test1(c)
    t2 = Test1(c)
    t2['num'] = ('testsrc1',7)
    
    subc = Catalog('sub-cat',c)
    ts1 = Test1(subc)
    ts1['num'] = ('testsrc2',8.5)
    ts2 = Test1(subc,num=('testsrc2',12.7),f=('testows',123.45))
    ts3 = Test1(subc)
    ts3['num'] = ('testsrc2',10.3)
    
    return c

def test_sed():
    from numpy.random import randn,rand
    from numpy import linspace
    from astropysics.spec import Spectrum
    from astropysics.phot import PhotObservation
    from astropysics.models import BlackbodyModel
    from astropysics.objcat import SEDField
    
    f = SEDField()
    scale = 1e-9
    
    f['s1'] = Spectrum(linspace(3000,8000,1024),scale*(randn(1024)/4+2.2),scale*rand(1024)/12)
    m = BlackbodyModel(T=3300)
    m.peak = 2.2
    x2 = linspace(7500,10000,512)
    err = randn(512)/12
    f['s2'] = Spectrum(x2,scale*(m(x2)+err),scale*err)
    
    f['o1'] = PhotObservation('BVRI',[13,12.1,11.8,11.5],.153)
    f['o2'] = PhotObservation('ugriz',randn(5,12)+12,rand(5,12)/3)
    
    return f