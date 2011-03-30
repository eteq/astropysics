#!/usr/bin/env python
from __future__ import division,with_statement
from astropysics.constants import pi
import numpy as np
from astropysics.objcat import StructuredFieldNode,Catalog,Field,CycleError,CycleWarning
from nose import tools

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
    
    @StructuredFieldNode.derivedFieldFunc
    def d1(val='val',d2='d2'):
        if d2 is not None:
            return val+np.exp(d2)
    
    @StructuredFieldNode.derivedFieldFunc(type=float) 
    def d2(d1='d1'):
        if d1 is not None: 
            return np.log(d1)


class Test3(StructuredFieldNode):
    a = Field('a',float)
    b = Field('b',float)
    
    @StructuredFieldNode.derivedFieldFunc(units='Msun',b='b')
    #b2 should fail because it doesn't exist, but b='b' above overrides
    def d(a='a',b='b2'): 
        return a - b/2
    
def test_deps():
    """
    Test DependentValue objects

    They must correctly deduce their values and follow the proper rules.
    """
    import warnings
    
    #supress warnings from initialization of cycle-producing objects
    warnings.simplefilter('ignore',CycleWarning)
    
    #check basic dependency results and units
    o2 = Test2(None)
    tools.assert_equal(o2.val(),4.2) #default value
    
    #these should both result in a cycle with each other
    tools.assert_raises(CycleError,o2.d1)
    tools.assert_raises(CycleError,o2.d2)

    #Try setting o2.d2 to various non-float objects
    def settoval(o,val):
        o[None] = val
    tools.assert_raises(TypeError,settoval,o2.d2,'astring')
    tools.assert_raises(TypeError,settoval,o2.d2,1)
    tools.assert_raises(TypeError,settoval,o2.d2,[4.0,2.3])
    tools.assert_raises(TypeError,settoval,o2.d2,np.array(1.5))
        
    o2.d2[None] = 1.5
    #should still fail b/c it is not set to the current value
    tools.assert_raises(CycleError,o2.d1)
    o2['d2'] = None #sets current to None/default
    tools.assert_almost_equal(o2.d1(),8.6816890703380647,12)
    
        
    #make sure derived arguments appear in the proper order - if a nad b are 
    #swapped, this gives the wrong answer.
    o3 = Test3(None)
    o3['a'] = (None,930016275284.81763)
    o3['b'] = (None,3000000000.0)
    #print 'should',o3.a()-o3.b()/2.0,'not',o3.b()-o3.a()/2.0
    tools.assert_equal(o3.d(),o3.a()-o3.b()/2.0)

    return o2,o3

def test_cat():
    """
    Test Catalog objects and related basic actions
    """
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
    """
    Test SEDField
    """
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