#!/usr/bin/env python
"""
This script runs tests of astropysics functionality that can be tested against
SOFA equivalents or semi-equivalents.  
"""
from __future__ import division

from astropysics import coords,obstools
from math import pi
import os
import numpy as np

tests = []

@tests.append
def trans_matricies(epoch,jd):
    icrs = coords.ICRSCoordinates(0,0,epoch=epoch)

    B = coords.ICRSCoordinates.frameBiasJ2000 
    P = coords.coordsys._precession_matrix_J2000_Capitaine(epoch)
    N = coords.coordsys._nutation_matrix(epoch)
    print 'B',B #disagrees in 6thish decimal?
    print 'P',P #disagrees in 4thish decimal?
    print 'N',N #good

    print 'BPN\n',icrs.getConversionMatrix(coords.EquatorialCoordinatesEquinox)
    print 'C\n',icrs.getConversionMatrix(coords.EquatorialCoordinatesCIRS)
    
@tests.append
def trans_coords(epoch,jd):
    icrs = coords.ICRSCoordinates('10d','20d',epoch=epoch)
    print icrs.ra.radians,icrs.dec.radians
    print icrs,'->',icrs.convert(coords.EquatorialCoordinatesCIRS)
    print icrs,'->',icrs.convert(coords.EquatorialCoordinatesEquinox)
    #print icrs,'->',icrs.convert(coords.ITRSCoordinates)
    
def compile_c():
    if not os.path.exists('sofatests'):
        res = os.system('./sofatests.build')
        if res != 0:
            raise ValueError('Could not build C sofatests.  Is SOFA installed, and are the paths in sofatests.build correct?')
        
def run_test(testname,testfunc,args):
    print '\nRunning PYTHON test',testname
    testfunc(*args)

    cmd = './sofatests %s %.16g %.16g'%(testname,epoch,jd)
    print 'Running C test',testname,'as',cmd  
    os.system(cmd)


tests = [(tfunc.func_name,tfunc) for tfunc in tests]
if __name__ == '__main__':
    from optparse import OptionParser
    
    op = OptionParser()
    op.add_option('-e','--epoch',default=None,help="the epoch to assume")
    op.add_option('-c','--compile',action='store_true',default=False,help="compile C tests even if they're already present")
    op.usage = '%prog [options] [testname1 testname2 ...]'
    ops,args = op.parse_args()
    
    if not os.path.exists('sofatests') or ops.compile:
        print 'Compiling:','./sofatests.build'
        res = os.system('./sofatests.build')
        if res != 0:
            raise ValueError('Could not build C sofatests.  Is SOFA installed, and are the paths in sofatests.build correct?')
    
    if ops.epoch is None:
        epoch = obstools.jd_to_epoch(None)
        jd = obstools.epoch_to_jd(epoch)
        print 'No Epoch given - current:',epoch,'JD:',jd
    else:
        epoch = float(ops.epoch)
        jd = obstools.epoch_to_jd(ops.epoch)
    
    if len(args)==0:
        print 'No Tests specified - running all.'
        for tn,tf in tests:
            run_test(tn,tf,(epoch,jd))
        print 'Done running tests!'
    else:
        testd = dict(tests)
        for arg in args:
            try:
                testd[arg]()
            except KeyError:
                print "Couldn't find python test %s, skipping"%arg

