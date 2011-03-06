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


#HELPERS
def _cartesian_transform(coord):
    from math import sin,cos

    lat = coord.lat.radians
    lng = coord.long.radians    
    
    sb = sin(lat)
    cb = cos(lat)
    sl = sin(lng)
    cl = cos(lng)
    
    #spherical w/ r=1 > cartesian
    x = cb*cl
    y = cb*sl
    z = sb
    return x,y,z

def _computematrix(intype,outtype,epoch):
    #(x,y,z)=(1,0,0)
    incx = intype(epoch=epoch)
    incx.lat.d = 0
    incx.long.d = 0
    outcx = incx.convert(outtype)
    
    #(x,y,z)=(0,1,0)
    incy = intype(epoch=epoch)
    incy.lat.d = 0
    incy.long.d = 90
    outcy = incy.convert(outtype)
    
    #(x,y,z)=(0,0,1)
    incz = intype(epoch=epoch)
    incz.lat.d = 90
    incz.long.d = 0
    outcz = incz.convert(outtype)
    
    Mt = (_cartesian_transform(outcx),_cartesian_transform(outcy),_cartesian_transform(outcz))
    return np.matrix(Mt).T
    
#END HELPERS

tests = []
   
@tests.append
def earth_rotation(epoch,jd):
    ERA = coords.earth_rotation_angle(jd,False)
    print 'ERA',ERA
    GAST = coords.greenwich_sidereal_time(jd,'simple')*pi/12
    print 'simple GAST (change to non-simple when ready)',GAST
    GMST = coords.greenwich_sidereal_time(jd,False)*pi/12
    print 'GMST', GMST
    
    def cmpfunc(cstdout,ERA,GAST,GMST):
        from warnings import warn
        
        ls = [l.split(':')[-1].strip() for l in cstdout.split('\n') if l.strip()!='']
        dERA = ERA-float(ls[0])
        dGAST = GAST-float(ls[1])
        dGMST = GMST-float(ls[2])
        
        for s in ('dERA','dGAST','dGMST'):
            val = locals()[s]
            if abs(val)>1e-8:
                warn('%s too large:%g !!!'%(s,val))
        
    return cmpfunc,ERA,GAST,GMST
    
@tests.append
def trans_matricies(epoch,jd):
    icrs = coords.ICRSCoordinates(0,0,epoch=epoch)

    B = coords.ICRSCoordinates.frameBiasJ2000 
    P = coords.coordsys._precession_matrix_J2000_Capitaine(epoch)
    N = coords.coordsys._nutation_matrix(epoch)
    print 'B',B #disagrees in 6thish decimal?
    print 'P',P #disagrees in 4thish decimal?
    print 'N',N #good
    
    BPN = _computematrix(coords.ICRSCoordinates,coords.EquatorialCoordinatesEquinox,epoch)
    print 'BPN\n',BPN
    C = _computematrix(coords.ICRSCoordinates,coords.CIRSCoordinates,epoch)
    print 'C\n',C
    
    def cmpfunc(cstdout,B,P,N,BPN,C):
        from warnings import warn
        
        ls = [l for l in cstdout.split('\n') if l.strip()!='']
        cB = np.matrix(eval(','.join(ls[1:4])))
        cP = np.matrix(eval(','.join(ls[5:8])))
        cN = np.matrix(eval(','.join(ls[9:12])))
        cBPN = np.matrix(eval(','.join(ls[13:16])))
        cC = np.matrix(eval(','.join(ls[17:20])))
        
        dB = np.abs(cB-B).max()
        dP = np.abs(cP-P).max()
        dN = np.abs(cN-N).max()
        dBPN = np.abs(cBPN-BPN).max()
        dC = np.abs(cC-C).max()
        
        for s in ('dB','dP','dN','dBPN','dC'):
            val = locals()[s]
            if val>5e-7:
                warn('%s too large:%g !!!'%(s,val))
    
    return cmpfunc,B,P,N,BPN,C
    
@tests.append
def trans_coords(epoch,jd):
    icrs = coords.ICRSCoordinates('10d','20d',epoch=epoch)
    print icrs
    #print "cartesian",_cartesian_transform(icrs)
    c = icrs.convert(coords.CIRSCoordinates)
    e = icrs.convert(coords.EquatorialCoordinatesEquinox)
    i = icrs.convert(coords.ITRSCoordinates)
    print 'ICRS->',c
    print 'ICRS->',e
    print 'ICRS->',i
    
    def cmpfunc(cstdout,c,e,i):
        from warnings import warn
        
        ldats = [l[(l.index(':')+2):].strip() for l in cstdout.split('\n') if l.strip()!='']
        racs,deccs = (ldats[0].split()[0])[3:],(ldats[0].split()[1])[4:]
        raes,deces = (ldats[1].split()[0])[3:],(ldats[1].split()[1])[4:]
        ilat,ilong = (ldats[2].split()[0])[4:],(ldats[2].split()[1])[5:]
        
        drac = coords.AngularCoordinate(racs).d-c.ra.d
        ddecc = coords.AngularCoordinate(deccs).d-c.dec.d
        drae = coords.AngularCoordinate(raes).d-e.ra.d
        ddece = coords.AngularCoordinate(deces).d-e.dec.d
        dlat = float(ilat)-i.lat.d
        dlong = float(ilong)-i.long.d
        
        for s in ('drac','ddecc','drae','ddece','dlat','dlong'):
            val = locals()[s]
            if abs(val)>5e-8:
                warn('%s too large:%g !!!'%(s,val))
        
        
    return cmpfunc,c,e,i
    
@tests.append
def earth_pv(epoch,jd):
    p,v = coords.ephems.earth_pos_vel(jd,barycentric=False,kms=False)
    pb,vb = coords.ephems.earth_pos_vel(jd,barycentric=True,kms=False)
    
    print 'Heliocentric pos:',p
    print 'Heliocentric vel:',v/365.25
    print 'SS Barycentric pos:',pb
    print 'SS Barycentric vel:',vb/365.25
    
    def cmpfunc(cstdout,p,v,pb,vb):
        from warnings import warn
        
        ldats = [l.split(':')[-1].strip() for l in cstdout.split('\n') if l.strip()!='']
        pc,vc,pbc,vbc = [np.array(l.split(','),dtype=float) for l in ldats]
        
        if np.abs(p-pc).max() > 1e-6:
            warn('earth_pv C and python do not match: p-pc=%g !!!'%(p-pc))
        if np.abs(v-vc).max() > 1e-6:
            warn('earth_pv C and python do not match: v-vc=%g !!!'%(v-vc))
        if np.abs(pb-pbc).max() > 1e-6:
            warn('earth_pv C and python do not match: pb-pbc=%g !!!'%(pb-pbc))
        if np.abs(vb-vbc).max() > 1e-6:
            warn('earth_pv C and python do not match: vb-vbc=%g !!!'%(vb-vbc))
             
    return cmpfunc,p,v/365.25,pb,vb/365.25
    
def compile_c():
    if not os.path.exists('sofatests'):
        res = os.system('./sofatests.build')
        if res != 0:
            raise ValueError('Could not build C sofatests.  Is SOFA installed, and are the paths in sofatests.build correct?')
        
def run_test(testname,testfunc,epoch,jd):
    """
    Run the python and C tests with the given name and python function - the
    name will be used to infer the C test.
    
    if the python functions return anything other than None, it will be assumed
    that the first return value is a function of the form f(cstdoput,*pyreturn).
    """
    import subprocess
    
    print '\n\n\nRunning PYTHON test',testname
    pyres = testfunc(epoch,jd)
    

    cmd = './sofatests %s %.16g %.16g'%(testname,epoch,jd)
    print '\nRunning C test',testname,'as','"'+cmd+'"'  
    
    if pyres is None:
        os.system(cmd)
    else:
        cmpfunc = pyres[0]
        
        sp = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        cstdout,cstderr = sp.communicate()
        if sp.returncode==0:
            print 'C code complete - output:\n',cstdout
        else:
            print 'Problem with C - retcode',sp.returncode,'- output:\n',cstdout,'\nstd err:',cstderr
        
        cmps = []
        cmps.append(cstdout)
        cmps.extend(pyres[1:])
        print 'Comparing python to C'
        cmpfunc(*cmps)
        print 'Comparison complete'


tests = [(tfunc.func_name,tfunc) for tfunc in tests]

def main(compile=False,epoch=None,teststorun=None):
    if not os.path.exists('sofatests') or compile:
        print 'Compiling:','./sofatests.build'
        res = os.system('./sofatests.build')
        if res != 0:
            raise ValueError('Could not build C sofatests.  Is SOFA installed, and are the paths in sofatests.build correct?')
    
    if epoch is None:
        epoch = obstools.jd_to_epoch(None)
        jd = obstools.epoch_to_jd(epoch)
        print 'No Epoch given - current:',epoch,'JD:',jd
    else:
        epoch = float(ops.epoch)
        jd = obstools.epoch_to_jd(ops.epoch)
    
    if teststorun is None:
        print 'No Tests specified - running all.'
        for tn,tf in tests:
            run_test(tn,tf,epoch,jd)
        print 'Done running tests!'
    else:
        testd = dict(tests)
        for arg in args:
            try:
                run_test(arg,testd[arg],epoch,jd)
            except KeyError:
                print "Couldn't find python test %s, skipping"%arg

if __name__ == '__main__':
    from optparse import OptionParser
    
    op = OptionParser()
    op.add_option('-e','--epoch',default=None,help="the epoch to assume")
    op.add_option('-c','--compile',action='store_true',default=False,help="compile C tests even if they're already present")
    op.usage = '%prog [options] [testname1 testname2 ...]'
    ops,args = op.parse_args()
    
    if len(args)==0:
        args = None
    main(ops.compile,ops.epoch,args)
