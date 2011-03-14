#!/usr/bin/env python
from __future__ import division,with_statement

from nose.tools import assert_almost_equal

def test_gal():
    """Cross-check Gal <-> Supergal <-> FK5 coordinate conversions.
    
    Implicitly also tests networkx conversion routing and matrix composition of
    transforms.
    
    Thanks to David Nicholls for the data set used for comparison.
    
    """
    from astropysics.coords.coordsys import SupergalacticCoordinates,\
                                        GalacticCoordinates,FK5Coordinates
                                   
    
    #data set computed with IDL glactc.pro and cross-checks with catalogs
    #RA,Dec,Glong,Glat,SGlong,SGlat
    s="""
    00:02:46.30,-52:46:18,319.1284,-62.7990,242.7085,-4.8166
    02:06:15.80,-60:56:24,287.5992,-53.9043,236.4422,-22.3149
    04:06:07.90,-52:40:06,261.9954,-45.9695,238.7820,-40.3614
    06:00:10.70,-31:47:14,237.7245,-24.0782,241.8464,-69.6481
    10:01:33.60,-06:31:30,245.9121,36.8999,110.4980,-43.4303
    12:00:47.40,-03:25:12,279.1791,57.0976,116.1007,-13.9687
    14:03:34.60,-27:16:47,322.0616,32.8979,147.5406,7.3568
    16:09:43.90,-00:06:55,11.5871,35.1849,133.7201,46.2550
    20:12:43.20,-03:54:22,38.8727,-19.8409,252.4600,62.5355
    22:07:50.90,-43:16:43,355.9298,-53.3561,240.8982,16.3463
    """.strip()

    fk5s = []
    fk2gals = []
    gals = []
    gal2sgals = []
    sgals = []
    fk2sgals = []

    for l in s.split('\n'):
        ls = l.strip().split(',')
        fk5s.append(FK5Coordinates(ls[0],ls[1],epoch=2000))
        gals.append(GalacticCoordinates(ls[2],ls[3]))
        fk2gals.append(fk5s[-1].convert(GalacticCoordinates))
        sgals.append(SupergalacticCoordinates(ls[4],ls[5]))
        gal2sgals.append(gals[-1].convert(SupergalacticCoordinates))
        fk2sgals.append(fk5s[-1].convert(SupergalacticCoordinates))
    
    for i in range(len(fk5s)):
        assert (gal2sgals[i]-sgals[i]).arcsec < 1,'Gal->SGal not within 1 arcsec:%f'%(gal2sgals[i]-sgals[i]).arcsec
        assert (fk2gals[i]-gals[i]).arcsec < 2,'FK5->Gal not within 2 arcsec:%f'%(fk2gals[i]-gals[i]).arcsec
        assert (fk2sgals[i]-sgals[i]).arcsec < 2,'FK5->SGal not within 2 arcsec:%f'%(fk2sgals[i]-sgals[i]).arcsec
        
    #now reverse the conversions just to make sure everything is symmetric
    
    for i in range(len(fk5s)):
        fksgalfk = (fk2sgals[i].convert(FK5Coordinates)-fk5s[i]).arcsec
        assert fksgalfk < 1e-9,'Fk5->SGal->FK5 too large:%g'%fksgalfk
        
        galsgalgal = (gal2sgals[i].convert(GalacticCoordinates)-gals[i]).arcsec
        assert galsgalgal < 1e-9,'Gal->SGal->Gal too large:%g'%galsgalgal
        
        fkgalfk = (fk2gals[i].convert(FK5Coordinates)-fk5s[i]).arcsec
        assert galsgalgal < 1e-9,'Fk5->Gal->Fk5 too large:%g'%galsgalgal
        
    return fk5s,fk2gals,gals,gal2sgals,sgals,fk2sgals

def test_main_eq_symm(rasdecs=None):
    """
    Test FK4<->FK5<->ICRS<->GCRS coordinate conversions.
    """
    from numpy import mgrid,array
    from astropysics.coords.coordsys import FK4Coordinates,FK5Coordinates, \
                                            ICRSCoordinates,GCRSCoordinates
    
    if rasdecs is None:
        rasdecs = (mgrid[0:360:6j,-80:80:5j]).reshape((2,6*5)).T
    
    gs = [GCRSCoordinates(ra,dec) for ra,dec in rasdecs]
    ics,f5s,f4s,f5s2,ics2,gs2 = [],[],[],[],[],[]
    
    for g in gs:
        ics.append(g.convert(ICRSCoordinates))
        f5s.append(ics[-1].convert(FK5Coordinates))
        f4s.append(f5s[-1].convert(FK4Coordinates))
        f5s2.append(f4s[-1].convert(FK5Coordinates))
        ics2.append(f5s2[-1].convert(ICRSCoordinates))
        gs2.append(ics2[-1].convert(GCRSCoordinates))
    
    gdiffs = []
    idiffs = []
    f5diffs = []
    for i in range(len(gs)):
        gdiff = (gs[i]-gs2[i]).arcsec
        idiff = (ics[i]-ics2[i]).arcsec
        f5diff = (f5s[i]-f5s2[i]).arcsec
        
        assert gdiff< 4e-10,'GCRS<-...->GCRS too large:%g'%gdiff
        assert idiff< 3e-10,'ICRS<-...->ICRS too large:%g'%idiff
        assert f5diff< 2e-10,'FK5<-...->FK5 too large:%g'%f5diff
        
        gdiffs.append(gdiff)
        idiffs.append(idiff)
        f5diffs.append(f5diff)
        
    return array(gdiffs),array(idiffs),array(f5diffs)
           
def test_cirs_eqx_symm(rasdecs=None):
    """
    Test GCRS<->ITRS and intermediate coordinate conversions.
    """
    from numpy import mgrid,array
    from astropysics.coords.coordsys import GCRSCoordinates,CIRSCoordinates, \
                                  EquatorialCoordinatesEquinox,ITRSCoordinates
    
    
    if rasdecs is None:
        rasdecs = (mgrid[0:360:6j,-80:80:5j]).reshape((2,6*5)).T
    
    gs = [GCRSCoordinates(ra,dec) for ra,dec in rasdecs]
    
    
    #through cirs
    cs,tcs,cs2,gs2 = [],[],[],[]
    for g in gs:
        cs.append(g.convert(CIRSCoordinates))
        tcs.append(cs[-1].convert(ITRSCoordinates))
        cs2.append(tcs[-1].convert(CIRSCoordinates))
        gs2.append(cs2[-1].convert(GCRSCoordinates))
    
    for i in range(len(gs)):
        gdiff = (gs2[i]-gs[i]).arcsec
    
    #through equinox
    eqs,tcs2,eqs2,gs3 = [],[],[],[]
    for g in gs:
        eqs.append(g.convert(EquatorialCoordinatesEquinox))
        tcs2.append(eqs[-1].convert(ITRSCoordinates))
        eqs2.append(tcs2[-1].convert(EquatorialCoordinatesEquinox))
        gs3.append(eqs2[-1].convert(GCRSCoordinates))
        
    gds1,gds2,tds,cds,eds = [],[],[],[],[]
    for i in range(len(gs)):
        gdiff1 = (gs2[i]-gs[i]).arcsec
        gdiff2 = (gs3[i]-gs[i]).arcsec
        
        tdiff = (tcs2[i]-tcs[i]).arcsec
        cdiff = (cs2[i]-cs[i]).arcsec
        ediff = (eqs2[i]-eqs[i]).arcsec
        
        assert gdiff1< 5e-10,'GCRS<-..CIRS..->GCRS too large:%g'%gdiff1
        assert cdiff< 5e-10,'CIRS->ITRS->CIRS too large:%g'%cdiff
        
        assert gdiff2< 5e-10,'GCRS<-..Equinox..->GCRS too large:%g'%gdiff2
        assert ediff< 5e-10,'Eq->ITRS->Eq too large:%g'%ediff
        
        #TODO:fix this difference when equinox->ITRS is fixed
        assert tdiff< 60,'GCRS->ITRS between CIRS and Eq too large:%g'%tdiff
        
        gds1.append(gdiff1)
        gds2.append(gdiff2)
        tds.append(tdiff)
        cds.append(cdiff)
        eds.append(ediff)
        
    return array(gds1),array(cds),array(gds2),array(eds),array(tds)

def test_cirs_eqx_ecl(rasdecs=None):
    """
    Test Ecliptic transforms between CIRS and Equinox.
    """
    from numpy import mgrid,array
    from astropysics.coords.coordsys import CIRSCoordinates, \
          EquatorialCoordinatesEquinox,EclipticCoordinatesCIRS,\
          EclipticCoordinatesEquinox,RectangularGeocentricEclipticCoordinates
    
    
    if rasdecs is None:
        rasdecs = (mgrid[0:360:6j,-80:80:5j]).reshape((2,6*5)).T
    
    cs = [CIRSCoordinates(ra,dec) for ra,dec in rasdecs]
    
    ecs,rgs,ecxs,eqxs,ecxs2,rgs2,ecs2,cs2 = [],[],[],[],[],[],[],[]
    for c in cs:
        ecs.append(c.convert(EclipticCoordinatesCIRS))
        rgs.append(ecs[-1].convert(RectangularGeocentricEclipticCoordinates))
        ecxs.append(rgs[-1].convert(EclipticCoordinatesEquinox))
        eqxs.append(ecxs[-1].convert(EquatorialCoordinatesEquinox))
        ecxs2.append(eqxs[-1].convert(EclipticCoordinatesEquinox))
        rgs2.append(ecxs2[-1].convert(RectangularGeocentricEclipticCoordinates))
        ecs2.append(rgs2[-1].convert(EclipticCoordinatesCIRS))
        cs2.append(ecs2[-1].convert(CIRSCoordinates))
        
    cds,ecds,rgds,ecxds = [],[],[],[]
    for i in range(len(cs)):
        cdiff = (cs2[i]-cs[i]).arcsec
        ecdiff = (ecs2[i]-ecs[i]).arcsec
        rgdiff = (rgs2[i]-rgs[i]).length
        ecxdiff = (ecxs2[i]-ecxs[i]).arcsec
        
        assert cdiff< 5e-10,'CIRS->...->CIRS too large:%g'%cdiff
        assert ecdiff< 5e-10,'EcCIRS->...->EcCIRS too large:%g'%ecdiff
        assert rgdiff< 2e-15,'RectEc->...->RectEc too large:%g'%rgdiff
        assert ecxdiff< 5e-10,'Eqx->...->Eqx too large:%g'%ecxdiff
        
        cds.append(cdiff)
        ecds.append(ecdiff)
        rgds.append(rgdiff)
        ecxds.append(ecxdiff)   
        
    return array(cds),array(ecds),array(rgds),array(ecxds)

def test_icrs_rect():
    """
    Test ICRSCoordinates <-> RectangularICRSCoordinates conversions.
    """
    from astropysics.coords.coordsys import RectangularICRSCoordinates,\
                                            ICRSCoordinates
    from numpy import  array,mgrid
    from numpy.random import randn
    from nose.tools import assert_almost_equal
    
#    ntests = 5
#    coords = randn(ntests,3)
    coords = mgrid[-1.5:1.5:5j,-1.5:1.5:5j,-1.5:1.5:5j].reshape((3,5*5*5)).T
    
    
    xds,yds,zds = [],[],[]
    for x,y,z in coords:
        if x==0 and y==0 and z==0:
            continue
        unit = 'au' if x<0 else 'pc'
        r = RectangularICRSCoordinates(x,y,z,unit=unit)
        c = r.convert(ICRSCoordinates)
        r2 = c.convert(RectangularICRSCoordinates)
        c2 = r2.convert(ICRSCoordinates)
        r3 = c2.convert(RectangularICRSCoordinates)
        r3.unit = unit
        
        #unit conversion only good to ~5 places
        assert_almost_equal(x,r3.x,5)
        assert_almost_equal(y,r3.y,5)
        assert_almost_equal(z,r3.z,5)
        
        xds.append(x-r3.x)
        yds.append(y-r3.y)
        zds.append(z-r3.z)
        
    return array(xds),array(yds),array(zds)
        
def test_gcrs_rect():
    """
    Test GCRSCoordinates <-> RectangularGCRSCoordinates conversions.
    """
    from astropysics.coords.coordsys import RectangularGCRSCoordinates,\
                                            GCRSCoordinates
    
    from numpy import  array,mgrid
    from numpy.random import randn
    from nose.tools import assert_almost_equal
    
#    ntests = 5
#    coords = randn(ntests,3)
    coords = mgrid[-1.5:1.5:5j,-1.5:1.5:5j,-1.5:1.5:5j].reshape((3,5*5*5)).T
    
    xds,yds,zds = [],[],[]
    for x,y,z in coords:
        if x==0 and y==0 and z==0:
            continue
        unit = 'au' if x<0 else 'pc'
        r = RectangularGCRSCoordinates(x,y,z,unit=unit)
        c = r.convert(GCRSCoordinates)
        r2 = c.convert(RectangularGCRSCoordinates)
        c2 = r2.convert(GCRSCoordinates)
        r3 = c2.convert(RectangularGCRSCoordinates)
        r3.unit = unit
        
        #unit conversion only good to ~5 places
        assert_almost_equal(x,r3.x,5)
        assert_almost_equal(y,r3.y,5)
        assert_almost_equal(z,r3.z,5)
        
        xds.append(x-r3.x)
        yds.append(y-r3.y)
        zds.append(z-r3.z)
        
    return array(xds),array(yds),array(zds)
        
def test_ecliptic_rect():
    """
    Test RectangularGeocentricEclipticCoordinates -> Ecliptic
    """
    from astropysics.coords.coordsys import EclipticCoordinatesCIRS,\
                                            EclipticCoordinatesEquinox,\
                                    RectangularGeocentricEclipticCoordinates
                                            
     
    from numpy.random import randn
    from nose.tools import assert_almost_equal
    
    ntests = 5
    coords = randn(ntests,3)
    
    for x,y,z in coords:
        r = RectangularGeocentricEclipticCoordinates(x,y,z)
        c1 = r.convert(EclipticCoordinatesCIRS)
        c2 = r.convert(EclipticCoordinatesEquinox)
        r1 = c1.convert(RectangularGeocentricEclipticCoordinates)
        r2 = c2.convert(RectangularGeocentricEclipticCoordinates)
        
        places = 7
        assert_almost_equal(x,r1.x,places)
        assert_almost_equal(y,r1.y,places)
        assert_almost_equal(z,r1.z,places)
        
        assert_almost_equal(x,r2.x,places)
        assert_almost_equal(y,r2.y,places)
        assert_almost_equal(z,r2.z,places)
        
def test_parallax(plot=False,icoords=None):
    from astropysics.coords.coordsys import ICRSCoordinates,GCRSCoordinates, \
                        RectangularGCRSCoordinates, RectangularICRSCoordinates
    from astropysics.constants import asecperrad
    from numpy import linspace,array,radians,mean,max
    
    
    if icoords is None:
        e = 23.439214393375188
        icoords = [(1,1),(270,90-e),(1,30),(180,-30),(80,-89)]
        #icoords = [(1,1),(45,0),(45,-30),(45,-45),(45,-60),(45,-75),(45,-89)]
        #icoords = [(360.1,0),(145,0),(145,30),(145,45),(145,60),(145,75),(145,89)]
        
    distpcs = [1 for i in icoords]
    epochs = linspace(2000,2001,50)
    
    asdiffs = []
    drasall = []
    ddecsall = []
    icsall = []
    gcsall = []
    ricsall = []
    rgcsall = []
    
    for (ra,dec),d in zip(icoords,distpcs):
        
        ics = []
        gcs = []
        rics = []
        rgcs = []
        
        for e in epochs:
            ics.append(ICRSCoordinates(ra,dec,distancepc=d,epoch=e))
            gcs.append(ics[-1].convert(GCRSCoordinates))
            rics.append(ics[-1].convert(RectangularICRSCoordinates))
            rgcs.append(gcs[-1].convert(RectangularGCRSCoordinates))
    
        asdiffs.append([(g-ics[0]).arcsec for g in gcs])
        drasall.append([(g.ra.r-ics[0].ra.r)*asecperrad for g in gcs])
        ddecsall.append([(g.dec.r-ics[0].dec.r)*asecperrad for g in gcs])
        
        icsall.append(ics)
        gcsall.append(gcs)
        ricsall.append(rics)
        rgcsall.append(rgcs)
        
        
    asdiffs = array(asdiffs)
    
    if plot:
        from matplotlib import pyplot as plt
        
        plt.figure(1)
        if plot != 'notclf':
            plt.clf()
            
        for asd,ics in zip(asdiffs,icsall):
            ic = ics[0]
            plt.plot(epochs-2000,asd,label='%.2f,%.2f'%(ic.ra.d,ic.dec.d))
            
        plt.xlabel('epoch - 2000')
        plt.ylabel('$\Delta_{\\rm ICRS,GCRS}$')
        plt.legend(loc=0)
        
        
        plt.figure(2)
        if plot != 'notclf':
            plt.clf()
            
        for dras,ddecs,ics in zip(drasall,ddecsall,icsall):
            ic = ics[0]
            plt.plot(dras,ddecs,label='%.2f,%.2f'%(ic.ra.d,ic.dec.d))
            
        plt.xlabel('$\Delta$RA')
        plt.ylabel('$\Delta$Dec')
        plt.xlim(-3,3)
        plt.ylim(-3,3)
        
        plt.legend(loc=0)
    
    assert max(asdiffs)<=1.05,'Object at 1 pc moves significantly more than 1 arcsec from center:%.3f'%max(asdiffs)
    
    return epochs,asdiffs,icsall,gcsall,ricsall,rgcsall
        