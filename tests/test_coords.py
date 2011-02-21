#!/usr/bin/env python
from __future__ import division,with_statement

from nose.tools import assert_almost_equal

def test_gal():
    """
    coords:Conversion to Galactic and Supergalactic coordinates. 
    
    Implicitly also tests networkx conversion routing and matrix composition of
    transforms.
    
    Thanks to David Nicholls for the data set.
    
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
        assert (gal2sgals[i]-sgals[i]).arcsec < 1,'Gal->SGal not within 1 arcsec'
        assert (fk2gals[i]-gals[i]).arcsec < 2,'FK5->Gal not within 2 arcsec'
        assert (fk2sgals[i]-sgals[i]).arcsec < 2,'FK5->SGal not within 2 arcsec'
        
    return fk5s,fk2gals,gals,gal2sgals,sgals,fk2sgals

def test_icrs_rect():
    """
    Tests conversion ICRSCoordinates <-> RectangularICRSCoordinates
    """
    from astropysics.coords.coordsys import RectangularICRSCoordinates,\
                                            ICRSCoordinates
     
    from numpy.random import randn
    from nose.tools import assert_almost_equal
    
    ntests = 5
    coords = randn(ntests,3)
    
    for x,y,z in coords:
        unit = 'au' if x<0 else 'pc'
        r = RectangularICRSCoordinates(x,y,z,unit=unit)
        c = r.convert(ICRSCoordinates)
        r2 = c.convert(RectangularICRSCoordinates)
        r2.unit = unit
        
        #unit conversion only good to ~5 places
        assert_almost_equal(x,r2.x,5)
        assert_almost_equal(y,r2.y,5)
        assert_almost_equal(z,r2.z,5)
        
def test_ecliptic_rect():
    """
    Tests conversion ecliptic <-> RectangularGeocentricEclipticCoordinates
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
        