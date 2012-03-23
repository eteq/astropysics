import astropysics.coords
import astropysics.obstools
import datetime,pytz


def greenwich():
    return astropysics.obstools.Site(lat=51.5,
                                    long=0,
                                    alt=0,
                                    tz=0,
                                    name="Greenwich"
                                    )
    
vernal_equinox_2012 = datetime.datetime(2012, 03, 20,
                                           5, 14,
                                           tzinfo = pytz.utc)

equatorial_on_sky_ve = astropysics.coords.coordsys.FK5Coordinates("17:6:32 +0.0 J2000.0")

def test_rst_functionality():
    site = greenwich()
    r,s,t = site.riseSetTransit(equatorial_on_sky_ve,
                                vernal_equinox_2012,
                                timeobj=True,utc=True
                                )
    
    equinox_day = vernal_equinox_2012.date()
    print 
    print equinox_day, t.date()
    assert t.date() == equinox_day



