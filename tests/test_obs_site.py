## Tim Staley , 2012

import astropysics.coords
import astropysics.obstools
import datetime, pytz
import unittest

from astropysics.coords.coordsys import FK5Coordinates

#---------------------------------------------------------------------------
#Test data: 

def greenwich():
    return astropysics.obstools.Site(lat=51.5,
                                    long=0,
                                    alt=0,
                                    tz=0,
                                    name="Greenwich"
                                    )

vernal_equinox_2012 = datetime.datetime(2012, 03, 20,
                                           5, 14,
                                           tzinfo=pytz.utc)

vernal_equinox_2012_p12h = vernal_equinox_2012 + datetime.timedelta(hours=12)


#Greenwich LST at V.E. (h:m:s) = 17:6:32  approx
greenwish_lst_at_ve = 17.1096992 #decimal hours
#12 hours later:
greenwish_lst_at_ve_p12h = 5.1425530145722549

##NB this is approximately zenith at the 2012 vernal equinox (at greenwich):
equatorial_on_sky_at_ve = FK5Coordinates("17:6:32 +0.0 J2000.0")

equatorial_off_sky_at_ve = FK5Coordinates("5:6:32 +0.0 J2000.0")

equatorial_sets_day_after_ve = FK5Coordinates("6:6:32 +0.0 J2000.0")

circumpolar_north_transit_at_ve = FK5Coordinates("17:6:32 +89.0 J2000.0")
circumpolar_north_transit_later = FK5Coordinates("5:6:32 +89.0 J2000.0")

never_visible_source = FK5Coordinates("17:6:32 -70.0 J2000.0")
#---------------------------------------------------------------------------


class TestLocalToSiderealConversions(unittest.TestCase):
    def setUp(self):
        self.site = greenwich()

    def test_local_sidereal_time(self):
        computed_lst_at_ve = self.site.localSiderialTime(vernal_equinox_2012)
        self.assertAlmostEqual(computed_lst_at_ve, greenwish_lst_at_ve)

        computed_lst_at_ve_p12 = self.site.localSiderialTime(vernal_equinox_2012 +
                                            datetime.timedelta(hours=12))
        self.assertAlmostEqual(computed_lst_at_ve_p12 , greenwish_lst_at_ve_p12h)

    def test_local_time_at_ve(self):
        computed_local_time_at_ve = self.site.localTime(greenwish_lst_at_ve,
                                                vernal_equinox_2012.date(),
                                                returntype='datetime',
                                                utc=True)
        self.assertEqual(computed_local_time_at_ve.date(),
                         vernal_equinox_2012.date())

        ct = computed_local_time_at_ve
        ve = vernal_equinox_2012
        delta = max(ct, ve) - min(ct, ve)
        self.assertTrue(abs(delta.seconds) < 30)

    def test_local_time_at_ve_p12h(self):
        computed_local_time_at_ve_p12 = self.site.localTime(greenwish_lst_at_ve_p12h,
                                                vernal_equinox_2012.date(),
                                                returntype='datetime',
                                                utc=True)
        self.assertEqual(computed_local_time_at_ve_p12.date(),
                         vernal_equinox_2012.date())

        ct = computed_local_time_at_ve_p12
        ve_p12 = vernal_equinox_2012_p12h
        delta = max(ct, ve_p12) - min(ct, ve_p12)
        self.assertTrue(abs(delta.seconds) < 30)


class TestRSTFunctionality(unittest.TestCase):
    def setUp(self):
        self.site = greenwich()

    def test_basic_rst_for_on_sky_equatorial_source(self):

        r, s, t = self.site.riseSetTransit(eqpos=equatorial_on_sky_at_ve,
                                          date=vernal_equinox_2012.date(),
                                          timeobj=True,
                                          utc=True
                                          )
        for result in r, s, t:
            self.assertIsInstance(result, datetime.datetime)

        #Always True (by specification):
        self.assertEqual(t.date(), vernal_equinox_2012.date())

        #True in this case:
        self.assertTrue(r.date() < vernal_equinox_2012.date())
        self.assertEqual(s.date(), vernal_equinox_2012.date())

        ve = vernal_equinox_2012
        delta = max(t, ve) - min(t, ve)
        self.assertTrue(delta.seconds < 30)

    def test_basic_rst_for_off_sky_equatorial_source(self):

        r, s, t = self.site.riseSetTransit(eqpos=equatorial_off_sky_at_ve,
                                          date=vernal_equinox_2012.date(),
                                          timeobj=True,
                                          utc=True
                                          )
        for result in r, s, t:
            self.assertIsInstance(result, datetime.datetime)

        #Always True (by specification):
        self.assertEqual(t.date(), vernal_equinox_2012.date())

        #True in this case:
        self.assertEqual(r.date() , vernal_equinox_2012.date())
        self.assertEqual(s.date(), vernal_equinox_2012.date())
#
        ve_p12 = vernal_equinox_2012_p12h
        delta = max(t, ve_p12) - min(t, ve_p12)
        self.assertTrue(delta.seconds < 140)

    def test_rst_for_source_setting_tomorrow(self):

        r, s, t = self.site.riseSetTransit(eqpos=equatorial_sets_day_after_ve,
                                          date=vernal_equinox_2012.date(),
                                          timeobj=True,
                                          utc=True
                                          )
        for result in r, s, t:
            self.assertIsInstance(result, datetime.datetime)

        #Always True (by specification):
        self.assertEqual(t.date(), vernal_equinox_2012.date())

        #True in this case:
        self.assertEqual(r.date() , vernal_equinox_2012.date())
        self.assertGreater(s.date(), vernal_equinox_2012.date())


    def regression_test_rst_for_equatorial_sources(self):
        #NB these were simply generated by desk-checking current 
        #(sensible looking) output.
        #So it's currently just a regression test.
        #TO DO: Implement a cross-check against pyephem / catalog results.

        previous_sensible_results = [
             ['2012-03-19 23:11:18.175698+00:00', # Rose yesterday evening 
              '2012-03-20 11:16:36.004857+00:00',
              '2012-03-20 05:13:57.090277+00:00'], #Transit 5am VE day
             ['2012-03-20 11:09:20.223480+00:00',
              '2012-03-20 23:14:38.052638+00:00',
              '2012-03-20 17:11:59.138059+00:00'], #Transit ~5pm VE day
             ['2012-03-20 12:09:10.394128+00:00',
              '2012-03-21 00:14:28.223287+00:00', #Sets day after VE
              '2012-03-20 18:11:49.308707+00:00']
             ]

        test_results = []
        for pos in (equatorial_on_sky_at_ve,
                     equatorial_off_sky_at_ve,
                     equatorial_sets_day_after_ve) :
            r, s, t = self.site.riseSetTransit(eqpos=pos,
                                               date=vernal_equinox_2012.date(),
                                               timeobj=True, utc=True
                                               )
            test_results.append([str(r), str(s), str(t)])
#        print
#        print test_results
#        print
        self.assertEqual(previous_sensible_results, test_results)

    def test_rst_for_circumpolar_source(self):
        for t in vernal_equinox_2012, vernal_equinox_2012_p12h:
            r, s, t = self.site.riseSetTransit(circumpolar_north_transit_at_ve,
                                        t.date(),
                                        timeobj=True, utc=True
                                        )

            equinox_day = vernal_equinox_2012.date()
            self.assertEqual(t.date(), equinox_day)
            self.assertEqual([r, s], [None, None])


    def test_rst_for_never_visible_source(self):
        for t in vernal_equinox_2012, vernal_equinox_2012_p12h:
            r, s, t = self.site.riseSetTransit(eqpos=never_visible_source,
                                                  date=t.date(),
                                                  timeobj=True, utc=True
                                                  )
            self.assertEqual([r, s, t], [None, None, None])





