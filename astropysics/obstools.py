#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""

========
obstools
========

The :mod:`obstools` module stores tools for oberving (pre- and post-) as well as
functioning as a module for various corrections and simple calculations that
don't have a better place to live.

Most implementations are for optical astronomy, as that is what the primary
author does.

Note that throughout this module it is assumed that UTC == UT1, so if leap
seconds are abandoned, this will need to be reworked...

Also note that some of these functions require the
:mod:`dateutil <http://pypi.python.org/pypi/python-dateutil>` package (it is
included with matplotlib)

.. todo:: examples/tutorials


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.obstools
   :parts: 1

Module API
----------

"""
#TODO: exposure time calculator (maybe in phot instead?)
#TODO: make Extinction classes spec.HasSpecUnits and follow models framework?
#TODO: return Pipeline support to Extinction 
#TODO: Instruments and telescopes for sites

from __future__ import division,with_statement
from .constants import pi
from .utils import PipelineElement,DataObjectRegistry
import numpy as np

try: 
    from dateutil import tz as tzmod
    tzoffset = tzmod.tzoffset
except ImportError:
    import datetime
    tzmod = None
    class tzoffset(datetime.tzinfo):
        """
        Backup class to do basic fixed-offset time zones if :mod:`dateutil.tz` 
        is missing.
        """
        
        def __init__(self,name,offset):
            if not isinstance(name,basestring):
                raise TypeError('name must be a string')
            self._name = name
            
            self._hoffset = offset
            
        def dst(self):
            return False
        
        def tzname(self):
            return self._name
        
        def utcoffset(self):
            from datetime import timedelta
            return timedelta(hours=self._hoffset)

def jd_to_calendar(jd,rounding=1000000,output='datetime',gregorian=None,mjd=False):
    """
    Convert a julian date to a calendar date and time.
    
    The julian date should be passed as the `jd` parameter, or None to get a 
    
    
    `rounding` determines a fix for floating-point errors. It specifies the
    number of milliseconds by which to round the result to the nearest second.
    If 1000000, no milliseconds are recorded. If larger, a ValueError is raised.
    
    
    
    `output` determines the format of the returned object and can be:
    
    * 'datetime'
        A list of :class:`datetime.datetime` objects in UTC will be returned. If
        the input is a scalar, a single object will be returned.
    * 'array'
        A Nx7 array will be returned of the form
        [(year,month,day,hr,min,sec,msec),...] unless the input was a scalar, in
        which case it will be a length-7 array.
    * 'fracarray'
        An Nx3 array (year,month,day) where day includes the decimal portion.
        
    If `gregorian` is True, the output will be in the Gregorian calendar.
    Otherwise, it will be Julian. If None, it will be assumed to switch over on
    October 4/15 1582.
    
    If `mjd` is True, the input is interpreted as a modified julian date instead
    of a standard julian date.
    
    
    **Examples**
    
    .. testsetup::
    
        from astropysics.obstools import jd_to_calendar
    
    .. doctest::
        
        >>> jd_to_calendar(2451545)
        datetime.datetime(2000, 1, 1, 12, 0, tzinfo=tzutc())
        >>> jd_to_calendar(2305812.5)
        datetime.datetime(1600, 12, 31, 0, 0, tzinfo=tzutc())
        >>> jd_to_calendar([2415020.5,2305447.5],output='array')
        array([[1900,    1,    1,    0,    0,    0,    0],
               [1600,    1,    1,    0,    0,    0,    0]])
        >>> jd_to_calendar(0.0,output='fracarray')
        array([[ -4.71200000e+03,   1.00000000e+00,   1.50000000e+00]])
        
    """
    import datetime
    from dateutil import tz
    
    if jd is None:
        jd = calendar_to_jd(datetime.datetime.now(tz.tzlocal()))
    
    jd = np.array(jd,copy=True,dtype=float)
    scalar = jd.shape == ()
    jd = jd.ravel()
    
    if mjd:
        jd += 2400000.5
    
    if rounding > 1000000:
        raise ValueError('rounding cannot exceed a second')
    elif rounding <= 0:
        jd += .5 
    else:
        rounding = int(rounding)
        roundingfrac = rounding/86400000000
        jd += .5 + roundingfrac 
        
    z = np.floor(jd).astype(int) 
    dec = jd - z #fractional piece
    
    #fix slight floating-point errors if they hapepn TOOD:check
    dgtr1 = dec>=1.0
    dec[dgtr1] -= 1.0
    z[dgtr1] += 1
    
    
    if gregorian is None:
        gregorian = 2299161
        
    if gregorian is True:
        alpha = ((z-1867216.25)/36524.25).astype(int)
        z += 1 + alpha - alpha//4
    elif gregorian is False:
        pass
    else:
        gmask = z >= gregorian
        alpha = ((z[gmask]-1867216.25)/36524.25).astype(int)
        z[gmask] += 1 + alpha - alpha//4
    
    b = z + 1524
    c = ((b-122.1)/365.25).astype(int)
    d = (365.25*c).astype(int)
    e = ((b-d)/30.6001).astype(int)
    
    day = b - d - (30.6001*e).astype(int)
    
    mmask = e<14
    month = e
    month[mmask] -= 1
    month[~mmask] -= 13
    year = c
    year[month>2] -= 4716
    year[month<=2] -= 4715
    
    if output == 'fracarray':
        dec = dec-roundingfrac
        dec[dec<0]=0
        return np.array((year,month,day+dec)).T
    
    if rounding == 1000000:
        secdec = dec*86400
        sec = secdec.astype(int)
        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr
        #sec[sec==secdec] -= 1
        msec = None
    else:
        msec = (dec*86400000000.).astype('int64') 
        if rounding > 0:
            div = (msec//1000000)*1000000
            toround = (msec - div)<(2*rounding)
            msec[toround] = div + rounding
            msec  -= rounding

        sec = msec//1000000
        msec -= 1000000*sec
        
        min = sec//60
        sec -= 60*min
        hr = min//60
        min -= 60*hr
    
    if output == 'datetime':
        tzi = tz.tzutc()
        if msec is None:
            ts = (year,month,day,hr%24,min%60,sec%60)
        else:
            ts = (year,month,day,hr%24,min%60,sec%60,msec%1000000)
        res = [datetime.datetime(*t,**dict(tzinfo=tzi)) for t in zip(*ts)]
    elif output == 'array':
        msec = np.zeros_like(sec) if msec is None else msec
        res = np.array([year,month,day,hr%24,min%60,sec%60,msec]).T
    else:
        raise ValueError('invlid output form '+str(output))
    if scalar:
        return res[0]
    else:
        return res
    
    
    

def calendar_to_jd(caltime,tz=None,gregorian=True,mjd=False):
    """
    Convert calendar value to julian date
    
    the input `caltime` can either be:
    
    * a sequence (yr,month,day,[hr,min,sec]). 
    * a sequence as above with one or more elements a sequence; a sequence will
      be returned.
    * a :class:`datetime.datetime` or :class:`datetime.date` object 
    * a sequence of such objects; a sequence will be returned.
    
    If datetime objects are given, and `tz` is None, values are converted to
    UTC based on the datetime.tzinfo objects (if present).  if `tz` is not None,
    the tzinfo in the datetime objects is ignored.
    
    If the time is unspecified, it is taken to be noon (i.e. Julian Date =
    Julian Day Number)
    
    If `tz` is a string, it is taken to be a timezone name that will be used
    to convert all the dates to UTC (requires :mod:`dateutil` package).  If `tz`
    it is a scalar, it is taken to be the hour offset of the timezone. 
    Or, if it is a :class:`datetime.tzinfo` object, that object will be used to
    do the conversion to UTC.  
    
    If `gregorian` is True, the input will be interpreted as in the Gregorian
    calendar. Otherwise, it will be Julian. If None, it will be assumed to
    switch over on October 4/15 1582.
    
    If `mjd` is True, a modified julian date is returned instead of a standard
    julian date.
    
    
    **Examples**
    
    .. testsetup::
    
        from astropysics.obstools import calendar_to_jd
        import datetime,dateutil
    
    .. doctest::
        
        >>> calendar_to_jd((2010,1,1))
        2455198.0
        >>> calendar_to_jd(datetime.datetime(2000,12,21,3,0,0))
        2451899.625
        >>> calendar_to_jd([2004,3,(5,6)])
        array([ 2453070.,  2453071.])
        >>> dates = [datetime.datetime(2004,3,5),datetime.datetime(2004,3,9)]
        >>> calendar_to_jd(dates)
        array([ 2453069.5,  2453073.5])
        >>> tz = dateutil.tz.tzoffset('2',3*3600)
        >>> calendar_to_jd((2010,1,1),tz)
        2455197.875

        
    """
    #Adapted from xidl  jdcnv.pro
    from datetime import datetime,date,tzinfo
    
    if caltime is None:
        from dateutil.tz import tzlocal
        datetimes = [datetime.now(tzlocal())]
        scalarout = True
    elif isinstance(caltime,datetime) or isinstance(caltime,date):
        datetimes = [caltime]
        scalarout = True
    elif all([isinstance(ct,datetime) or isinstance(ct,date) for ct in caltime]):
        datetimes = caltime
        scalarout = False
    else:
        datetimes = None
        caltime = list(caltime)
        if not (3 <= len(caltime) < 8):
            raise ValueError('caltime input sequence is invalid size')
        while len(caltime) < 7:
            if len(caltime) == 3:
                #make hours 12
                caltime.append(12*np.ones_like(caltime[-1]))
            else:
                caltime.append(np.zeros_like(caltime[-1]))
        yr,month,day,hr,min,sec,msec = caltime
        scalarout = all([np.shape(v) is tuple() for v in caltime])
        
    #if input objects are datetime objects, generate arrays
    if datetimes is not None:
        yr,month,day,hr,min,sec,msec = [],[],[],[],[],[],[]
        for dt in datetimes:
            if not hasattr(dt,'hour'):
                dt = datetime(dt.year,dt.month,dt.day,12)
            
            if tz is None:
                off = dt.utcoffset()
                if off is not None:
                    dt = dt - off
                
            yr.append(dt.year)
            month.append(dt.month)
            day.append(dt.day)
            hr.append(dt.hour)
            min.append(dt.minute)
            sec.append(dt.second)
            msec.append(dt.microsecond)
                
                
    
    yr = np.array(yr,dtype='int64',copy=False).ravel()
    month = np.array(month,dtype='int64',copy=False).ravel()
    day = np.array(day,dtype='int64',copy=False).ravel()
    hr = np.array(hr,dtype=float,copy=False).ravel()
    min = np.array(min,dtype=float,copy=False).ravel()
    sec = np.array(sec,dtype=float,copy=False).ravel()
    msec = np.array(msec,dtype=float,copy=False).ravel()
    
    #do tz conversion if tz is provided  
    if isinstance(tz,basestring) or isinstance(tz,tzinfo):
        if isinstance(tz,basestring):
            from dateutil import tz
            tzi = tz.gettz(tz)
        else:
            tzi = tz
        
        utcoffset = []
        for t in zip(yr,month,day,hr,min,sec,msec):
            #microsecond from float component of seconds
            
            dt = datetime(*[int(ti) for ti in t],**dict(tzinfo=tzi))
            utcdt = dt.utcoffset()
            if utcdt is None:
                utcoffset.append(0)
            else:
                utcoffset.append(utcdt.days*24 + (utcdt.seconds + utcdt.microseconds*1e-6)/3600)
    else:
        utcoffset = tz
            
#    ly = ((month-14)/12).astype(int) #In leap years, -1 for Jan, Feb, else 0
#    jdn = day - 32075l + 1461l*(yr+4800l+ly)//4
    
#    jdn += 367l*(month - 2-ly*12)//12 - 3*((yr+4900l+ly)//100)//4
    
#    res = jdn + (hr/24.0) + min/1440.0 + sec/86400.0 - 0.5

    #this algorithm from meeus 2ed
    m3 = month < 3
    yr[m3] -= 1
    month[m3] += 12
        
    cen = yr//100
    
    if gregorian is None:
        gregorian = (1582,10,4)
    if gregorian is True:
        gregoffset = 2 - cen + cen//4
    elif gregorian is False:
        gregoffset = 0
    else:
        gregoffset = 2 - cen + cen//4
        gmask = (yr>gregorian[0])&(month>gregorian[1])&(day>gregorian[2])
        gregoffset[~gmask] = 0
    
        
    jdn = (365.25*(yr+4716)).astype(int) + \
          (30.6001*(month + 1)).astype(int) + \
               day + gregoffset - 1524.5
    res = jdn + hr/24.0 + min/1440.0 + sec/86400.0
    
    if mjd:
        res -= 2400000.5
    
    if np.any(utcoffset):
        res -= np.array(utcoffset)/24.0
    
    if scalarout:
        return res[0]
    else:
        return res
    

def jd_to_epoch(jd,julian=True,asstring=False):
    """
    Convert a Julian Date to a Julian or Besselian Epoch.
    
    If `julian` is True, a Julian Epoch will be used (the year is exactly 365.25
    days long). Otherwise, the epoch will be Besselian (assuming a tropical year
    of 365.242198781 days). 
    
    if `asstring` is True, a string of the form 'J2000.0' will be returned. If
    it is an integer, the number sets the number of significant figures in the
    output string Otherwise, a scalar is returned (an int if a whole year, float
    if not).
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epj.html
    """
    if julian:
        epoch = 2000.0 + (jd - 2451545.0)/365.25
    else:
        epoch = 1900 + (jd - 2415020.31352)/365.242198781
        
    if round(epoch)==epoch:
        epoch = int(epoch) 
        
    if asstring:
        if asstring is not True:
            fmt = ('J' if julian else 'B')+'{0:.'+str(int(asstring))+'}'
        else: 
            fmt = 'J{0}' if julian else 'B{0}'
        return fmt.format(epoch)
    else:
        return epoch

def epoch_to_jd(epoch,julian=True):
    """
    Convert a Julian or Besselian Epoch to a Julian Day.
    
    `jepoch` can be a string (if the string has a B or J at the beginning, the
    `julian` argument is ignored), or a float.
    
    If `julian` is True, a Julian Epoch will be used (the year is exactly 365.25
    days long). Otherwise, the epoch will be Besselian (assuming a tropical year
    of 365.242198781 days). 
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epj.html
    """
    if isinstance(epoch,basestring):
        if epoch[0]=='J':
            julian = True
            epoch = epoch[1:]
        elif epoch[0]=='B':
            julian = False
            epoch = epoch[1:]
        epoch = float(epoch)
    
    if julian:
        return (epoch - 2000)*365.25 + 2451545.0
    else:
        return (epoch - 1900)*365.242198781 + 2415020.31352
    
def greenwich_sidereal_time(jd,apparent=True):
    """
    Returns the sidereal time at the 0 degrees longitude on a given Julian Date
    `jd` (or dates if `jd` is an array). If `apparent` is True, the apparent
    sidereal time is returned, otherwise, it is the mean.
    """    
    #algorithm described on USNO web site http://aa.usno.navy.mil/faq/docs/GAST.php
    jd0 = np.round(jd-.5)+.5
    h = (jd - jd0) * 24.0
    d = jd - 2451545.0
    d0 = jd0 - 2451545.0
    t = d/36525
    
    #mean sidereal time @ greenwich
    gmst = 6.697374558 + 0.06570982441908*d0 + 0.000026*t**2 + 1.00273790935*h
           #- 1.72e-9*t**3 #left off as precision to t^3 is unneeded
   
    if apparent:
        eps =  np.radians(23.4393 - 0.0000004*d) #obliquity
        L = np.radians(280.47 + 0.98565*d) #mean longitude of the sun
        omega = np.radians(125.04 - 0.052954*d) #longitude of ascending node of moon
        dpsi = -0.000319*np.sin(omega) - 0.000024*np.sin(2*L) #nutation longitude
        return (gmst + dpsi*np.cos(eps))%24.0
    else:
        return gmst%24.0 
    
    

class Site(object):
    """
    This class represents a location on Earth from which skies are observable.
    
    lat/long are coords.AngularCoordinate objects, altitude in meters
    """
#    tznames = {'EST':-5,'CST':-6,'MST':-7,'PST':-8,
#               'EDT':-4,'CDT':-5,'MDT':-6,'PDT':-7}
    def __init__(self,lat,long,alt=0,tz=None,name=None):
        """
        generate a site by specifying latitude and longitude (as 
        coords.AngularCoordinate objects or initializers for one), optionally
        providing altitude (in meters), time zone (either as a timezone name
        provided by the system or as an offset from UTC), and/or a site name.
        """
        self.latitude = lat
        self.latitude.range = (-90,90)
        self.longitude = long
        self.longitude.range = (-180,180)
        self.altitude = alt
        if tz is None:
            self.tz = self._tzFromLong(self._long)
        elif isinstance(tz,basestring) and tzmod is not None:
            self.tz = tzmod.gettz(tz)
            if self.tz is None:
                raise ValueError('unrecognized time zone string '+tz)
        else:
            self.tz = tzoffset(str(tz),int(tz*60*60))
        if name is None:
            name = 'Default Site'
        self.name = name
        self._currjd = None
    
    @staticmethod
    def _tzFromLong(long):
        offset =  int(np.round(long.degrees/15))
        return tzoffset(str(offset),offset*60*60)
    
    def _getLatitude(self):
        return self._lat
    def _setLatitude(self,val):
        from .coords import AngularCoordinate
        from operator import isSequenceType
        if isinstance(val,AngularCoordinate):
            self._lat = val
        elif isSequenceType(val) and not isinstance(val,basestring):
            self._lat = AngularCoordinate(*val)
        else:
            self._lat = AngularCoordinate(val)
    latitude = property(_getLatitude,_setLatitude,doc='Latitude of the site in degrees')
    
    def _getLongitude(self):
        return self._long
    def _setLongitude(self,val):
        from .coords import AngularCoordinate
        from operator import isSequenceType
        if isinstance(val,AngularCoordinate):
            self._long = val
        elif isSequenceType(val) and not isinstance(val,basestring):
            self._long = AngularCoordinate(*val)
        else:
            self._long = AngularCoordinate(val)
    longitude = property(_getLongitude,_setLongitude,doc='Longitude of the site in degrees')
    
    def _getAltitude(self):
        return self._alt
    def _setAltitude(self,val):
        if val is None:
            self._alt = None
        else:
            self._alt = float(val)
    altitude = property(_getAltitude,_setAltitude,doc='Altitude of the site in meters')
    
    def _getCurrentobsjd(self):
        if self._currjd is None:
            from datetime import datetime
            return calendar_to_jd(datetime.utcnow(),tz=None)
        else:
            return self._currjd
    def _setCurrentobsjd(self,val):
        if val is None:
            self._currjd = None
        else:
            if np.isscalar(val):
                self._currjd = val
            else:
                self._currjd = calendar_to_jd(val)
    currentobsjd = property(_getCurrentobsjd,_setCurrentobsjd,doc="""
    Date and time to use for computing time-dependent values.  If set to None,
    the jd at the instant of calling will be used.  It can also be set as
    datetime objects or (yr,mon,day,hr,min,sec) tuples.
    """)
    
    def localSiderialTime(self,*args,**kwargs):
        """
        Compute the local siderial time given an input civil time. The various
        input forms are used to determine the interpretation of the civil time:
        
        * localSiderialTime()
            current local siderial time for this Site or uses the value of the
            :attr:`currentobsjd` property.
        * localSiderialTime(JD)
            input argument is julian date UT1
        * localSdierialTime(:class:`datetime.date`)
            compute the local siderial time for midnight on the given date
        * localSiderialTime(:class:`datetime.datetime`)
            the datetime object specifies the local time. If it has tzinfo, the
            object's time zone will be used, otherwise the :class:`Site's<Site>`
        * localSiderialTime(time,year,month,day)
            input arguments determine local time - time is in hours
        * localSiderialTime(year,month,day,hr,min,sec)
            local time - hours and minutes will be interpreted as integers
          
        *keywords*  
        
        * `apparent`
            if True (default) the returned time will be local apparent sidereal
            time (i.e. nutation terms included), otherwise it will be local mean
            sidereal time
        * `returntype`
            a string that determines the form of the returned LST as described
            below
        
        returns the local siderial time in a format that depends on the
        `returntype` keyword. It can be:
    
        * None/'hours' (default)
            LST in decimal hours
        * 'string'
            LST as a hh:mm:ss.s 
        * 'datetime'
            a :class:`datetime.time` object 
       
        """
        #guts of calculation adapted from xidl 
        import datetime
        
        rettype = kwargs.pop('returntype',None)
        apparent = kwargs.pop('apparent',True)
        if len(kwargs)>0:
            raise TypeError('got unexpected argument '+kwargs.keys()[0])
        if len(args)==0:
            jd = self.currentobsjd
        elif len(args)==1:
            if hasattr(args[0],'year'):
                if hasattr(args[0],'hour'):
                    if args[0].tzinfo is None:
                        dtobj = datetime.datetime(args[0].year,args[0].month,
                                        args[0].day,args[0].hour,args[0].minute,
                                        args[0].second,args[0].microsecond,self.tz)
                    else:
                        dtobj = args[0]
                else: #only date provided
                    dtobj = datetime.datetime(args[0].year,args[0].month,
                                              args[0].day,tzinfo=self.tz)
                jd = calendar_to_jd(dtobj,tz=None)
            else:
                jd = args[0]
        elif len(args) == 4:
            time,year,month,day = args
            hr = int(np.floor(time))
            min = int(np.floor(60*(time - hr)))
            sec = int(np.floor(60*(60*(time-hr) - min)))
            msec = int(np.floor(1e6*(60*(60*(time-hr) - min) - sec)))
            jd = calendar_to_jd(datetime.datetime(year,month,day,hr,min,sec,msec,self.tz),tz=None)
        elif len(args) == 6:
            year,month,day,hr,min,sec = args
            msec = int(1e6*(sec - np.floor(sec)))
            sec = int(np.floor(sec))
            jd = calendar_to_jd(datetime.datetime(year,month,day,hr,min,sec,msec,self.tz),tz=None)
        else:
            raise TypeError('invalid number of input arguments')
        
        lst = (greenwich_sidereal_time(jd,apparent) + self._long.d/15)%24.0 
        
#        #from idl astro ct2lst.pro         
#        jd2000 = 2451545.0
#        t0 = jd - jd2000
#        t = t0//36525 #TODO:check if true-div
        
#        #Compute GMST in seconds.  constants from Meeus 1ed, pg 84
#        c1,c2,c3,c4 = 280.46061837,360.98564736629,0.000387933,38710000.0
#        theta = c1 + (c2 * t0) + t**2*(c3 - t/ c4 )
        
#        #TODO: Add in mean->apparent corrections
        
#        #Compute LST in hours.
#        lst = np.array((theta + self._long.d)/15.0) % 24.0
        
        if lst.shape == tuple():
            lst = lst.ravel()[0]
        
        if rettype is None or rettype == 'hours':
            return lst
        elif rettype == 'string':
            hr = int(lst)
            min = 60*(lst - hr)
            sec = 60*(min - int(min))
            min = int(min)
            return '%i:%i:%f'%(hr,min,sec)
        elif rettype == 'datetime':
            hr = int(lst)
            min = 60*(lst - hr)
            sec = 60*(min - int(min))
            min = int(min)
            msec = int(1e6*(sec-int(sec)))
            sec = int(sec)
            return datetime.time(hr,min,sec,msec)
        else:
            raise ValueError('invalid returntype argument')
        
    def localTime(self,lsts,date=None,apparent=True,returntype=None,utc=False):
        """
        Computes the local civil time given a particular local siderial time. 
        
        `lsts` are the input local times, and may be either a scalar or array,
        and must be in decimal hours. Althernatively, it can be a
        :class:`datetime.datetime` object, in which case the date will be
        inferred from this object and the `date` argument will be ignored.
        
        `date` should be a :class:`datetime.date` object or a (year,month,day)
        tuple.  If None, the current date will be assumed as inferred from 
        :attr:`Site.currentobsjd`
        
        if `lsts` is a :class:`datetime.datetime` object, the object will be
        interpreted as the local siderial time with the corresponding date (any
        timezone info will be ignored), and the `date` argument will be ignored.
                        
        if `apparent` is True, the inputs are assumed to be in local apparent
        siderial time (i.e. nutation terms included), otherwise local mean
        siderial time.
                              
                              
        returns the local time in a format that depends on the `returntype`
        keyword. It can be:
        
        * None/'hours' (default)
            local time in decimal hours
        * 'string'
            local time as a hh:mm:ss.s 
        * 'datetime'
            a :class:`datetime.time` object with the appropriate tzinfo
            
        If `utc` is True, the time will be converted to UTC before being
        returned.
        
        """
        import datetime
        from operator import isSequenceType
        
        if isinstance(lsts,datetime.datetime):
            utcoffset = lsts.replace(tzinfo=self.tz).utcoffset()
            date = lsts.date()
            lsts = datetime.time()
            lsts = lsts.hour+lsts.minute/60+(lsts.second+lsts.microsecond*1e6)/3600
        else:
            jds,dt = self._processDate(date)
            utcoffset = dt.replace(tzinfo=self.tz).utcoffset()
            if isSequenceType(dt):
                raise ValueError('must provide only one date for localTime')
            date = dt.date()
                
        lsts = np.array(lsts,copy=False)
        scalarout = False
        if lsts.shape == tuple():
            scalarout = True
            lsts = lsts.ravel()
        
            
        lst0 = self.localSiderialTime(date)
        lthrs = (lsts - lst0)%24
        dayoffs = np.floor(lsts - lst0/24)
        
        lthrs /= 1.0027378507871321 #correct for siderial day != civil day
        
        if utc:
            lthrs = (lthrs - utcoffset.days*24 - utcoffset.seconds/3600)%24
        
        if returntype is None or returntype == 'hours':
            res = lthrs
        elif returntype == 'string':
            res = []
            for lthr in lthrs: 
                hr = int(lthr)
                min = 60*(lthr - hr)
                sec = 60*(min - int(min))
                min = int(min)
                res.append('%i:%i:%f'%(hr,min,sec))
        elif returntype == 'datetime':
            res = []
            for i,dayoff in zip(lthrs,dayoffs): 
                hr = int(lthr)
                min = 60*(lthr - hr)
                sec = 60*(min - int(min))
                min = int(min)
                msec = int(1e6*(sec-int(sec)))
                sec = int(sec)
                
                res.append(datetime.datetime(date.year,date.month,date.day+dayoff,hr,min,sec,msec,self.tz))
        else:
            raise ValueError('invalid returntype argument')
        
        if scalarout:
            return res[0]
        else:
            return res
       
    
    def equatorialToHorizontal(self,eqpos,lsts,epoch=None):
        """
        Generates a list of horizontal positions (or just one) from a provided
        equatorial position and local siderial time (or sequence of LSTs) in 
        decimal hours and a latitude in degrees.  If  `epoch` is not None, it 
        will be used to set the epoch in the equatorial system.
        
        Note that this generally should *not* be used to compute the observed 
        horizontal coordinates directly - :meth:`apparentCoordinates` includes all 
        the corrections and formatting.  This method is purely for the 
        coordinate conversion.
        """
        from .coords import EquatorialCoordinates,HorizontalCoordinates
        
        latitude = self.latitude.d
        
        if epoch is not None:
            #make copy so as not to change the epoch of the input
            eqpos = EquatorialCoordinates(eqpos)
            eqpos.epoch = epoch
            
        lsts = np.array(lsts,copy=False)
        singleout = lsts.shape == tuple()
        
        HA = lsts.ravel() - eqpos.ra.hours 
        sHA = np.sin(pi*HA/12)
        cHA = np.cos(pi*HA/12)
        
        sdec = np.sin(eqpos.dec.radians)
        cdec = np.cos(eqpos.dec.radians)
        slat = np.sin(np.radians(latitude))
        clat = np.cos(np.radians(latitude))
        
        alts = np.arcsin(slat*sdec+clat*cdec*cHA)
        calts = np.cos(alts)
        azs = np.arctan2(-cdec*sHA,clat*sdec-slat*cdec*cHA)%(2*pi)
        
        if eqpos.decerr is not None or eqpos.raerr is not None:
            decerr = eqpos.decerr.radians if eqpos.decerr is not None else 0
            raerr = eqpos.raerr.radians if eqpos.raerr is not None else 0
            
            dcosalt = np.cos(alts)
            daltdH = -clat*cdec*sHA/dcosalt
            daltddec = (slat*cdec-clat*sdec*cHA)/dcosalt
            
            dalts = ((daltdH*raerr)**2 + (daltddec*decerr)**2)**0.5
            
            #error propogation computed with sympy following standard rules
            dtanaz = 1+np.tan(azs)**2
            dazdH = (cHA*cdec/(cHA*cdec*slat-clat*sdec) \
                  + cdec*cdec*sHA*sHA*slat*(cHA*cdec*slat-clat*sdec)**-2) \
                      /dtanaz
            dazddec = ((sHA*sdec)/(clat*sdec - cHA*cdec*slat) \
                     + ((cdec*clat+cHA*sdec*slat)*cdec*sHA)*(cHA*cdec*slat-clat*sdec)**-2) \
                      /dtanaz
            
            dazs = ((dazdH*raerr)**2 + (dazddec*decerr)**2)**0.5
        else:
            dazs = dalts = None
                
        if singleout:
            if dazs is None:
                return HorizontalCoordinates(*np.degrees((alts[0],azs[0])))
            else:
                return HorizontalCoordinates(*np.degrees((alts[0],azs[0],dazs[0],dalts[0])))
        else:
            if dazs is None:
                return [HorizontalCoordinates(alt,az) for alt,az in np.degrees((alts,azs)).T]
            else:
                return [HorizontalCoordinates(alt,az,daz,dalt) for alt,az,daz,dalt in np.degrees((alts,azs,dazs,dalts)).T]
            
    def riseSetTransit(self,eqpos,date=None,alt=-.5667,timeobj=False,utc=False):
        """
        Computes the rise, set, and transit times of a provided equatorial 
        position in local time.  
        
        `alt` determines the altitude to be considered as risen or set in 
        degrees. Default is for approximate rise/set including refraction.
        
        `date` determines the date at which to do the computation. See
        :func:`calendar_to_jd` for acceptable formats
        
        *returns*
        (rise,set,transit) as :class:datetime.time objects if `timeobj` is True
        or if it is False, they are in decimal hours.  If the object is 
        circumpolar, rise and set are both None.  If it is never visible,
        rise,set, and transit are all None
        """
        import datetime
        from dateutil import tz
        from math import sin,cos,radians,acos
        
        alt = radians(alt)
        
        transit = self.localTime(eqpos.ra.hours,date,utc=utc)
        
        #TODO:iterative algorithm for rise/set
        lat = self.latitude.radians
        dec = eqpos.dec.radians
        #local hour angle for alt
        coslha = (sin(alt) - sin(lat)*sin(dec))/(cos(lat)*cos(dec))
        if coslha<-1:
            #circumpolar
            rise=set=None
        elif coslha>1:
            #never visible
            rise=set=transit=None
        else:
            lha = acos(coslha)*12/pi
            lha /= 1.0027378507871321 #correction for siderial day != calendar day
            rise = (transit - lha)%24 
            set = (transit + lha)%24 #TODO:iterative algorithm
        
        if timeobj:
            def hrtotime(t):
                if t is None:
                    return None
                else:
                    hr = np.floor(t)
                    t = (t-hr)*60
                    min = np.floor(t)
                    t = (t-min)*60
                    sec = np.floor(t)
                    t = (t-sec)*1e6
                    msec = np.floor(t)
                    return datetime.time(int(hr),int(min),int(sec),int(msec),
                                    tzinfo=tz.tzutc() if utc else tz.tzlocal())
            return hrtotime(rise),hrtotime(set),hrtotime(transit)
        else:
            return rise,set,transit
        
        
    def apparentCoordinates(self,coords,datetime=None,precess=True,refraction=True):
        """
        computes the positions in horizontal coordinates of an object with the 
        provided fixed coordinates at the requested time(s).
        
        `datetime` can be a sequence of datetime objects or similar representations,
        or it can be a sequnce of JD's.  If None, it uses the current time or 
        the value of the :attr:`currentobsjd` property.
        
        If `precess` is True, the position will be precessed to the epoch of 
        the observation (almost always the right thing to do)
        
        If `refraction` is True, an added correction to the altitude due to
        atmospheric refraction at STP (formula from Meeus ch 16) is included. If
        `refraction` is a (non-0) float, it will be taken as the temperature at
        which to perform the refraction calculation. If it evaluates to False,
        no refraction correction is performed
        
        
        *returns*
        a sequence of :class:`astropysics.coords.HorizontalCoordinates` objects
        """
        from .coords import HorizontalCoordinates
        from operator import isSequenceType
        
        if datetime is None:
            return self.equatorialToHorizontal(coords,self.localSiderialTime())
        elif hasattr(datetime,'year') or (isSequenceType(datetime) and hasattr(datetime[0],'year')):
            jd = calendar_to_jd(datetime,self.tz).ravel()
        else:
            jd = np.array(datetime,copy=False)
            if len(jd.shape)>1:
                jd = np.array([calendar_to_jd(v,self.tz) for v in jd])
        lsts = self.localSiderialTime(jd)
        
        if precess:
            res = self.equatorialToHorizontal(coords,lsts,epoch=jd_to_epoch(jd[0]))
        else:
            res = self.equatorialToHorizontal(coords,lsts)
            
        if refraction:
            from math import tan,radians
            if not hasattr(self,'_cached_pressure') or self._cached_pressure is None:
                from math import exp
                from.constants import g0,Rb
                t0 = 273 #K
                M = 28.9644 #g/mol
                self._cached_pressure = P = exp(g0*M*self.altitude/(Rb*t0))
            else:
                P = self._cached_pressure
                
            if refraction is True:
                T = 273
            else:
                T  = float(refraction)
            
            h = res.alt.radians
            R = 1.02/tan(h+(10.3/(h+5.11))) #additive correction in arcmin
            #for inverse problem of apparent h->true/airless h, use:
            #R = 1/tan(h0+(7.31/(h0+4.4)))
            
            res.alt._decval += radians(R/60.)
        
        return res
        
    def _processDate(self,date):
        """
        utitily function to convert a date in a variety of formats to a 
        tuple jds,datetimes
        """
        if date is None:
            jd = self.currentobsjd
        elif np.isscalar(date):
            jd = date
        else:
            jd = calendar_to_jd(date)
            
        return jd,jd_to_calendar(jd)
        
    def nightTable(self,coord,date=None,strtablename=None,localtime=True,
                            hrrange=(18,6,25)):
        """
        tabulates altitude, azimuth, and airmass values for the provided fixed
        position on a particular date, specified either as a datetime.date 
        object, a (yr,mon,day) tuple, or a julain date (will be rounded to 
        nearest)
        
        if `date` is None, the current date is used
        
        If `localtime` is True, the hours output (and input) will be in local
        time for this Site. Otherwise, it is UTC.
        
        `hrrange` determines the size of the table - it should be a 3-tuple
        (starthr,endhr,n) where starthr is on the day specified in date and 
        endhr is date + 1 day
        
        if `strtablename` is True, a string is returned with a printable table 
        of observing data.  If `strtable` is a string, it will be used as the 
        title for the table.  Otherwise, a record array is returned with the 
        hour(UTC), alt, az, and airmass
        """        
        from .astropysics import coords
        import datetime
        
        #for objects that can get a position with no argument
        if hasattr(coord,'equatorialCoordinates'):
            coord = coord.equatorialCoordinates()
        
        if date is None:
            jd = self.currentobsjd
        elif np.isscalar(date):
            jd = date
        else:
            jd = calendar_to_jd(date)
            
        date = jd_to_calendar(jd).date()
        jd0 = np.floor(jd)
        
        if localtime:
            dt = datetime.datetime.combine(date,datetime.time(12,tzinfo=self.tz))
            offs = dt.utcoffset()
            utcoffset = offs.days*24+offs.seconds/3600
        else:
            utcoffset = 0
        
        starthr,endhr,n = hrrange
        starthr -= utcoffset #local->UTC
        endhr -= utcoffset #local->UTC
        
        startjd = jd0 + (starthr - 12)/24
        endjd = jd0 + (endhr + 12)/24
        jds = np.linspace(startjd,endjd,n)
        
        timehr = (jds-np.round(np.mean(jds))+.5)*24+utcoffset #UTC hr
        
        
        alt,az = coords.objects_to_coordinate_arrays(self.apparentCoordinates(coord,jds,precess=False))
        
        z = 90 - alt
        airmass = 1/np.cos(np.radians(z))
        
        ra = np.rec.fromarrays((timehr,alt,az,airmass),names = 'hour,alt,az,airmass')

        
        if strtablename is not None:
            rise,set,transit = self.riseSetTransit(coord,date=date,timeobj=True)
            
            lines = []
            
            if isinstance(strtablename,basestring):
                lines.append('Object:'+str(strtablename))
            lines.append('{0} {1}'.format(coord.ra.getHmsStr(),coord.dec.getDmsStr()))
            lines.append(str(date))
            if transit is None:
                lines.append('Not visible from this site!')
            elif rise is None:
                lines.append('Circumpolar')
                lines.append('Transit : {0}'.format(transit))
            else:
                lines.append('Rise : {0}'.format(rise))
                lines.append('Set : {0}'.format(set))
                lines.append('Transit : {0}'.format(transit))
            lines.append('')
            lines.append('time\talt\taz\tairmass')
            lines.append('-'*(len(lines[-1])+lines[-1].count('\t')*4))
            
            currdate = date+datetime.timedelta(int(np.floor(starthr/24)))
            lines.append('\t'+str(currdate))
            
            oldhr24 = int(np.floor(ra[0].hour))
            for r in ra:
                hr = int(np.floor(r.hour))
                if hr%24 < oldhr24:
                    currdate = currdate+datetime.timedelta(1)
                    lines.append('\t'+str(currdate))
                min = int(np.round((r.hour-hr)*60))
                if min==60:
                    hr+=1
                    min = 0
                alt = r.alt
                az = r.az
                am = r.airmass if r.airmass>0 else '...'
                line = '{0:2}:{1:02}\t{2:.3}\t{3:.4}\t{4:.3}'.format(hr%24,min,alt,az,am)
                lines.append(line)
                oldhr24 = hr%24
            return '\n'.join(lines)
        else:
            rise,set,transit = self.riseSetTransit(coord,date=date,timeobj=False)
            ra.sitedate = date
            ra.utcoffset = utcoffset
            ra.rise = rise
            ra.set = set
            ra.transit = transit
            return ra
        
    def nightPlot(self,coords,date=None,plottype='altam',onlynight=False,
                       moon=True,sun=True,clf=True,utc=False,colors=None,
                       plotkwargs=None):
                           
        """
        Generates plots of important observability quantities for the provided
        coordinate objects for a single night.
        
        `coords` should be a :class:`astropysics.coords.LatLongCoordinates`, a
        sequence of such objects, or a dictionary mapping the name of an object
        to the object itself. These names will be used to label the plot.
        
        `plottype` can be one of:
        
        * 'altam' : a plot of time vs. altitude with a secondary axis for sec(z)
        * 'am' : a plot of time vs. airmass/sec(z)
        * 'altaz': a plot of azimuth vs. altitude
        * 'sky': polar projection of the path on the sky
        
        If `onlynight` is True, the plot will be shrunk to only show the times 
        between sunset and sunrise.  Otherwise, an entire day/night will be 
        plotted.
        
        If `moon` is True, the path of the Moon will be plotted on relevant
        plots. If it is a dictionary, they will be passed as kwargs to
        :func:`matplotlib.pyplot.plot` for the moon.
        
        If `sun` is True, shaded regions for 18 degree,12 degree, and
        sunrise/set will be added where relevant.
        
        If `clf` is True, the figure is cleared before the observing plot is
        made.
        
        If `utc` is True, the times will be in UTC instead of local time.
        
        `colors` should be a sequence of :mod:`matplotlib` color specifiers, or
        None to use the default color cycle.
        
        `plotkwargs` will be provided as a keyword dictionary to
        :func:`matplotlib.pyplot.plot`, unless it is None
        """
        import matplotlib.pyplot as plt
        from operator import isMappingType,isSequenceType
        from .coords import LatLongCoordinates,Sun,Moon
        from .plotting import add_mapped_axis
        
        if isMappingType(coords):
            names = coords.keys()
            coords = coords.values()
            nonames = False
        elif not isSequenceType(coords):
            coords = [coords]
            names = ['' for c in coords]
            nonames = True
        else:
            names = ['' for c in coords]
            nonames = True
            
        if colors:
            plt.gca().set_color_cycle(colors)
            
        if isMappingType(plotkwargs):
            plotkwargs = [plotkwargs for c in coords]    
        elif plotkwargs is None:
            plotkwargs = [None for c in coords]
        
        inter = plt.isinteractive()
        try:
            plt.ioff()
            if clf:
                plt.clf()
            oldright = None
            if plottype == 'altam' or plottype == 'alt' or plottype == 'am':   
                if plottype == 'altam':
                    oldright = plt.gcf().subplotpars.right
                    plt.subplots_adjust(right=0.86)
                    
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.nightTable(c,date,hrrange=(12,12,100),localtime=True)
                    if kw is None:
                        kw = {}
                    kw.setdefault('label',n)
                    kw.setdefault('zorder',3)
                    kw.setdefault('lw',2)
                    if plottype == 'am':
                        ammask = ra.airmass>0
                        x = ra.hour[ammask]
                        y = ra.airmass[ammask]
                    else:
                        x = ra.hour
                        y = ra.alt
                    plt.plot(x,y,**kw)
                
                
                plt.title(str(ra.sitedate))
                 
                
                if 'alt' in plottype:
                    if plt.ylim()[0] < 0:
                        lowery = 0
                    else:
                        lowery = plt.ylim()[0]
                    uppery = plt.ylim()[1] if plt.ylim()[1]<90 else 90
                    plt.ylim(lowery,uppery)
                    plt.ylabel(r'${\rm altitude} [{\rm degrees}]$')
                elif 'am' == plottype:
                    if plt.ylim()[0] < 1:
                        lowery = 1
                    else:
                        lowery = plt.ylim()[0]
                    uppery = plt.ylim()[1]
                    plt.ylim(lowery,uppery)
                    if lowery==1:
                        yticks = list(plt.yticks()[0])
                        yticks.insert(0,1)
                        plt.yticks(yticks)
                        plt.ylim(uppery,lowery)
                    plt.ylabel(r'${\rm airmass} / \sec(z)$')
                    
                if not nonames:
                    plt.legend(loc=0)
                    
                if sun:
                    seqp = Sun(ra.sitedate).equatorialCoordinates()
                    rise,set,t = self.riseSetTransit(seqp,date,0)
                    rise12,set12,t12 = self.riseSetTransit(seqp,date,-12)
                    rise18,set18,t18 = self.riseSetTransit(seqp,date,-18)
                    #need to correct them back to the previous day
                    set -= 24
                    set12 -= 24
                    set18 -= 24
                    
                    xl,xu = plt.xlim()
                    yl,yu = plt.ylim()
                    
                    grey1 = (0.5,0.5,0.5)  
                    grey2 = (0.65,0.65,0.65)
                    grey3 = (0.8,0.8,0.8)
                    plt.fill_between((set,set12),lowery,uppery,lw=0,color=grey3,zorder=1)
                    plt.fill_between((set12,set18),lowery,uppery,lw=0,color=grey2,zorder=1)
                    plt.fill_between((set18,rise18),lowery,uppery,lw=0,color=grey1,zorder=1)
                    plt.fill_between((rise18,rise12),lowery,uppery,lw=0,color=grey2,zorder=1)
                    plt.fill_between((rise12,rise),lowery,uppery,lw=0,color=grey3,zorder=1)
                    
                    
                    
                    plt.xlim(xl,xu)
                    plt.ylim(yl,yu)
                    
                if moon:
                    xls = plt.xlim()
                    yls = plt.ylim()
                    
                    m = Moon(ra.sitedate)
                    ram = self.nightTable(m,date,hrrange=(12,12,100),localtime=True)
                    if not isMappingType(moon):
                        moon = {}
                    moon.setdefault('label','Moon')
                    moon.setdefault('zorder',2)
                    moon.setdefault('ls','--')
                    moon.setdefault('lw',1)
                    moon.setdefault('c','k')
                    if plottype == 'am':
                        ammask = ram.airmass>0
                        plt.plot(ram.hour[ammask],ram.airmass[ammask],**moon)
                    else:
                        plt.plot(ram.hour,ram.alt,**moon)
                    
                    plt.xlim(*xls)
                    plt.ylim(*yls)
                    
                if plottype == 'altam':
                    firstax = plt.gca()
                    yticks = plt.yticks()[0]
                    newticks = list(yticks[1:])
                    newticks.insert(0,(yticks[0]+newticks[0])/2)
                    newticks.insert(0,(yticks[0]+newticks[0])/2)
                    plt.twinx()
                    plt.ylim(lowery,uppery)
                    plt.yticks(newticks,['%2.2f'%(1/np.cos(np.radians(90-yt))) for yt in newticks])
                    plt.ylabel(r'$\sec(z)$')
                    plt.axes(firstax)
                    
                xtcks = np.round(np.arange(13)*(x[-1]-x[0])/12+x[0]).astype(int)
                if utc:
                    plt.xticks(xtcks,['%i'%np.round(xt%24) for xt in xtcks])
                    plt.xlabel(r'$\rm UTC (hours)$')
                else:
                    print 'here'
                    plt.xticks(xtcks,[xt%24 for xt in xtcks])
                    plt.xlabel(r'$\rm Local Time (hours)$')
                    
                if onlynight:
                    plt.xlim(xtcks[0]/2,xtcks[-1]/2)
                else:
                    plt.xlim(xtcks[0],xtcks[-1])  
                    
                    
            elif plottype == 'altaz':
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.nightTable(c,date,hrrange=(0,0,100))
                    
                    if kw is None:
                        plt.plot(ra.az,ra.alt,label=n)
                    else:
                        kw['label'] = n
                        plt.plot(ra.az,ra.alt,**kw)
                        
                plt.xlim(0,360)
                plt.xticks(np.arange(9)*360/8)
                plt.ylim(0,90)
                
                plt.title(str(ra.sitedate))
                plt.xlabel(r'${\rm azimuth} [{\rm degrees}]$')
                plt.ylabel(r'${\rm altitude} [{\rm degrees}]$')
                
                if not nonames:
                    plt.legend(loc=0)
                
            elif plottype == 'sky':
                for n,c,kw in zip(names,coords,plotkwargs):
                    
                    ra = self.nightTable(c,date,hrrange=(0,0,100))
                    if kw is None:
                        plt.polar(np.radians(ra.az),90-ra.alt,label=n)
                    else:
                        kw['label'] = n
                        plt.polar(np.radians(ra.az),90-ra.alt,**kw)
                
                xticks = [0,45,90,135,180,225,270,315]      
                xtlabs = ['N',r'$45^\circ$','E',r'$135^\circ$','S',r'$225^\circ$','W',r'$315^\circ$']  
                plt.xticks(np.radians(xticks),xtlabs)
                plt.ylim(0,90)
                yticks = [15,30,45,60,75]
                plt.yticks(yticks,[r'${0}^\circ$'.format(int(90-yt)) for yt in yticks])
                
                plt.title(str(ra.sitedate))
                        
                if not nonames:
                    plt.legend(loc=0)
            else:
                raise ValueError('unrecognized plottype {0}'.format(plottype))
            if inter:
                plt.show()
                plt.draw()
            if oldright is not None:
                plt.gcf().subplotpars.right = oldright
        finally:
            plt.interactive(inter)
            
    def yearPlot(self,coords,startdate=None,n=13,months=12,sun=True,moon=True,
                      utc=False,clf=True,colors=None):
        """
        Plots the location transit and rise/set times of object(s) over ~year
        timescales.
        
        `coords` should be a :class:`astropysics.coords.LatLongCoordinates`, a
        sequence of such objects, or a dictionary mapping the name of an object
        to the object itself. These names will be used to label the plot.
        
        `startdate` is a :class:`datetime.date` object or a date tuple
        (year,month,day) that is used as the start of the plot.
        
        `n` specifies the number of points to include of the objects
        
        `months` is the number of months to draw the plot for.        
        
        If `sun` is True, shaded regions for 18 degree,12 degree, and
        sunrise/set will be included.
        
        If `moon` is True, the path of the Moon will be plotted.
        
        If `utc` is True, the times will be in UTC instead of local time.
        
        If `clf` is True, the figure is cleared before the observing plot is
        made.
        
        `colors` should be a sequence of :mod:`matplotlib` color specifiers, or
        None to use the default color cycle.
        """
                      
        import matplotlib.pyplot as plt
        import datetime
        from matplotlib.dates import MonthLocator,WeekdayLocator,DateFormatter, \
                                     MONDAY,DayLocator,YearLocator
        from operator import isMappingType,isSequenceType
        from .coords import Sun,Moon
        
        if isMappingType(coords):
            names = coords.keys()
            coords = coords.values()
        elif not isSequenceType(coords):
            coords = [coords]
            names = ['' for c in coords]
        else:
            names = ['' for c in coords]
            
        jdstart,startdt = self._processDate(startdate)
        jds = jdstart+np.linspace(0,365.25*months/12,n)
        jd1 = jds - calendar_to_jd((1,1,1))
        
        if utc:
            utco = startdt.replace(tzinfo=self.tz).utcoffset()
            center = -utco.days*24 - utco.seconds/3600
        else:
            center = 0
        cp12 = center + 12
        
        inter = plt.isinteractive()
        try:
            plt.ioff()
            
            if clf:
                plt.clf()
                
            if colors:
                plt.gca().set_color_cycle(colors)
            for c,nm in zip(coords,names):
                rst = [self.riseSetTransit(c,jd_to_calendar(jd),0,utc=utc) for jd in jds]
                rise,set,transit = np.array(rst).T
                transit[transit>cp12] -= 24 #do this to get the plot to cross over night time
                
                c = plt.gca()._get_lines.color_cycle.next()
                
                tline = plt.plot_date(jd1,transit,label=nm,color=c)[0]
                rise[rise>transit]-=24
                set[set<transit]+=24
                plt.errorbar(jd1,transit,(set-transit,transit-rise),ecolor=c,fmt=None)
                
                
            if sun:
                jdsun = jdstart+np.linspace(0,365.25*months/12,365)
                jd1sun = jdsun - calendar_to_jd((1,1,1))
                sun = Sun()
                    
                rst = []
                rst12 = []
                rst18 = []
                for jd in jdsun:
                    sun.jd = jd
                    rst.append(self.riseSetTransit(sun.equatorialCoordinates(),jd_to_calendar(jd),0,utc=utc))
                    rst12.append(self.riseSetTransit(sun.equatorialCoordinates(),jd_to_calendar(jd),-12,utc=utc))
                    rst18.append(self.riseSetTransit(sun.equatorialCoordinates(),jd_to_calendar(jd),-18,utc=utc))
                    
                rise,set,t = np.array(rst).T
                rise12,set12,t12 = np.array(rst12).T
                rise18,set18,t18 = np.array(rst18).T
                for a in (rise,set,rise12,set12,rise18,set18):
                    a[a>cp12] -= 24
                    
                grey1 = (0.5,0.5,0.5)  
                grey2 = (0.65,0.65,0.65)
                grey3 = (0.8,0.8,0.8)
                plt.fill_between(jd1sun,rise18,set18,color=grey1)
                plt.fill_between(jd1sun,rise18,rise12,color=grey2)
                plt.fill_between(jd1sun,rise12,rise,color=grey3)
                plt.fill_between(jd1sun,set12,set18,color=grey2)
                plt.fill_between(jd1sun,set,set12,color=grey3)
            
            if moon:
                jdmoon = jdstart+np.linspace(0,365.25*months/12,30*months)
                jd1moon = jdmoon - calendar_to_jd((1,1,1))
                moon = Moon()
                
                rst = []
                for jd in jdmoon:
                    moon.jd = jd
                    rst.append(self.riseSetTransit(moon.equatorialCoordinates(),jd_to_calendar(jd),0,utc=utc))
                rise,set,t = np.array(rst).T
                
                t[t>cp12] -= 24
                transitionmask = np.roll(t,1)>t
                transitionmask[0] = transitionmask[-1] = False
                
                lasti = 0
                for i in np.where(transitionmask)[0]:
                    lastplot = plt.plot(jd1moon[lasti:i],t[lasti:i],'--k')[0]
                    lasti = i
                if lasti == len(t):
                    lastplot._label = 'Moon'
                else:
                    plt.plot(jd1moon[lasti:],t[lasti:],'--k',label='Moon')
            if utc:
                plt.ylabel('UTC (hours)')
            else:    
                plt.ylabel('Local Time (hours)')
            plt.xlabel('Month (small ticks on Mondays)')
            
                
            
            if months <=3:
                plt.gca().xaxis.set_major_locator(MonthLocator())
                plt.gca().xaxis.set_major_formatter(DateFormatter("%b"))
                plt.gca().xaxis.set_minor_locator(DayLocator((6,10,15,20,25)))
                plt.gca().xaxis.set_minor_formatter(DateFormatter('$%d$'))
            elif 36 >= months > 12:
                plt.gca().xaxis.set_major_locator(MonthLocator(interval=3))
                plt.gca().xaxis.set_major_formatter(DateFormatter("%b '%y"))
                plt.gcf().autofmt_xdate()
            elif months > 36:
                plt.gca().xaxis.set_major_locator(YearLocator())
                plt.gca().xaxis.set_major_formatter(DateFormatter("%Y"))
                plt.gcf().autofmt_xdate()
            else:
                plt.gca().xaxis.set_major_locator(MonthLocator())
                plt.gca().xaxis.set_major_formatter(DateFormatter("%b '%y"))
                plt.gca().xaxis.set_minor_locator(WeekdayLocator(MONDAY))
                plt.gcf().autofmt_xdate()
            
            plt.xlim(jd1[0],jd1[-1])
            plt.ylim(center-12,center+12)
            ytcks = np.round(np.arange(13)*2-12+center).astype(int)
            plt.yticks(ytcks,ytcks%24)
            
            if any([nm!='' for nm in names]):
                plt.legend(loc=0)

            if inter:
                plt.show()
                plt.draw()
        finally:
            plt.interactive(inter)
        

def __loadobsdb(sitereg):
    from .io import _get_package_data
    obsdb = _get_package_data('obsdb.dat')
    from .coords import AngularCoordinate
    
    obs = None
    for l in obsdb.split('\n'):
        ls = l.strip()
        if len(ls)==0 or ls[0]=='#':
            continue
        k,v = [ss.strip() for ss in l.split('=')]
        if k == 'observatory':
            if obs is not None:
                o = Site(lat,long,alt,tz,name)
                sitereg[obs] = o
            obs = v.replace('"','')
            name = long = lat = alt = tz = None
        elif k == 'name':
            name = v.replace('"','')
        elif k == 'longitude':
            vs = v.split(':')
            dec = float(vs[0])
            if len(vs)>1:
                dec += float(vs[1])/60
            #longitudes here are deg W instead of (-180, 180)
            dec*=-1
            if dec <=-180:
                dec += 360
            long = AngularCoordinate(dec)
        elif k == 'latitude':
            vs = v.split(':')
            dec = float(vs[0])
            if len(vs)>1:
                dec += float(vs[1])/60
            lat = AngularCoordinate(dec)
        elif k == 'altitude':
            alt = float(v)
        elif k == 'timezone':
            #time zones are also flipped
            tz = -1*float(v)
    if obs is not None:
        o = Site(lat,long,alt,tz,name)
        sitereg[obs] = o


sites = DataObjectRegistry('sites',Site)
try:
    sites['uciobs'] = Site(33.63614044191056,-117.83079922199249,80,'US/Pacific','UC Irvine Observatory')
except ValueError: #in case US/Pacific is not present for some reason
    sites['uciobs'] = Site(33.63614044191056,-117.83079922199249,80,-8,'UC Irvine Observatory')
sites['greenwich'] = Site('51d28m38s',0,7,0,'Royal Observatory,Greenwich')
__loadobsdb(sites)
#<-----------------Attenuation/Reddening and dust-related---------------------->

class Extinction(PipelineElement):
    """
    This is the base class for extinction-law objects. Extinction laws can be
    passed in as a function to the initializer, or subclasses should override 
    the function f with the preferred law as f(self,lambda), and their 
    initializers should set the zero point
    
    Note that functions are interpreted as magnitude extinction laws ... if 
    optical depth is desired, f should return (-2.5/log(10))*tau(lambda)
    
    A0 is the normalization factor that gets multiplied into the reddening law.
    """
    def  __init__(self,f=None,A0=1):
        if f is not None:
            if not callable(f):
                self.f = f
            else:
                raise ValueError('function must be a callable')
        self.A0 = A0
        
        self._plbuffer = None
    
    def f(self,lamb):
        raise NotImplementedError('must specify an extinction function ')
    
    def __call__(self,*args,**kwargs):
        return self.A0*self.f(*args,**kwargs)
    
    def correctPhotometry(self,mags,band):
        """
        Uses the extinction law to correct a magnitude (or array of magnitudes)
        
        bands is either a string specifying the band, or a wavelength to use
        """
        from .phot import bandwl
        #TODO:better band support
        if isinstance(band,basestring):
            wl = bandwl[band]
        else:
            wl = band
        
        return mags-self(wl)
        
    def Alambda(self,band):
        """
        determines the extinction for this extinction law in a given band or bands
        
        band can be a wavelength or a string specifying a band 
        """
        from operator import isSequenceType
        if isSequenceType(band) and not isinstance(band,basestring):
            return -1*np.array([self.correctPhotometry(0,b) for b in band])
        else:
            return -1*self.correctPhotometry(0,band)
        
    def correctColor(self,colors,bands):
        """
        Uses the supplied extinction law to correct a color (or array of colors)
        where the color is in the specified bands
        
        bands is a length-2 sequence with either a band name or a wavelength for
        the band, or of the form 'bandname1-bandname2' or 'E(band1-band2)'
        """
        
        if isinstance(bands,basestring):
            b1,b2=bands.replace('E(','').replace(')','').split(',')
        else:
            b1,b2=bands
            
        from .phot import bandwl
        #TODO:better band support
        wl1 = bandwl[b1] if isinstance(b1,basestring) else b1
        wl2 = bandwl[b2] if isinstance(b1,basestring) else b2
    
        return colors-self(wl1)+self(wl2)
    
    def correctSpectrum(self,spec,newspec=True):
        """
        Uses the supplied extinction law to correct a spectrum for extinction.
        
        if newspec is True, a copy of the supplied spectrum will have the 
        extinction correction applied
        
        returns the corrected spectrum
        """
    
        if newspec:
            spec = spec.copy()
            
        oldunit = spec.unit
        spec.unit = 'wavelength-angstrom'
        corr = 10**(self(spec.x)/2.5)
        spec.flux *= corr
        spec.err *= corr
        
        spec.unit = oldunit
        return spec
    
    __builtinlines ={
            'Ha':6562.82,
            'Hb':4861.33,
            'Hg':4340.46,
            'Hd':4101.74,
            'He':3970.07
            }
    def correctFlux(self,flux,lamb):
        if isinstance(lamb,basestring) and lamb in Extinction.__builtinlines:
            lamb = Extinction.__builtinlines[lamb]
        return flux*10**(self(lamb)/2.5)
    
    
    __balmerratios={
            'Hab':(2.86,6562.82,4861.33),
            'Hag':(6.16,6562.82,4340.46),
            'Had':(11.21,6562.82,4101.74),
            'Hae':(18.16,6562.82,3970.07),
            'Hbg':(2.15,4861.33,4340.46),
            'Hbd':(3.91,4861.33,4101.74),
            'Hbe':(6.33,4861.33,3970.07),
            'Hgd':(1.82,4340.46,4101.74),
            'Hge':(2.95,4340.46,3970.07),
            'Hde':(1.62,4101.74,3970.07)
            }
    def computeA0FromFluxRatio(self,measured,expected,lambda1=None,lambda2=None,filterfunc=None):
        """
        This derives the normalization of the Extinction function from provided
        ratios for theoretically expected fluxes.  If multiple measurements are
        provided, the mean will be used
        
        measured is the measured line ratio (possibly an array), while expected
        is either the expected line ratios, or a string specifying the
        appropriate balmer flux ratio as "Hab","Hde",etc. (for Halpha/Hbeta or
        Hdelta/Hepsilon) to assume case B recombination fluxes (see e.g. 
        Osterbrock 2007).
        
        lambda1 and lambda2 are the wavelengths for the ratios F1/F2, or None if
        a string is provided 
    
        filterfunc is a function to be applied as the last step - it can 
        either be used to contract an array to (e.g. np.mean), or filter out
        invalid values (e.g. lambda x:x[np.isfinite(x)]).  Default
        does nothing
        
        returns A0,standard deviation of measurements
        """
        if isinstance(expected,basestring):
            R = measured
            balmertuple = Extinction.__balmerratios[expected]
            R0=balmertuple[0]
            lambda1,lambda2=balmertuple[1],balmertuple[2]
        elif isinstance(expected[0],basestring):
            if not np.all([isinstance(e,basestring) for e in expected]):
                raise ValueError('expected must be only transitions or only numerical')
            R0,lambda1,lambda2=[],[],[]
            for e in expected:
                balmertuple = Extinction.__balmerratios[e]
                R0.append(balmertuple[0])
                lambda1.append(balmertuple[1])
                lambda2.append(balmertuple[2])
            R0,lambda1,lambda2=np.array(R0),np.array(lambda1),np.array(lambda2)    
        else:
            if lambda1 and lambda2:
                R,R0 = np.array(measured,copy=False),np.array(expected,copy=False)
                lambda1,lambda2 = np.array(lambda1,copy=False),np.array(lambda2,copy=False)
            else:
                raise ValueError('need to provide wavelengths if not specifying transitions')
        
        A0 = -2.5*np.log10(R/R0)/(self.f(lambda1)-self.f(lambda2))
        
        if filterfunc is None:
            filterfunc = lambda x:x
        self.A0 = A0 = filterfunc(A0)
        return A0,np.std(A0)
    
    #PipelineElement methods
#    def _plFeed(self,data,src):
#        from .utils import  PipelineError
#        from .spec import Spectrum
#        if self._plbuffer is None:
#            self._plbuffer = {'in':[],'out':[]} 
    
#        if isinstance(data,Spectrum):
#            self._plbuffer['in'].append(('spec',data))
#        else:
#            raise PipelineError('unrecognized Extinction correction input data')
        
#    def _plProcess(self):
#        from .utils import  PipelineError
#        if self._plbuffer is None:
#            self._plbuffer = {'in':[],'out':[]} 
        
#        type,data = self._plbuffer['in'].pop(0)
#        try:
#            if type=='spec':
#                newspec = self.correctSpectrum(data)
#                self._plbuffer['out'].append(newspec)
#            else:
#                assert False,'Impossible point - code error in Extinction pipeline'
#        except:
#            #TODO:check this
#            self.insert(0,spec(type,data))
            
#    def _plExtract(self):
#        if self._plbuffer is None:
#            self._plbuffer = {'in':[],'out':[]} 
            
#        if len(self._plbuffer['out']<1):
#            return None
#        else:
#            return self._plbuffer['out'].pop(0)
        
#    def plClear(self):
#        self._plbuffer = None
    
class CalzettiExtinction(Extinction):
    """
    This is the average extinction law derived in Calzetti et al. 1994:
    x=1/lambda in mu^-1
    Q(x)=-2.156+1.509*x-0.198*x**2+0.011*x**3
    """
    _poly=np.poly1d((0.011,-0.198,1.509,-2.156)) 
    def __init__(self,A0=1):
        super(CalzettiExtinction,self).__init__(A0=A0)
        
    def f(self,lamb):
        return (-2.5/np.log(10))*self._poly(1e4/lamb)
    
class _EBmVExtinction(Extinction):
    """
    Base class for Extinction classes that get normalization from E(B-V)
    """
    
    def __init__(self,EBmV=1,Rv=3.1):
        super(_EBmVExtinction,self).__init__(f=None,A0=1)
        self.Rv=Rv
        self.EBmV=EBmV
    
    def _getEBmV(self):
        from .phot import bandwl
        return self.A0*self.f(bandwl['V'])/self.Rv
    def _setEBmV(self,val):
        from .phot import bandwl
        Av=self.Rv*val
        self.A0=Av/self.f(bandwl['V'])
    EBmV = property(_getEBmV,_setEBmV)
    
    def f(self,lamb):
        raise NotImplementedError
    
class FMExtinction(_EBmVExtinction):
    """
    Base class for Extinction classes that use the form from 
    Fitzpatrick & Massa 90
    """
    
    def __init__(self,C1,C2,C3,C4,x0,gamma,EBmV=1,Rv=3.1):
        self.x0 = x0
        self.gamma = gamma
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
        
        super(FMExtinction,self).__init__(EBmV=EBmV,Rv=Rv)
        
    def f(self,lamb):
        x=1e4/np.array(lamb,copy=False)
        C1,C2,C3,C4 = self.C1,self.C2,self.C3,self.C4
        gamma,x0 = self.gamma,self.x0
        
        xsq=x*x
        D=xsq*((xsq-x0*x0)**2+xsq*gamma*gamma)**-2
        FMf=C1+C2*x+C3*D
        
        if np.isscalar(FMf):
            if x>=5.9:
                FMf+=C4*(0.5392*(x-5.9)**2+0.05644*(x-5.9)**3)
        else:
            C4m=x>=5.9
            FMf[C4m]+=C4*(0.5392*(x[C4m]-5.9)**2+0.05644*(x[C4m]-5.9)**3)
        
        return FMf+self.Rv #EBmV is the normalization and is multiplied in at the end
    
class CardelliExtinction(_EBmVExtinction):
    """
    Milky Way Extinction law from Cardelli et al. 1989
    """
    def f(self,lamb):
        scalar=np.isscalar(lamb)
        x=1e4/np.array(lamb,ndmin=1) #CCM x is 1/microns
        a,b=np.ndarray(x.shape,x.dtype),np.ndarray(x.shape,x.dtype)
        
        if any((x<0.3)|(10<x)):
            raise ValueError('some wavelengths outside CCM 89 extinction curve range')
        
        irs=(0.3 <= x) & (x <= 1.1)
        opts = (1.1 <= x) & (x <= 3.3)
        nuv1s = (3.3 <= x) & (x <= 5.9)
        nuv2s = (5.9 <= x) & (x <= 8)
        fuvs = (8 <= x) & (x <= 10)
        
        #TODO:pre-compute polys
        
        #CCM Infrared
        a[irs]=.574*x[irs]**1.61
        b[irs]=-0.527*x[irs]**1.61
        
        #CCM NIR/optical
        a[opts]=np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
        b[opts]=np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)
        
        #CCM NUV
        y=x[nuv1s]-5.9
        Fa=-.04473*y**2-.009779*y**3
        Fb=-.2130*y**2-.1207*y**3
        a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)+Fa
        b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)+Fb
        
        a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)
        b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)
        
        #CCM FUV
        a[fuvs]=np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
        b[fuvs]=np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)
        
        AloAv = a+b/self.Rv
        
        if scalar:
            return AloAv[0]
        else:
            return AloAv
    
class LMCExtinction(FMExtinction):
    """
    LMC Extinction law from Gordon et al. 2003 LMC Average Sample
    """
    def __init__(self,EBmV=.3,Rv=3.41):
        super(LMCExtinction,self).__init__(-.890,0.998,2.719,0.400,4.579,0.934,EBmV,Rv)
    
class SMCExtinction(FMExtinction):
    """
    SMC Extinction law from Gordon et al. 2003 SMC Bar Sample
    """
    def __init__(self,EBmV=.2,Rv=2.74):
        super(SMCExtinction,self).__init__(-4.959,2.264,0.389,0.461,4.6,1,EBmV,Rv)

def get_SFD_dust(long,lat,dustmap='ebv',interpolate=True):
    """
    Gets map values from Schlegel, Finkbeiner, and Davis 1998 extinction maps.
    
    `dustmap` can either be a filename (if '%s' appears in the string, it will be
    replaced with 'ngp' or 'sgp'), or one of:
    
    * 'i100' 
        100-micron map in MJy/Sr
    * 'x'
        X-map, temperature-correction factor
    * 't'
        Temperature map in degrees Kelvin for n=2 emissivity
    * 'ebv'
        E(B-V) in magnitudes
    * 'mask'
        Mask values 
        
    For these forms, the files are assumed to lie in the current directory.
    
    Input coordinates are in degrees of galactic latiude and logitude - they can
    be scalars or arrays.
    
    if `interpolate` is an integer, it can be used to specify the order of the
    interpolating polynomial
    """
    from numpy import sin,cos,round,isscalar,array,ndarray,ones_like
    from pyfits import open
    
    if type(dustmap) is not str:
        raise ValueError('dustmap is not a string')
    dml=dustmap.lower()
    if dml == 'ebv' or dml == 'eb-v' or dml == 'e(b-v)' :
        dustmapfn='SFD_dust_4096_%s.fits'
    elif dml == 'i100':
        dustmapfn='SFD_i100_4096_%s.fits'
    elif dml == 'x':
        dustmapfn='SFD_xmap_%s.fits'
    elif dml == 't':
        dustmapfn='SFD_temp_%s.fits'
    elif dml == 'mask':
        dustmapfn='SFD_mask_4096_%s.fits'
    else:
        dustmapfn=dustmap
    
    if isscalar(long):
        l=array([long])*pi/180
    else:
        l=array(long)*pi/180
    if isscalar(lat):
        b=array([lat])*pi/180
    else:
        b=array(lat)*pi/180
        
    if not len(l)==len(b):
        raise ValueError('input coordinate arrays are of different length')
    
    
    
    if '%s' not in dustmapfn:
        f=open(dustmapfn)
        try:
            mapds=[f[0].data]
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        
        polename=dustmapfn.split('.')[0].split('_')[-1].lower()
        if polename=='ngp':
            n=[1]
            if sum(b > 0) > 0:
                print 'using ngp file when lat < 0 present... put %s wherever "ngp" or "sgp" should go in filename'
        elif polename=='sgp':
            n=[-1]
            if sum(b < 0) > 0:
                print 'using sgp file when lat > 0 present... put %s wherever "ngp" or "sgp" should go in filename'
        else:
            raise ValueError("couldn't determine South/North from filename - should have 'sgp' or 'ngp in it somewhere")
        masks = [ones_like(b).astype(bool)]
    else: #need to do things seperately for north and south files
        nmask = b >= 0
        smask = ~nmask
        
        masks = [nmask,smask]
        ns = [1,-1]
        
        mapds=[]
        f=open(dustmapfn%'ngp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
        f=open(dustmapfn%'sgp')
        try:
            mapds.append(f[0].data)
        finally:
            f.close()
        assert mapds[-1].shape[0] == mapds[-1].shape[1],'map dimensions not equal - incorrect map file?'
    
    retvals=[]
    for n,mapd,m in zip(ns,mapds,masks):
        #project from galactic longitude/latitude to lambert pixels (see SFD98)
        npix=mapd.shape[0]
        
        x=npix/2*cos(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        y=-npix/2*n*sin(l[m])*(1-n*sin(b[m]))**0.5+npix/2-0.5
        #now remap indecies - numpy arrays have y and x convention switched from SFD98 appendix
        x,y=y,x
        
        if interpolate:
            from scipy.ndimage import map_coordinates
            if type(interpolate) is int:
                retvals.append(map_coordinates(mapd,[x,y],order=interpolate))
            else:
                retvals.append(map_coordinates(mapd,[x,y]))
        else:
            x=round(x).astype(int)
            y=round(y).astype(int)
            retvals.append(mapd[x,y])
            
            
    
        
    if isscalar(long) or isscalar(lat):
        for r in retvals:
            if len(r)>0:
                return r[0]
        assert False,'None of the return value arrays were populated - incorrect inputs?'
    else:
        #now recombine the possibly two arrays from above into one that looks like  the original
        retval=ndarray(l.shape)
        for m,val in zip(masks,retvals):
            retval[m] = val
        return retval
        
    
def get_dust_radec(ra,dec,dustmap,interpolate=True):
    from .coords import equatorial_to_galactic
    l,b = equatorial_to_galactic(ra,dec)
    return get_SFD_dust(l,b,dustmap,interpolate)


  
#DEPRECATED!
def extinction_correction(lineflux,linewl,EBmV,Rv=3.1,exttype='MW'):
    """
    
    Extinction correct a la Cardelli et al 89 from the supplied line data and
    a given E(B-V) along the sightline 
    
    inputs may be numpy arrays
    
    linewl is in angstroms, lineflux in erg s^-1 cm^-2
    
    if the input lineflux is None (or NaN but NOT False or 0) , Alambda is returned instead
    """
    from numpy import array,logical_and,logical_or,polyval,log10,isscalar,where
    from .phot import bandwl
    
    from warnings import warn
    warn('extinction_correction function is deprecated - use Extinction class instead',DeprecationWarning)
    
    if exttype=='LMC':
        eo = LMCExtinction(EBmV=EBmV,Rv=Rv)
        return eo.correctFlux(lineflux,linewl)
    elif exttype == 'SMC':
        eo = SMCExtinction(EBmV=EBmV,Rv=Rv)
        return eo.correctFlux(lineflux,linewl)
    
    if isinstance(linewl,basestring):
        linewl=bandwl[linewl]
        
    if lineflux is None or lineflux is False:
        lineflux=np.nan
        
    
    lineflux=np.array(lineflux,ndmin=1,dtype=float)
    linewl=np.array(linewl,ndmin=1)
    EBmV=np.array(EBmV,ndmin=1)
    n=np.max((lineflux.size,linewl.size,EBmV.size))
    
    if n!=1:
        if lineflux.size == 1:
            lineflux = np.ones(n)*lineflux[0]
        elif lineflux.size != n:
            raise ValueError("lineflux is incorrect length")
        if linewl.size == 1:
            linewl = np.ones(n)*linewl[0]
        elif linewl.size != n:
            raise ValueError("linewl is incorrect length")
        if EBmV.size == 1:
            EBmV = np.ones(n)*EBmV[0]
        elif EBmV.size != n:
            raise ValueError("EBmV is incorrect length")
        
    x=1e4/linewl #CCM x is 1/microns
    
    a=array(x)
    b=array(x)
    
    if any(logical_or(x<0.3,10<x)):
        raise ValueError('some wavelengths outside CCM 89 extinction curve range')
    
    irs=where(logical_and(0.3 <= x,x <= 1.1))
    opts=where(logical_and(1.1 <= x,x <= 3.3))
    nuv1s=where(logical_and(3.3 <= x,x <= 5.9))
    nuv2s=where(logical_and(5.9 <= x,x <= 8))
    fuvs=where(logical_and(8 <= x,x <= 10))
    
    #CCM Infrared
    a[irs]=.574*x[irs]**1.61
    b[irs]=-0.527*x[irs]**1.61
    
    #CCM NIR/optical
    a[opts]=polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
    b[opts]=polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)
    
    #CCM NUV
    y=x[nuv1s]-5.9
    Fa=-.04473*y**2-.009779*y**3
    Fb=-.2130*y**2-.1207*y**3
    a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)+Fa
    b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)+Fb
    
    a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)
    b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)
    
    #CCM FUV
    a[fuvs]=polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
    b[fuvs]=polyval((.374,-.42,4.257,13.67),x[fuvs]-8)
    
    AloAv=a+b/Rv #Al/Av
    Al=AloAv*Rv*EBmV #from Rv=Av/E(B-V)
    
    
    finalval=Al
    magi= ~np.isnan(lineflux)
    realmag=-2.5*log10(lineflux[magi])-Al[magi]
    finalval[magi]=10**(realmag/-2.5)
    
    if n==1:
        return finalval[0]
    else:
        return finalval
        
#DEPRECATED!
def extinction_from_flux_ratio(frobs,frexpect,outlambda=None,Rv=3.1,tol=1e-4):
    """
    frobs is the observed flux ratio f1/f2
    
    frexpect is either a string code specifying a hydrogen transition
    (e.g. 'Hab' is Halpha/Hbeta, 'Hde' is Hdelta/Hepsilon, from Ostriker 2E)    
    or a tuple of the form (expected f1/f2,lambda1,lambda2) wl in angstroms
    
    outputs E(B-V) to a tolerance specified by tol if outlambda is 0/False/None,
    otherwise outputs Alambda (tol still determines E(B-V)) (outlambda can be
    UBVRI or ugriz as from B&M)
    
    frobs can be an array, but other values cannot
    """
    from scipy.optimize import fmin
    
    from warnings import warn
    warn('extinction_from_flux_ratio function is deprecated - use Extinct class instead',DeprecationWarning)
    
    scalarout=np.isscalar(frobs)
    
    frobs=np.array(frobs,ndmin=1)
    
    hd={
    'Hab':(2.86,6562.82,4861.33),
    'Hag':(6.16,6562.82,4340.46),
    'Had':(11.21,6562.82,4101.74),
    'Hae':(18.16,6562.82,3970.07),
    'Hbg':(2.15,4861.33,4340.46),
    'Hbd':(3.91,4861.33,4101.74),
    'Hbe':(6.33,4861.33,3970.07),
    'Hgd':(1.82,4340.46,4101.74),
    'Hge':(2.95,4340.46,3970.07),
    'Hde':(1.62,4101.74,3970.07)
    }
    if isinstance(frexpect,basestring):
        frexpect=hd[frexpect]
    
    fr0,l1,l2=frexpect
    
    A2mA1=-2.5*np.log10(fr0/frobs)
    
    fres = lambda x,Y: np.abs(Y - (extinction_correction(None,l2,x,Rv)-extinction_correction(None,l1,x,Rv)))
    
    EBmV=np.ndarray(A2mA1.shape)
    EBmVr=EBmV.ravel()
    for i,A in enumerate(A2mA1.ravel()):
        EBmVr[i]=fmin(fres,0,args=(A,),xtol=tol,disp=0)[0]
        
    if scalarout:
        EBmV=EBmV[0]
        
    if outlambda:
        if isinstance(outlambda,basestring):
            from .phot import bandwl
            outlambda=bandwl[outlambda]
        return extinction_correction(None,outlambda,EBmV,Rv)
    else:
        return EBmV
    

