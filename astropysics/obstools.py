#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu) 

"""
This module stores tools for oberving (pre- and post-) as well as functioning as
a module for various corrections and calculations that don't have a better place
to live.

The focus is currently on optical astronomy, as that is what the primary author
does.

Note that some of these functions require the
:mod:`dateutil<http://pypi.python.org/pypi/python-dateutil>` package (it is
included with matplotlib)
"""
#TODO: add info to Observatory class, particularly atmospheric extinction
#TODO: exposure time calculator (maybe in phot instead?)
#TODO: make Extinction classes spec.HasSpecUnits and follow models framework?
#TODO: return Pipeline support to Extinction 

from __future__ import division,with_statement
from .constants import pi
import numpy as np

from .utils import PipelineElement,DataObjectRegistry
#from .models import FunctionModel1DAuto as _FunctionModel1DAuto

def jd_to_gregorian(jd,bceaction=None,rounding=1e-5):
    """
    Convert a julian date to a gregorian calendar date and time.
    
    
    `bceaction` indicates what to do if the result is a BCE year.   datetime 
    generally only supports positive years, so the following options apply:
    
    * 'raise' or None
        raise an exception
    * 'neg'
        for any BCE years, convert to positive years (delta operations in
        datetime will not be correct)
    * 'negreturn'
        same as 'neg', but changes the return type to a tuple with the datetime
        objects as the first element and a boolean array that is True for the
        objects that are BCE as the second element.
    * a scalar
        add this number to the year
    
    
    `rounding` determines a fix for floating-point errors - if the seconds are
    within the request number of seconds of an exact second, the output will be
    set to that second
    
    returns :class:`datetime.datetime` objects in UTC. If `bceaction` is
    'negreturn', it will instead be a (datetime,bcemask) tuple
    """
    import datetime
    
    jdn = np.array(jd,copy=False)
    scalar = jdn.shape == ()
    jdn = jdn.ravel()
    
    #implementation from http://www.astro.uu.nl/~strous/AA/en/reken/juliaansedag.html
    x2 = np.floor(jdn - 1721119.5)
    c2 = (4*x2 + 3)//146097
    x1 = x2 - (146097*c2)//4
    c1 = (100*x1 + 99)//36525
    x0 = x1 - (36525*c1)//100
    
    yr = 100*c2 + c1
    mon = (5*x0 + 461)//153
    day = x0 - (153*mon-457)//5 + 1
    
    mmask = mon>12
    mon[mmask] -= 12
    yr[mmask] += 1
    
    
    time = (jdn - np.floor(jdn) - 0.5) % 1
    if time == 0:
        hr = min = sec = msec = np.zeros_like(yr)
    elif time == 0.5:
        min = sec = msec = np.zeros_like(yr)
        hr = 12*np.ones_like(yr)
    else:
        hr = np.floor(time*24).astype(int)
        time = time*24 - hr
        min = np.floor(time*60).astype(int)
        time = time*60 - min
        sec = np.floor(time*60).astype(int)
        time = time*60 - sec
        msec = np.floor(time*1e6).astype(int)
        
        #rounding fix for floating point errors
        if rounding:
            msecundermask = (1-time) < rounding
            msecovermask = time < rounding
            msec[msecundermask|msecovermask] = 0
            sec[msecundermask] += 1
            min[sec==60] += 1
            sec[sec==60] = 0
            hr[min==60] += 1
            min[min==60] = 0
            hr[hr==24] = 0
            #date is always right
        
    
    try:
        from dateutil import tz
        tzi = tz.tzutc()
    except ImportError:
        tzi = None
    
    if bceaction is not None:
        if bceaction == 'raise':
            pass
        elif bceaction == 'neg' or bceaction == 'return':
            bcemask = yr<datetime.MINYEAR
            yr[bcemask] -= 1 #correct by 1 because JD=0 -> 1 BCE
            yr[bcemask]*= -1 
                
        else:
            yr += bceaction
    
    if scalar:
        res = datetime.datetime(yr[0],mon[0],day[0],hr[0],min[0],sec[0],msec[0],tzi)
    else:
        res = [datetime.datetime(*t,tzinfo=tzi) for t in zip(yr,mon,day,hr,min,sec,msec)]
        
    if bceaction == 'return':
        return res,bcemask
    else:
        return res
    

def gregorian_to_jd(gtime,tz=None):
    """
    Convert gregorian calendar value to julian date
    
    the input `gtime` can either be a sequence (yr,month,day,[hr,min,sec]), 
    where each element may be  or a :class:`datetime.datetime` object
    
    if datetime objects are given, and `tz` is None, values are converted to
    UTC based on the datetime.tzinfo objects (if present).  if `tz` is not None,
    the tzinfo in the datetime objects is ignored.
    
    If `tz` is a string, it is taken to be a timezone name that will be used
    to convert all the dates to UTC (requires :mod:`dateutil` package).  If `tz`
    it is a scalar, it is taken to be the hour offset of the timezone. 
    Or, if it is a :class:`datetime.tzinfo` object, that object will be used to
    do the conversion to UTC.  
    """
    #Adapted from xidl  jdcnv.pro
    from datetime import datetime,date,tzinfo
    
    if isinstance(gtime,datetime) or isinstance(gtime,date):
        datetimes = [gtime]
        scalarout = True
    elif all([isinstance(gt,datetime) for gt in gtime]):
        datetimes = gtime
        scalarout = False
    else:
        datetimes = None
        gtime = list(gtime)
        if not (3 <= len(gtime) < 7):
            raise ValueError('gtime input sequence is invalid size')
        while len(gtime) < 6:
            gtime.append(np.zeros_like(gtime[-1]))
        yr,month,day,hr,min,sec = gtime
        scalarout = False #already a scalar form
        
    #if input objects are datetime objects, generate arrays
    if datetimes is not None:
        yr,month,day,hr,min,sec = [],[],[],[],[],[]
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
            sec.append(dt.second+dt.microsecond/1e6)
                
                
    
    yr = np.array(yr,dtype='int64',copy=False).ravel()
    month = np.array(month,dtype='int64',copy=False).ravel()
    day = np.array(day,dtype='int64',copy=False).ravel()
    hr = np.array(hr,dtype=float,copy=False).ravel()
    min = np.array(min,dtype=float,copy=False).ravel()
    sec = np.array(sec,dtype=float,copy=False).ravel()
    
    #do tz conversion if tz is provided  
    if isinstance(tz,basestring) or isinstance(tz,tzinfo):
        if isinstance(tz,basestring):
            from dateutil import tz
            tzi = tz.gettz(tz)
        else:
            tzi = tz
        
        utcoffset = []
        for t in zip(yr,month,day,hr,min,sec):
            #microsecond from float component of seconds
            t = list(t)
            t.append(int((t[-1]-np.floor(t[-1]))*1e6))
            dt = datetime(*t,tzinfo=tzi)
            utcdt = dt.utcoffset()
            if utcdt is None:
                utcoffset.append(0)
            else:
                utcoffset.append(utcdt.days*24 + (utcdt.seconds + utcdt.microseconds*1e-6)/3600)
    else:
        utcoffset = tz
            
    ly = int((month-14)/12)		#In leap years, -1 for Jan, Feb, else 0
    jdn = day - 32075l + 1461l*(yr+4800l+ly)//4
    jdn += 367l*(month - 2-ly*12)//12 - 3*((yr+4900l+ly)//100)//4
    
    res = jdn + (hr/24.0) + min/1440.0 + sec/86400.0 - 0.5
    
    if np.any(utcoffset):
        res -= np.array(utcoffset)/24.0
    
    if scalarout:
        return res[0]
    else:
        return res
    
    
def besselian_epoch_to_jd(bepoch):
    """
    Convert a Besselian epoch to Julian Date, assuming a tropical year of 
    365.242198781 days.
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epb.html
    """
    return (bepoch - 1900)*365.242198781 + 2415020.31352

def jd_to_besselian_epoch(jd):
    """
    Convert a Julian Date to a Besselian epoch, assuming a tropical year of 
    365.242198781 days
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epb.html
    """
    return 1900 + (jd - 2415020.31352)/365.242198781

def jd_to_epoch(jd):
    """
    Convert a Julian Date to a Julian Epoch, assuming the year is exactly
    365.25 days long
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epj.html
    """
    return 2000.0 + (jd - 2451545.0)/365.25

def epoch_to_jd(jepoch):
    """
    Convert a Julian Epoch to a Julian Day, assuming the year is exactly
    365.25 days long
    
    :Reference: http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/epj.html
    """
    return (jepoch - 2000)*365.25 + 2451545.0
    
    

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
        elif isinstance(tz,basestring):
            from dateutil import tz  as tzmod
            self.tz = tzmod.gettz(tz)
            if self.tz is None:
                raise ValueError('unrecognized time zone string '+tz)
        else:
            from dateutil import tz as tzmod
            self.tz = tzmod.tzoffset(str(tz),int(tz*60*60))
        if name is None:
            name = 'Default Site'
        self.name = name
        self._currjd = None
    
    @staticmethod
    def _tzFromLong(long):
        from dateutil import tz
        offset =  int(np.round(long.degrees/15))
        return tz.tzoffset(str(offset),offset*60*60)
    
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
            return gregorian_to_jd(datetime.utcnow(),tz=None)
        else:
            return self._currjd
    def _setCurrentobsjd(self,val):
        if val is None:
            self._currjd = None
        else:
            if np.isscalar(val):
                self._currjd = val
            else:
                self._currjd = gregorian_to_jd(val)
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
                jd = gregorian_to_jd(dtobj,tz=None)
            else:
                jd = args[0]
        elif len(args) == 4:
            time,year,month,day = args
            hr = int(np.floor(time))
            min = int(np.floor(60*(time - hr)))
            sec = int(np.floor(60*(60*(time-hr) - min)))
            msec = int(np.floor(1e6*(60*(60*(time-hr) - min) - sec)))
            jd = gregorian_to_jd(datetime.datetime(year,month,day,hr,min,sec,msec,self.tz),tz=None)
        elif len(args) == 6:
            year,month,day,hr,min,sec = args
            msec = int(1e6*(sec - np.floor(sec)))
            sec = int(np.floor(sec))
            jd = gregorian_to_jd(datetime.datetime(year,month,day,hr,min,sec,msec,self.tz),tz=None)
        else:
            raise TypeError('invalid number of input arguments')
        
        #algorithm described on USNO web site http://aa.usno.navy.mil/faq/docs/GAST.php
        jd0 = np.round(jd-.5)+.5
        h = (jd - jd0) * 24.0
        d = jd - 2451545.0
        d0 = jd0 - 2451545.0
        t = d/36525
        
        gmst = 6.697374558 + 0.06570982441908*d0 + 1.00273790935*h + 0.000026*t**2
       
        if apparent:
            eps =  np.radians(23.4393 - 0.0000004*d) #obliquity
            L = np.radians(280.47 + 0.98565*d) #mean longitude of the sun
            omega = np.radians(125.04 - 0.052954*d) #longitude ofascending node of moon
            dpsi = -0.000319*np.sin(omega) - 0.000024*np.sin(2*L) #nutation longitude
            lst = (gmst + dpsi*np.cos(eps) + self._long.d/15)%24.0
        else:
            lst = (gmst + self._long.d/15)%24.0 
        
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
        
    def localTime(self,lsts,date=None,apparent=True,returntype=None):
        """
        Computes the local civil time given a particular local siderial time. 
        
        `lsts` are the input local times, and may be either a scalar or array,
        and must be in decimal hours.
        
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
            a :class:`datetime.time` object with the local tzinfo
        
        """
        import datetime
        
        if isinstance(lsts,datetime.datetime):
            date = lsts.date()
            lsts = datetime.time()
            lsts = lsts.hour+lsts.minute/60+(lsts.second+lsts.microsecond*1e6)/3600
        else:
            if date is None:
                date = jd_to_gregorian(self.currentobsjd).date()
            elif hasattr(date,'date'):
                date = date.date()
            else:
                date = datetime.date(*date)
                
        lsts = np.array(lsts,copy=False)
        scalarout = False
        if lsts.shape == tuple():
            scalarout = True
            lsts = lsts.ravel()
        
            
        lst0 = self.localSiderialTime(date)
        lthrs = (lsts - lst0)%24
        dayoffs = np.floor(lsts - lst0/24)
        
        lthrs /= 1.0027378507871321 #correct for siderial day != civil day
        
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
        horizontal coordinates directly - :meth:`apparentPosition` includes all 
        the corrections and formatting.  This method is purely for the 
        coordinate conversion.
        """
        from .coords import EquatorialPosition,HorizontalPosition
        
        latitude = self.latitude.d
        
        if epoch is not None:
            #make copy so as not to change the epoch of the input
            eqpos = EquatorialPosition(eqpos)
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
                return HorizontalPosition(*np.degrees((alts[0],azs[0])))
            else:
                return HorizontalPosition(*np.degrees((alts[0],azs[0],dazs[0],dalts[0])))
        else:
            if dazs is None:
                return [HorizontalPosition(alt,az) for alt,az in np.degrees((alts,azs)).T]
            else:
                return [HorizontalPosition(alt,az,daz,dalt) for alt,az,daz,dalt in np.degrees((alts,azs,dazs,dalts)).T]
            
    def riseSetTransit(self,eqpos,alt=0,timeobj=False):
        """
        computes the rise, set, and transit times of a provided equatorial 
        position.  
        
        `alt` determines the altitude to be considered as risen or set in 
        degrees
        
        *returns*
        (rise,set,transit) as :class:datetime.time objects if `timeobj` is True
        or if it is False, they are in decimal hours.  If the object is 
        circumpolar, rise and set are both None.  If it is never visible,
        rise,set, and transit are all None
        """
        import datetime
        from math import sin,cos,radians,acos
        
        alt = radians(alt)
        
        transit = self.localTime(eqpos.ra.hours)
        
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
                    return datetime.time(int(hr),int(min),int(sec),int(msec))
            return hrtotime(rise),hrtotime(set),hrtotime(transit)
        else:
            return rise,set,transit
        
        
    def apparentPosition(self,coords,datetime=None,precess=True):
        """
        computes the positions in horizontal coordinates of an object with the 
        provided fixed coordinates at the requested time(s).
        
        `datetime` can be a sequence of datetime objects or similar representations,
        or it can be a sequnce of JD's.  If None, it uses the current time or 
        the value of the :attr:`currentobsjd` property.
        
        If `precess` is True, the position will be precessed to the epoch of 
        the observation
        
        
        *returns*
        a sequence of :class:`astropysics.coords.HorizontalPosition` objects
        """
        from .coords import HorizontalPosition
        from operator import isSequenceType
        
        if datetime is None:
            return self.equatorialToHorizontal(coords,self.localSiderialTime())
        elif hasattr(datetime,'year') or (isSequenceType(datetime) and hasattr(datetime[0],'year')):
            jd = gregorian_to_jd(datetime,self.tz).ravel()
        else:
            jd = np.array(datetime,copy=False)
            if len(jd.shape)>1:
                jd = np.array([gregorian_to_jd(v,self.tz) for v in jd])
        lsts = self.localSiderialTime(jd)
        
        if precess:
            return self.equatorialToHorizontal(coords,lsts,epoch=jd_to_epoch(jd[0]))
        else:
            return self.equatorialToHorizontal(coords,lsts)
        
    def observingTable(self,coord,date=None,strtablename=None,hrrange=(18,6,25)):
        """
        tabulates altitude, azimuth, and airmass values for the provided fixed
        position on a particular date, specified either as a datetime.date 
        object, a (yr,mon,day) tuple, or a julain date (will be rounded to 
        nearest)
        
        if `date` is None, the current date is used
        
        `hrrange` determines the size of the table - it should be a 3-tuple
        (starthr,endhr,n) where starthr is on the day specified in date and 
        endhr is date + 1 day
        
        if `strtablename` is True, a string is returned with a printable table 
        of observing data.  If `strtable` is a string, it will be used as the 
        title for the table.  Otherwise, a record array is returned with the 
        hour, alt, az, and airmass
        """        
        from .astropysics import coords
        import datetime
        
        if date is None:
            jd = self.currentobsjd
        elif np.isscalar(date):
            jd = date
        else:
            jd = gregorian_to_jd(date)
            
        date = jd_to_gregorian(jd).date()
        jd0 = np.floor(jd)
        
        starthr,endhr,n = hrrange
        startjd = jd0 + (starthr - 12)/24
        endjd = jd0 + (endhr + 12)/24
        jds = np.linspace(startjd,endjd,n)
        
        timehr = (jds-np.round(np.mean(jds))+.5)*24
        
        alt,az = coords.objects_to_coordinate_arrays(self.apparentPosition(coord,jds,precess=False))
        
        z = 90 - alt
        airmass = 1/np.cos(np.radians(z))
        
        ra = np.rec.fromarrays((timehr,alt,az,airmass),names = 'hour,alt,az,airmass')

        
        if strtablename is not None:
            rise,set,transit = self.riseSetTransit(coord,timeobj=True)
            
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
            rise,set,transit = self.riseSetTransit(coord,timeobj=False)
            ra._sitedate = date
            ra.rise = rise
            ra.set = set
            ra.transit = transit
            return ra
        
    def observingPlot(self,coords,date=None,names=None,clf=True,plottype='altam',
                           plotkwargs=None):
        """
        generates plots of important observability quantities for the provided
        coordinates.
        
        `coords` should be an :class:`astropysics.coords.LatLongPosition`
        
        `plottype` can be one of:
        
        * 'altam' : a plot of time vs. altitude with a secondary axis for sec(z)
        * 'am' : a plot of time vs. airmass/sec(z)
        * 'altaz': a plot of azimuth vs. altitude
        * 'sky': polar projection of the path on the sky
        
        """
        import matplotlib.pyplot as plt
        from operator import isMappingType
        from .coords import LatLongPosition
        from .plotting import add_mapped_axis
        
        if isinstance(coords,LatLongPosition):
            coords = [coords]
            
        nonames = False
        if names is not None:
            if isinstance(names,basestring):
                names = [names]
            if len(names) != len(coords):
                raise ValueError('names do not match coords')
        else:
            nonames = True
            names = ['' for c in coords]
            
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
            if plottype == 'altam' or plottype == 'alt':   
                if plottype == 'altam':
                    oldright = plt.gcf().subplotpars.right
                    plt.subplots_adjust(right=0.86)
                    
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.observingTable(c,date,hrrange=(0,0,100))
                    if kw is None:
                        plt.plot(ra.hour,ra.alt,label=n)
                    else:
                        kw['label'] = n
                        plt.plot(ra.hour,ra.alt,**kw)
                    
                plt.xlim(0,24)
                
                if plt.ylim()[0] < 0:
                    lowery = 0
                else:
                    lowery = plt.ylim()[0]
                uppery = plt.ylim()[1] if plt.ylim()[1]<90 else 90
                plt.ylim(lowery,uppery)
                
                if not nonames:
                    plt.legend(loc=0)
                
                plt.title(str(ra._sitedate))
                plt.xlabel(r'$\rm hours$')
                plt.ylabel(r'${\rm altitude} [{\rm degrees}]$')
                if plottype == 'altam':
                    yticks = plt.yticks()[0]
                    newticks = list(yticks[1:])
                    newticks.insert(0,(yticks[0]+newticks[0])/2)
                    newticks.insert(0,(yticks[0]+newticks[0])/2)
                    plt.twinx()
                    plt.ylim(lowery,uppery)
                    plt.yticks(newticks,['%2.2f'%(1/np.cos(np.radians(90-yt))) for yt in newticks])
                    plt.ylabel(r'$\sec(z)$')
                    
                plt.xticks(np.arange(13)*2)
                    
            elif plottype == 'am':
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.observingTable(c,date,hrrange=(0,0,100))
                    
                    ammask = ra.airmass>0
                    if kw is None:
                        plt.plot(ra.hour[ammask],ra.airmass[ammask],label=n)
                    else:
                        kw['label'] = n
                        plt.plot(ra.hour[ammask],ra.airmass[ammask],**kw)
                        
                plt.xlim(0,24)
                plt.xticks(np.arange(13)*2)
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
                    plt.ylim(lowery,uppery)
                
                plt.title(str(ra._sitedate))
                plt.xlabel(r'$\rm hours$')
                plt.ylabel(r'${\rm airmass} / \sec(z)$')
                
                if not nonames:
                    plt.legend(loc=0)
                    
            elif plottype == 'altaz':
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.observingTable(c,date,hrrange=(0,0,100))
                    
                    if kw is None:
                        plt.plot(ra.az,ra.alt,label=n)
                    else:
                        kw['label'] = n
                        plt.plot(ra.az,ra.alt,**kw)
                        
                plt.xlim(0,360)
                plt.xticks(np.arange(9)*360/8)
                plt.ylim(0,90)
                
                plt.title(str(ra._sitedate))
                plt.xlabel(r'${\rm azimuth} [{\rm degrees}]$')
                plt.ylabel(r'${\rm altitude} [{\rm degrees}]$')
                
                if not nonames:
                    plt.legend(loc=0)
                
            elif plottype == 'sky':
                for n,c,kw in zip(names,coords,plotkwargs):
                    ra = self.observingTable(c,date,hrrange=(0,0,100))
                    
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
                
                plt.title(str(ra._sitedate))
                        
                if not nonames:
                    plt.legend(loc=0)
            else:
                raise ValueError('unrecognized plottype {0}'.format(plottype))
            
            plt.show()
            plt.draw()
            if oldright is not None:
                plt.gcf().subplotpars.right = oldright
        finally:
            plt.interactive(inter)
        
        
        
class Observatory(Site):
    """
    Represents an observatory/Telescope with associated information such as
    platform limits or extinction measurements
    """
    


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
                o = Observatory(lat,long,alt,tz,name)
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
        o = Observatory(lat,long,alt,tz,name)
        sitereg[obs] = o


sites = DataObjectRegistry('sites',Site)
sites['uciobs'] = Observatory(33.63614044191056,-117.83079922199249,80,'PST','UC Irvine Observatory')
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
    from .phot import bandwl
    __lambdaV=bandwl['V']
    del bandwl
    
    def __init__(self,EBmV=1,Rv=3.1):
        super(_EBmVExtinction,self).__init__(f=None,A0=1)
        self.Rv=Rv
        self.EBmV=EBmV
    
    def _getEBmV(self):
        return self.A0*self.f(self.__lambdaV)/self.Rv
    def _setEBmV(self,val):
        Av=self.Rv*val
        self.A0=Av/self.f(self.__lambdaV)
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
        x=1e4/np.atleast_1d(lamb) #CCM x is 1/microns
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
        
    
    lineflux=np.atleast_1d(lineflux).astype(float)
    linewl=np.atleast_1d(linewl)
    EBmV=np.atleast_1d(EBmV)
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
    
    frobs=np.atleast_1d(frobs)
    
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
    
