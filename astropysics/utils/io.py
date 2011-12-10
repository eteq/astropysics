#Copyright 2009 Erik Tollerud
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
The :mod:`io` module contains classes and functions for loading and saving data
in various relevant formats used in astronomy, as well as convinience functions
for retrieval of built-in data.

.. todo:: examples


Classes and Inheritance Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. inheritance-diagram:: astropysics.utils.io
   :parts: 1

Module API
^^^^^^^^^^

"""

from __future__ import division,with_statement
from .gen import add_docs_and_sig as _add_docs_and_sig
from .gen import add_docs as _add_docs
from ..config import get_config
import numpy as np

_io_config = get_config('io')

def funpickle(fileorname,number=0,usecPickle=True):
    """
    Unpickle a pickled object from a specified file and return the contents.

    :param fileorname: The file from which to unpickle objects
    :type fileorname: a file name string or a :class:`file` object
    :param number:
        The number of objects to unpickle - if <1, returns a single object.
    :type number: int
    :param usecPickle:
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster).
    :type usecPickle: bool

    :returns: A list of length given by `number` or a single object if number<1

    """
    if usecPickle:
        import cPickle as pickle
    else:
        import pickle

    if isinstance(fileorname,basestring):
        f = open(fileorname,'r')
        close = True
    else:
        f = fileorname
        close = False

    try:
        if number > 0:
            res = []
            for i in range(number):
                res.append(pickle.load(f))
        elif number < 0:
            res = []
            eof = False
            while not eof:
                try:
                    res.append(pickle.load(f))
                except EOFError:
                    eof = True
        else: #number==0
            res = pickle.load(f)
    finally:
        if close:
            f.close()

    return res

def fpickle(object,fileorname,usecPickle=True,protocol=None,append=False):
    """
    Pickle an object to a specified file.

    :param object: the python object to pickle
    :param fileorname: The file from which to unpickle objects
    :type fileorname: a file name string or a :class:`file` object
    :param usecPickle:
        If True, the :mod:`cPickle` module is to be used in place of
        :mod:`pickle` (cPickle is faster).
    :type usecPickle: bool
    :param protocol:
        Pickle protocol to use - see the :mod:`pickle` module for details on
        these options. If None, the most recent protocol will be used.
    :type protocol: int or None
    :param append:
        If True, the object is appended to the end of the file, otherwise the
        file will be overwritten (if a file object is given instead of a
        file name, this has no effect).
    :type append: bool

    """

    if usecPickle:
        import cPickle as pickle
    else:
        import pickle

    if protocol is None:
        protocol = pickle.HIGHEST_PROTOCOL

    if isinstance(fileorname,basestring):
        f = open(fileorname,'a' if append else 'w')
        close = True
    else:
        f = fileorname
        close = False

    try:
        pickle.dump(object,f,protocol=protocol)
    finally:
        if close:
            f.close()



class _PyfitsOpener(object):
    """
    A context manager for use with :func:`open_with_pyfits` (see that docstring
    for details).
    """
    def __init__(self, *args,**kwargs):
        import pyfits
        self.f = pyfits.open(*args,**kwargs)
    def __enter__(self):
        return self.f
    def __exit__(self, *exc_info):
        self.f.close()

def open_with_pyfits(*args,**kwargs):
    """
    This is a convinience function for :mod:`pyfits`, allowing the following
    usage in python 2.5 or above (in 2.5, ``from __future__ import
    with_statement`` is needed at the beginning of the file)::

        from astropysics.utils.io import open_with_pyfits

        with open_with_pyfits('filename') as f:
            h = f[0].header
            d = f[0].data

            # ... do something more with the file...

        # at this indent level, the pyfits file is now closed. It will also be
        # closed if an exception is thrown.

    The arguments are the same as for :func:`pyfits.open`.

    """
    return _PyfitsOpener(*args,**kwargs)

#<-----------------------Data retrieval and caching---------------------------->

def get_package_data(dataname):
    """
    Use this function to load data files distributed with the astropysics
    source code.

    :param str dataname:
        The name of a file in the package data directory.
    :returns: The content of the requested file as a string.

    .. note::
        The data accessed by this function is distinct from the data accessed
        via :func:`get_data`. Package data is crucial basic data included in the
        astropysics source distribution, while standard data is for larger or
        optional data files that are downloaded as needed.

    """
    from .. import __name__ as rootname
    from .. import __file__ as rootfile
    from pkgutil import get_loader
    from os.path import dirname
    path = dirname(rootfile)+'/data/'+dataname
    return get_loader(rootname).get_data(path)

def _readrem(remote,reportprogress=False):
    """
    Reads the provided remote url and returns the result, possible reporting
    on progress.

    Progress will be approximately once every percent, or less often for files
    under 100 kb. If the file size is unknown, messages progress will slowly
    become less frequent.

    :param remote: A :class:`urllib.addinfourl` object with the remote url.
    :param reportprogress:
        See :func:`set_data_download_reporter`

    :returns: The read data of the remote, e.g. remote.read()
    """
    from math import ceil

    if reportprogress:
        if reportprogress is True:
            def reportprogress(progress,started,ended,url):
                from sys import stdout

                if started is False:
                    stdout.write('Downloading file from '+url)
                    stdout.write('\nDownload progress:  0%')
                    stdout.flush()
                elif ended is True:
                    stdout.write('\n')
                    stdout.flush()
                elif type(progress) is float:
                    perc = progress*100.
                    stdout.write('\rDownload progress:%3i%%'%perc)
                    stdout.flush()
                else: #type must be int
                    kb = progress/1000.
                    mb = kb/1000. if kb > 1000 else None
                    gb = mb/1000. if mb > 1000 else None
                    if gb is not None:
                        stdout.write('\rDownload progress:%.1f GB'%gb)
                    elif gb is not None:
                        stdout.write('\rDownload progress:%.1f MB'%mb)
                    else:
                        stdout.write('\rDownload progress:%.1f kB'%kb)
                    stdout.flush()

        if not callable(reportprogress):
            raise TypeError('Invalid reportprogress argument %s'%reportprogress)

        if 'content-length' in remote.info():
            totbytes = int(remote.info()['content-length'])
            if totbytes > 102400:
                nbytesper = int(ceil(totbytes/100.))
            else:
                nbytesper = 1024
        else:
            totbytes = 1
            nbytesper = 1024

        res = []
        currres = True
        donebytes = 0
        reportprogress(0,False,False,remote.url)
        while currres!='':
            currres = remote.read(nbytesper)
            res.append(currres)
            donebytes += len(currres)
            fracorbytes = donebytes/totbytes
            reportprogress(fracorbytes,True,False,remote.url)
        reportprogress(fracorbytes,True,True,remote.url)
        return ''.join(res)

    else:
        return remote.read()

#TODO: Document these in future config docs
_data_store = _io_config.get('data_store',True)
_data_reporter = _io_config.get('data_reporter',True)

def get_data(dataurl,asfile=False,localfn=None):
    """
    Retrieves a data file from a remote source (usually the internet), and
    optionally caches that data locally. See :func:`set_data_store` for control
    of caching behavior, and :func:`set_data_download_reporter` for control of
    download progress reporting.

    :param str dataurl:
        The URL of the data to be retrieved. If not a URL (i.e. does not have a
        ':' somewhere), this will be interpreted as a local file name inside the
        astrpysics data directory.
    :param bool asfile:
        If True, a file-like object is returned that can be used to access the
        data. Otherwise, a string with the full content of the file is returned.
    :param localfn:
        The filename to use for saving (or loading, if the file is present)
        locally. If it is None, the filename will be inferred from the URL
        if possible. This file name is always relative to the astropysics data
        directory (see :func:`astropysics.config.get_data_dir`). This has no
        effect if :func:`set_data_store` is set to False.

    :returns: A file-like object or a string (see `asfile`)

    :raises urllib2.URLError:
        If the `dataurl` request is a URL and cannot be found.
    :raises IOError:
        If the dataurl is requested as a local data file and not found.

    .. note::
        The data accessed by this function is distinct from the data accessed
        via :func:`get_package_data`. Package data is crucial basic data
        included in the astropysics source distribution, while standard data is
        for larger or optional data files that are downloaded as needed.

    """
    import os,urlparse,urllib2,shelve,contextlib
    from ..config import get_data_dir

    global _data_store,_data_reporter

    #if not a URL, just look locally
    if urllib2.splittype(dataurl)[0] is None:
        fn = os.path.join(get_data_dir(),dataurl)
        if asfile:
            return open(fn)
        else:
            with open(fn) as f:
                return f.read()

    store = _data_store
    reportprogress = _data_reporter

    if store:
        datadir = get_data_dir()
        with contextlib.closing(shelve.open(os.path.join(datadir,'urlmap'))) as url2fn:
            #decide if the url needs to be retrieved
            if dataurl not in url2fn or store=='refresh' or \
              not os.path.exists(os.path.join(datadir,url2fn[dataurl])):
                with contextlib.closing(urllib2.urlopen(dataurl)) as remote:
                    rinfo = remote.info()
                    urldata = None
                    if localfn is None:
                        #figure out the local name from the data/URL
                        if 'Content-Disposition' in rinfo:
                            #often URLs that redirect to a download provide the fielname in the header info
                            localfn = rinfo['Content-Disposition'].split('filename=')[1]
                        else:
                            #otherwise fallback on the url filename
                            localfn = urlparse.urlsplit(dataurl)[2].split('/')[-1]
                            if localfn.strip()=='':
                                #url does not have a path, so try to use the title of the page
                                localfn = None
                                if 'html' in rinfo['content-type']:
                                    urldata = _readrem(remote,reportprogress)
                                    starttag = urldata.find('<title>')
                                    endtag = urldata.find('<title>')
                                    if startag>-1 and endtag>-1:
                                        localfn = urldata[starttag:endtag]
                        if localfn is None:
                            raise urllib2.URLError('Could not determine local file name for the URL '+dataurl)
                    if urldata is None:
                        urldata = _readrem(remote,reportprogress)
                #now save the downloaded data
                with open(os.path.join(datadir,localfn),'w') as f:
                    f.write(urldata)
                url2fn[dataurl] = localfn

            handle = open(os.path.join(datadir,url2fn[dataurl]),'r')
    else:
        handle = urllib2.urlopen(dataurl)

    if asfile:
        return handle
    else:
        return handle.read()

def set_data_store(store=_data_store):
    """
    Sets the behavior of the caching mechanism when :func:`get_data` is called.
    The value can subsequently be retrieved via :func:`get_data_store`.

    :param store:
        If True, when :func:`get_data` is called the data file will only be
        downloaded if that URL is accessed for the first time. If False, the
        data will always be retrieved from the remote source and not saved. If
        it is the string 'refresh', the file is downloaded regardless of whether
        or not it is present, but the downloaded version will be used in future
        calls where it is True.

    :raises TypeError: If `store` is an inappropriate type

    """
    if isinstance(store,basestring):
        if store!='refresh':
            raise TypeError('invalid data store string '+store)
    elif store is not True and store is not False:
        raise TypeError('data store must be True, False, or a string')

    global _data_store
    _data_store = store

def get_data_store():
    """
    Returns the current behavior of the caching mechanism for :func:`get_data`.
    See :func:`set_data_store` for details of the possible values

    :returns:
        The current value of the data store. Possible values are listed in
        :func:`set_data_store`.

    """
    return _data_store

def set_data_download_reporter(reportprogress=_data_reporter):
    """
    Sets the behavior of the download progress reporter for :func:`get_data`.
    The value can subsequently be retrieved via
    :func:`get_data_download_reporter`.

    :param reportprogress:
        If False, no progress reports go out while file is downloading. If True,
        a progress message is printed at the command line. Otherwise, this must
        be a callable.

    .. note::
        If `reportprogress` is a callable, it will be called as
        func(progress,started,finished,url), where `progress` is either a float
        between 0 and 1 indicating the fraction of the file downloaded (if the
        file size is known), or the number of bytes downloaded if the file size
        is unknown (1 or larger). `started` is True if the download has begun,
        and `finished` is True when the download has finished. The function will
        always be called at the beginning as func(0,False,False) and once at the
        end as func(?,True,True). `url` is a string with the target URL.

    :raises TypeError: If `reportprogress` is an inappropriate type

    .. note::
        If :func:`set_data_store` has been set to False, no download progress
        will be reported, as the data will be retrieved immediately without
        download to a local file.

    """
    if not (callable(reportprogress) or \
            reportprogress is True or reportprogress is False):
        raise TypeError('data store must be True, False, or a callable')

    global _data_reporter
    _data_reporter = reportprogress

def get_data_download_reporter():
    """
    Returns the current behavior of the download progress reporter for
    :func:`get_data`. See :func:`set_data_download_reporter` for details of the
    possible values.

    :returns:
        The current value of the data store. Possible values are listed in
        :func:`set_data_download_reporter`.

    """
    return _data_reporter


#<-----------------------General IO utilities---------------------------------->

def loadtxt_text_fields(fn,fieldline=1,asrecarray=True,updatedict=None,**kwargs):
    """
    Load a text file into a structured numpy array where the field names and
    types are inferred from a line (typically the first) in the file.

    It must begin with a comment (specified by the 'comments' keyword) and have
    fields (sperated by the `delimiter`, default whitespace) composed of a
    colon-seperated field name and valid numpy dtype (e.g. 'f' for floats or 'i'
    for integers).

    If a third colon-seperated component is present, it is a python expression
    (with no spaces) that can be used to derive the value for that column given
    values from the other columns.  numpy functions can be used with the prefix
    'np'.

    An example field line might be::

        #data1:f derived:f:data1+data2**2 data2:i

    This will result in a 3-field record array with fields 'data1', 'derived',
    and 'data2', where data1 and data2 are the columns of the input file.

    :param fn: The file name to load.
    :type fn: string
    :param fieldline: The line number that has the field information.
    :type fieldline: int
    :param bool asrecarray:
        If True, returns a :class:`numpy.recarray`. Otherwise, returns a regular
        :class:`numpy.ndarray` with a structured dtype.
    :param updatedict:
        A dictionary that will be updated with the loaded data such that the
        keys are the field names and the values are arrays with the data. If
        None, this will be ignored.
    :type updatedict: dict or None

    extra keywords are passed into :func:`numpy.loadtxt`

    :returns: A numpy record array or regular array with data from the text file.

    .. note::
        The `updatedict` parameter is seful for injecting the loaded file into
        the local namespace: ``loadtxt_text_fields(...,updatedict=locals())``.

    """

    fieldline -= 1 #1-based to 0-based
    typedelim = ':'

    if 'dtype' in kwargs:
        raise ValueError("can't use field lines")

    comments = kwargs['comments'] if 'comments' in kwargs else '#'
    delimiter = kwargs['delimiter'] if 'delimiter' in kwargs else None

    with open(fn,'r') as f:
        for i,l in enumerate(f):
            if i >= fieldline:
                l = l.replace(comments,'')
                fields = l.split(delimiter)
                break
    try:
        dtype = []
        loadtype = []
        exprd = {}
        for fi in fields:
            t = tuple(fi.split(typedelim))
            if len(t)==2:
                dtype.append(t)
                loadtype.append(t)
            elif len(t)==3:
                dtype.append(t[:2])
                exprd[t[0]] = t[2]
            else:
                raise ValueError('Invalid field line %s'%fi)

        loadtype = np.dtype(loadtype)
        dtype = np.dtype(dtype)

    except TypeError: #figure out where the problem was
        for fi in fields:
            t = tuple(fi.split(typedelim))
            try:
                np.dtype(t[1])
            except TypeError:
                raise TypeError('dtype code {1} invalid for field {0}'.format(*t))
            #if something else went wrong, re-raise
            raise

    loadarr = np.loadtxt(fn,dtype=loadtype,**kwargs)
    arr = np.empty(loadarr.size,dtype=dtype)

    ldict = {}
    for fin in loadarr.dtype.names:
        arr[fin] = ldict[fin] = loadarr[fin]

    for fin,expr in exprd.iteritems():
        arr[fin] = eval(expr,globals(),ldict)

    if updatedict is not None:
        d = dict([(fi,arr[fi]) for fi in arr.dtype.names])
        updatedict.update(d)

    if asrecarray:
        return arr.view(np.recarray)
    else:
        return arr


class FixedColumnDataParser(object):
    """
    Parses a data file composed of lines with a fixed set of columns
    with the same number of bytes in each

    """


    def __init__(self,skiprows=0,firstcolindx=1,commentchars='#'):
        """
        :param skiprows: The number of rows to skip initially
        :type skiprows: int
        :param indexing:
            The value of the first column as used in the definition of the upper
            and lower edges of fields - that is, if 1, the "1-indexing"
            convntion is used, or if 0, standard python "0-indexing" is used.
        :type oneindexed: positive int
        :param commentchars:
            The comment character - any lines that start with these characters
            will be ignored.
        :type commentchars: string (interpreted as sequence of characters)
        """
        self.cols = {}
        self.skiprows = skiprows
        self.commentchars = list(commentchars)
        self.firstcolindx = firstcolindx

    def _overlapcheck(self,lower,upper,exclude=None):
        for n,(l,u,f,conv) in self.cols.iteritems():
            if n != exclude and lower <= u and upper >= l:
                raise ValueError('input range %i-%i overlaps on column %s'%(lower,upper,n))



    def addColumn(self,name,lower,upper,format=None,converters=None):
        """
        :param name: Name of the column
        :type name: string
        :param lower:
            lowest character index for this column *including* this one.
            Convention for the value of first column is controlled by
            :attr:`firstcol` attribute of the object - default is 1-indexing
        :type lower: int
        :param upper:
            highest character index for this column *including* this one.
            Convention for the value of first column is controlled by
            :attr:`firstcol` attribute of the object - default is 1-indexing
        :type upper: int
        :param format:
            Format for this column, or None to infer from the contents.
        :type format: input to :class:`numpy.dtype` or None
        :param converters:
            A dictionary mapping one string to another - all entries in the
            parsed data file with a key will be converted to the corresponding
            value. If None, no conversion will take place.
        :type converters: string dict or None

        """
        if format is not None: #check that its a valid format specifier
            format = np.dtype(format)
        self._overlapcheck(lower,upper,name)
        self.cols[name] = (lower,upper,format,converters)

    def addColumnsFromFile(self,fn,columnlinestart=None,sep=None,useskiprows=True,
                                maxcols=100):
        """
        Adds columns parsed from a file (typically a data file that will
        afterwards be read).

        The parsed file is expected to have lines that begin with the
        `columnlinestart` argument and the rest of the line should be able to be
        split (using `sep`) into three or four columns in order
        name,lower,upper,format .


        :param fn: File name of the file to parse
        :type fn: string
        :param columnlinestart:
            String that indicates the line is a column line, or None to just use
            the first `maxcols` columns.
        :type columnlinestart: string or None
        :param sep:
            The seperator string to split the line or None for whitespace.
        :type sep: string or None
        :param useskiprows:
            If True, the :attr:`skiprows` attribute of this object will be used
            to determine how many rows to skip.
        :type useskiprows: bool
        :param maxcols:
            Maximum number of columns to parse before quitting (or None to be
            unlimited).
        :type maxcols: int or None

        :returns: number of columns added

        **Examples**

        For a file that starts like::

            ID 0 3 int
            data 4 10 float
            moredata 11 15 int

        use addColumnsFromFile('filename',columnlinestart=None,sep=None,maxcols=3,firstcolindex=0)

        And for a file with format specifier::

            #... ID,1,4,int
            #... data,5,11,float
            #... moredata,12,16,int

        use addColumnsFromFile('filename',columnlinestart='#... ',sep=',')


        """
        #used below to properly re-adjust string dtypes
        def maybe_convert_string(linesplit):
            if linesplit[3] in ('a','S','string','str'):
                chrs = linesplit[2]-linesplit[1]+1
                return 'S'+str(chrs)
            else:
                return linesplit[3]

        with open(fn) as f:
            if useskiprows:
                for i in range(self.skiprows):
                    f.readline()
            icols = 0
            for l in f:
                ls = l.strip()
                if columnlinestart is None:
                    linesplit = ls.split(sep) if sep else ls.split()
                    linesplit[1] = int(linesplit[1])
                    linesplit[2] = int(linesplit[2])
                    if len(linesplit)>=4:
                        linesplit[3] = maybe_convert_string(linesplit)
                    if len(linesplit)>=5:
                        linesplit[4] = sep.join(linesplit[4:])
                        linesplit = linesplit[:5]
                        if isinstance(linesplit[4],basestring):
                            linesplit[4] = eval(linesplit[4])
                            if not isinstance(linesplit[4],dict):
                                raise TypeError('5th column is not a dict string')
                    self.addColumn(*linesplit)
                    icols += 1
                elif ls.startswith(columnlinestart):
                    lsr = ls.replace(columnlinestart,'')
                    linesplit = (lsr.split(sep) if sep else lsr.split())
                    linesplit[1] = int(linesplit[1])
                    linesplit[2] = int(linesplit[2])
                    if len(linesplit)>=4:
                        linesplit[3] = maybe_convert_string(linesplit)
                    if len(linesplit)>=5:
                        linesplit[4] = sep.join(linesplit[4:])
                        linesplit = linesplit[:5]
                        if isinstance(linesplit[4],basestring):
                            linesplit[4] = eval(linesplit[4])
                            if not isinstance(linesplit[4],dict):
                                raise TypeError('5th column is not a dict string')
                    self.addColumn(*linesplit)
                    icols += 1
                if maxcols is not None and icols >= maxcols:
                    break
        return icols

    def delColumn(self,name):
        """
        Removes a column by name

        :param name: name of the column to remove
        :type name: string
        """
        del self.cols[name]

    def parseFile(self,fn,maskedarray=False):
        """
        Parse a file that follows this object's format.

        :param fn: File name of the file to parse
        :type fn: string
        :params maskedarray:
            If True, the function returns a masked array, otherwise a record
            array.
        :type maskedarray: bool

        :returns:
            If `maskedarray` is False, a tuple (recarr,masks) where recarr is a
            :class:`numpy.core.records.recarray` and masks is a dictionary
            mapping column names to masks. The masks are True if the value is
            valid, and False if not.  If `maskedarray` is True, a
            :class:`numpy.ma.core.MaskedArray` is returned.
        """
        lists = dict([(k,list()) for k in self.cols])
        addi = -int(self.firstcolindx)

        if len(self.cols) == 0:
            raise IndexError('No columns defined in FixedColumnDataParser')

        with open(fn) as f:
            for i,row in enumerate(f):
                rs = row.strip()
                validrow = i >= self.skiprows and rs != '' and rs[0] not in self.commentchars
                if validrow:
                    for n,(l,u,f,convs) in self.cols.iteritems():
                        if convs is None:
                            li,ui = l+addi,u+addi+1
                            lists[n].append(row[li:ui])
                        else:
                            li,ui = l+addi,u+addi+1
                            val = row[li:ui]
                            lists[n].append(convs.get(val,val))

        sorti = np.argsort([l for l,u,f,convs in self.cols.values()])
        sortnames = np.array(self.cols.keys())[sorti]

        alist = []
        masks = []
        for n in sortnames:
            lst = [e if e.strip()!='' else None for e in lists[n]]
            f = self.cols[n][2]
            masks.append((n,np.array([l is not None and l is not '' for l in lst])))

            if f is not None:
                print 'converting field',n,'to type',f
                if f.kind == 'i':
                    arr = np.array([0 if l is None else l for l in lst])
                else:
                    arr = np.array(lst)
                alist.append(arr.astype(f))
            else:
                arr = np.array(lst)
                alist.append(arr)

        recarr = np.rec.fromarrays(alist,names=','.join(sortnames))
        masks = dict(masks)

        if maskedarray:
            marr = np.ma.MaskedArray(recarr)
            for n,m in masks.iteritems():
                marr[n].mask = ~m
            return marr
        else:
            return recarr,masks

    def writeFile(self,fn,data,masks=None,colstart='# '):
        """
        Writes a data array out to a data file using this format.

        :param fn: The file to write to.
        :type fn: string
        :param data:
            The data array to write. Must match this objects column names and
            formats.
        :type data:
            Structured :class:`numpy.array` or
            :class:`numpy.core.records.recarray`.
        :param masks:
            If supplied, maps field names to boolean arrays where True indicates
            the data value *should* be written. If None, all data is written.
        :type masks: None or dict
        :param colstart: the string to use to indicate a column specifier
        :type colstart: string

        :except TypeError: If the data dtype doesn't match the columns.
        """
        for n,(l,u,f) in self.cols.iteritems():
            for k,(dt,c) in data.dtype.fields.iteritems():
                if k==n:
                    if dt != f:
                        raise TypeError('Data dtype %s does not match column type %s'%(dt,f))
                    break
            else:
                raise TypeError('Data field %s not a column'%n)

        sortedcols = sorted([(l,n) for n,(l,u,f) in self.cols.iteritems()])
        sortedcols = [e[1] for e in sortedcols]


        with open(fn,'w') as f:
            for n in sortedcols:
                l,u,f = self.cols[n]
                f.write(colstart+' '.join((n,l,u,f))+'\n')
            f.write('\n')

            for rec in data:
                oldu=0
                for n in sortedcols:
                    l,u,f = self.cols[n]
                    spaces = l-oldu-1
                    chrs = u-l+1
                    if spaces>0:
                        f.write(' '*spaces)
                    f.write(str(rec[n])[:(u-l+1)])
                f.write('\n')



@_add_docs_and_sig(FixedColumnDataParser.parseFile,FixedColumnDataParser.addColumnsFromFile)
@_add_docs(FixedColumnDataParser)
def loadtxt_fixed_column_fields(fn,fncol=None,skiprows=0,comments='#',
                                columnlinestart='#:',columnsep=' ',maxcols=None,
                                oneindexed=True,maskedarray=False,
                                updatedict=None):
    """
    Loads a fixed column data file using :class:`FixedColumnDataParser`.

    :param fncol:
        The file from which to load the column data. If None, it will assumed to
        be the same as `fn`.

    :param updatedict:
        A dictionary that will be updated with the loaded data such that the
        keys are the field names and the values are arrays with the data. If
        None, it will be ignored.

    For other arguments see extra documentation below and
    :class:`FixedColumnDataParser`.

    Additional keywords are passed into :func:`numpy.loadtxt`.

    :returns:
        :class:`~numpy.recarray` or a regular :class:`~numpy.ndarray` with data
        loaded from the text file.

    To match the default format, data files to load should look like::

        #:col1 1 4 int
        #:col2 6 6 bool
        #:col3 8 10 S3
        #:col4 12 16 float {'------':''}
        # could have a comment here if you want
        1239 1 abc 12.325
        2489 0 zyx ------
        9526 1 qwe 89632.


    :meth:`FixedColumnDataParser.parseFile` describes the return types and
    `maskedarray` parameter.
    {docstr:parseFile}

    For details on arguments controlling how columns are inferred, See
    :meth:`FixedColumnDataParser.addColumnsFromFile`.
    {docstr:addColumnsFromFile}

    """
    if fncol is None:
        fncol = fn
    fcp = FixedColumnDataParser(skiprows,oneindexed,comments)
    fcp.addColumnsFromFile(fncol,columnlinestart,columnsep,True,maxcols)

    arr = fcp.parseFile(fn,maskedarray=maskedarray)
    if updatedict is not None:
        d = dict([(fi,arr[fi]) for fi in arr.dtype.names])
        updatedict.update(d)
    return arr

def load_tipsy_format(fn):
    """
    This function loads a file in the Tipsy ASCII format
    (http://www-hpcc.astro.washington.edu/tipsy/man/readascii.html)
    and outputs a dictionary with entries for grouped data
    """
    dout = {}
    with open(fn) as f:
        ntotal,ngas,nstar = [int(e) for e in f.readline().split()]
        ndark = ntotal-nstar-ngas
        dims = int(f.readline().strip())
        time = float(f.readline().strip())
    dout['ntotal'] = ntotal
    dout['ndark'] = ndark
    dout['nstar'] = nstar
    dout['ngas'] = ngas
    dout['dims'] = dims
    dout['time'] = time

    A = np.loadtxt(fn,skiprows=3)
    i = 0
    j = ntotal
    dout['mass'] = A[i:j]
    dout['mass_gas'] = dout['mass'][:ngas]
    dout['mass_dark'] = dout['mass'][ngas:ndark]
    dout['mass_star'] = dout['mass'][(ngas+ndark):]

    pos = []
    for d in range(dims):
        i = j
        j += ntotal
        pos.append(A[i:j])

    dout['pos'] = np.array(pos,copy=False)
    dout['pos_gas'] = dout['pos'][:,:ngas]
    dout['pos_dark'] = dout['pos'][:,ngas:ndark]
    dout['pos_star'] = dout['pos'][:,(ngas+ndark):]

    vel = []
    for d in range(dims):
        i = j
        j += ntotal
        vel.append(A[i:j])
    dout['vel'] = np.array(vel,copy=False)
    dout['vel_gas'] = dout['vel'][:,:ngas]
    dout['vel_dark'] = dout['vel'][:,ngas:ndark]
    dout['vel_star'] = dout['vel'][:,(ngas+ndark):]

    if ndark>0:
        i = j
        j += ndark
        dout['softening_dark'] = A[i:j]

    if nstar>0:
        i = j
        j += nstar
        dout['softening_star'] = A[i:j]

    if ngas>0:
        i = j
        j += ngas
        dout['density'] = A[i:j]
        i = j
        j += ngas
        dout['temperature'] = A[i:j]
        i = j
        j += ngas
        dout['sph_smoothing'] = A[i:j]
        i = j
        j += ngas
        dout['metals_gas'] = A[i:j]

    if nstar>0:
        i = j
        j += nstar
        dout['metals_star'] = A[i:j]
        i = j
        j += nstar
        dout['formation_time'] = A[i:j]

    i=j
    j+= ntotal
    dout['potential_en'] = A[i:j]
    dout['potential_en_gas'] = dout['potential_en'][:ngas]
    dout['potential_en_dark'] = dout['potential_en'][ngas:ndark]
    dout['potential_en_star'] = dout['potential_en'][(ngas+ndark):]

    if j!=A.size:
        from warnings import warn
        warn('Tipsy file expected to have %i entries, but found %i!'%(j,A.size))

    return dout

#<------------------------VOTable related-------------------------------------->
#add this back in if anything actually uses vo.table in the future
#try:
#    import vo.table
#except ImportError:
#    from warnings import warn
#    warn('vo.table not found - VOTable processing limited to VOTableReader class')

class VOTableReader(object):
    """
    This class represents a VOTable. Currently, it is read-only, and will
    probably not be enhanced due to the existence of Michael Droettboom's
    vo.table package (http://trac6.assembla.com/astrolib).
    """

    dtypemap = {"boolean":'b',
                "bit":None,
                "unsignedByte":'u1',
                "short":'i2',
                "int":'i4',
                "long":'i8',
                "char":'a',
                "unicodeChar":'U',
                "float":'f4',
                "double":'f8',
                "floatComplex":'c8',
                "doubleComplex":'c16'} #maps VOTable data types to numpy

    def __init__(self,s,filename=True):
        """
        instantiate a VOTable object from an XML VOTable

        If filename is True, the input string will be interpreted as
        a filename for a VOTable, otherwise s will be interpreted as
        an XML-formatted string with the VOTable data
        """
        from xml.dom import pulldom
        if filename:
            events = pulldom.parse(s)
        else:
            events = pulldom.parseString(s)

        self._tables = {} #keys are table names
        self._masks = {} #keys are table names
        self._resnames = {} #keys are table names
        self._fields = {} #keys are table names
        self._pars = {} #keys are table names
        self._parexts = {}
        self._tabdescs = {} #keys are table names
        self._resdescs = {} #keys are resource names from resnames
        self._res = []
        self.description=''

        namesep = ':'

        rcounter=0
        tcounter=0
        voparams=[]
        inres=None
        intab=None
        indat=None


        for ev,n in events:
            if ev == 'START_ELEMENT':
                nm = n.tagName
                if nm == 'RESOURCE':
                    if inres:
                        raise NotImplementedError('no support for nested resources')

                    rcounter+=1
                    resparams = []

                    if n.attributes.has_key('name'):
                        inres = n.attributes['name'].value
                    else:
                        inres = 'res%i'%rcounter

                elif nm == 'TABLE':
                    if not inres:
                        raise RuntimeError('table outside of resource - invalid VOTable?')

                    tcounter+=1
                    tabparams = []
                    fields = []

                    if intab:
                        raise RuntimeError('nested tables - invalid VOTable?')

                    if n.attributes.has_key('ref'):
                        raise NotImplementedError('table refs not yet implemented')

                    if n.attributes.has_key('name'):
                        intab = inres+namesep+n.attributes['name'].value
                    else:
                        intab = 'tab%i'%tcounter
                    if intab in self._tables:
                        intab = intab+'_'+str(tcounter)

                    if n.attributes.has_key('nrows'):
                        nrows = n.attributes['nrows'].value
                    else:
                        nrows = None

                elif nm == 'DATA':
                    if not intab:
                        raise RuntimeError('Data not in a table - invalid VOTable?')
                    params = []
                    params.extend(voparams)
                    params.extend(resparams)
                    params.extend(tabparams)

                    indat = True

                #data types
                elif nm == 'TABLEDATA':
                    events.expandNode(n)
                    array,mask = self._processTabledata(n,fields,params,nrows)

                elif nm == 'BINARY':
                    raise NotImplementedError('Binary data not implemented')

                elif nm == 'FITS':
                    raise NotImplementedError('FITS data not implemented')

                elif nm == 'PARAM':
                    events.expandNode(n)
                    if inres and not intab:
                        resparams.append(self._extractParam(n))
                    elif intab:
                        tabparams.append(self._extractParam(n))
                    else:
                        voparams.append(self._extractParam(n))
                elif nm == 'FIELD':
                    if not intab:
                        raise RuntimeError('nested tables - invalid VOTable?')
                    events.expandNode(n)
                    fields.append(self._extractField(n))
                elif nm == 'GROUP':
                    raise NotImplementedError('Groups not implemented')
                elif nm == 'DESCRIPTION':
                    events.expandNode(n)
                    n.normalize()
                    desc = n.firstChild.nodeValue
                    if inres:
                        if intab:
                            self._tabdescs[intab] = desc
                        else:
                            self._resdescs[inres] = desc
                    else:
                        self.description = desc


            elif ev == 'END_ELEMENT':
                nm = n.tagName
                if nm == 'RESOURCE':
                    inres = None
                elif nm == 'TABLE':

                    self._resnames[intab] = inres
                    self._fields[intab] = fields
                    self._applyParams(intab,params)
                    self._tables[intab] = array
                    self._masks[intab] = mask
                    intab = None
                    del array,mask,params,fields #do this to insure nothing funny happens in parsing - perhaps remove later?
                elif nm == 'DATA':
                    indat = False
            if ev == 'CHARACTERS':
                pass

    def _applyParams(self,intab,params):
        self._pars[intab] = dict([(str(p[0]),p[1]) for p in params])
        self._parexts[intab] = dict([(str(p[0]),p[2]) for p in params])

    def getTableNames(self):
        return self._tables.keys()

    def getTableResource(self,table=0):
        nm = self._tableToName(table)
        return self._resnames[nm]

    def getTableParams(self,table=0):
        nm = self._tableToName(table)
        fs = self._pars[nm]

    def getTableParamExtras(self,table=0):
        nm = self._tableToName(table)
        fs = self._parexts[nm]

    def getTableFieldNames(self,table=0):
        nm = self._tableToName(table)
        fs = self._fields[nm]
        return [str(f[0]) for f in fs]

    def getTableFieldDtypes(self,table=0):
        nm = self._tableToName(table)
        fs = self._fields[nm]
        return dict([(str(f[0]),str(f[1])) for f in fs])

    def getTableFieldExtras(self,table=0):
        nm = self._tableToName(table)
        fs = self._fields[nm]
        return dict([(str(f[0]),f[2]) for f in fs])

    def getTableArray(self,table=0):
        nm = self._tableToName(table)
        return self._tables[nm]

    def getTableMask(self,table=0):
        nm = self._tableToName(table)
        return self._masks[nm]

    def _tableToName(self,table):
        if isinstance(table,basestring):
            if table not in self._tables:
                raise KeyError('table %s not found'%table)
            return table
        else:
            i = int(table)
            return self._tables.keys()[i]

    def _extractField(self,n):
        n.normalize()
        name = n.attributes['name'].value
        desc = n.getElementsByTagName('DESCRIPTION')
        if len(desc) == 0:
            desc = ''
        elif len(desc) == 1:
            desc = desc[0].firstChild.nodeValue
        else:
            raise RuntimeError('multiple DESCRIPTIONs found in field %s - invalid VOTable?'%name)

        dtype = self.dtypemap[n.attributes['datatype'].value]
        if n.attributes.has_key('arraysize'):
            szs = n.attributes['arraysize'].value.strip()
            if dtype == 'a' or dtype == 'U':
                if 'x' in szs:
                    raise NotImplementedError('multidimensional strings not yet supported')
                elif szs == '*':
                    raise NotImplementedError('unlimited length strings not yet supported')
                dtype = dtype+szs.replace('*','') #fixed length strings
            else:
                raise NotImplementedError('array primitives not yet supported')

        #val = n.getElementsByTagName('VALUES') #not using this

        extrad={'DESCRIPTION':desc}
        keys = n.attributes.keys()
        for k in ('name','arraysize','datatype'):
            if k in keys:
                keys.remove(k)
        for k in keys:
            extrad[k] = n.attributes[k].value
        if len(extrad)==0:
            extrad = None

        return name,dtype,extrad

    def _extractParam(self,n):
        n.normalize()
        name = n.attributes['name'].value
        desc = n.getElementsByTagName('DESCRIPTION')
        if len(desc) == 0:
            desc = ''
        elif len(desc) == 1:
            desc = desc[0].firstChild.nodeValue
        else:
            raise RuntimeError('multiple DESCRIPTIONs found in param %s - invalid VOTable?'%name)
        dtype = self.dtypemap[n.attributes['datatype'].value]
        val = n.attributes['value'].value
        npval = np.array(val,dtype=dtype)
        #val = n.getElementsByTagName('VALUES') #not using this

        extrad={'DESCRIPTION':desc}
        keys = n.attributes.keys()
        for k in ('name','arraysize','datatype'):
            if k in keys:
                keys.remove(k)
        for k in keys:
            extrad[k] = n.attributes[k].value
        if len(extrad)==0:
            extrad = None

        return name,val,extrad

    def _processTabledata(self,n,fields,params,nrows):
        n.normalize()
        dt = np.dtype([(str(f[0]),str(f[1]))for f in fields])
        if not nrows:
            nrows = 0
            for c in n.childNodes:
                if c.nodeName == 'TR':
                    nrows+=1

        arr = np.ndarray(nrows,dtype=dt)
        mask = np.ones((nrows,len(dt)),dtype=bool)

        i = 0
        for c in n.childNodes:
            if c.nodeName == 'TR':
                j = 0
                for d in c.childNodes:
                    if d.nodeName == 'TD':
                        if d.hasChildNodes():
                            arr[i][j] = d.firstChild.nodeValue
                        else:
                            arr[i][j] = 0
                            mask[i][j] = False
                    elif c.nodeType == 3 and c.nodeValue.strip()=='':
                        pass #ignore whitespace text
                    else:
                        raise RuntimeError('non-TD inside TR - invalid VOTable?')
                    j+=1
                i+=1
            elif c.nodeType == 3 and c.nodeValue.strip()=='':
                pass #ignore whitespace text
            else:
                raise RuntimeError('non-TR inside Tabledata - invalid VOTable?')

        return arr,mask

class VOTable(VOTableReader): #name for backwards-compatibility
    def __init__(*args,**kwargs):
        from warnings import warn
        warn('VOTable class name deprecated - use VOTableReader',DeprecationWarning)
        return VOTableReader.__init__(*args,**kwargs)


#<--------------------------Spectrum loaders----------------------------------->

def load_wcs_spectrum(fn,fluxext=None,errext=None,hdrext=None,errtype='err'):
    """
    Loads a linear spectrum from a fits file with WCS keywords 'CD1_1' or
    'CDELT1' and 'CRVAL_1'.

    :param fn: File name of fits file with spectrum to load.
    :type fn: string
    :param fluxext:
        Fits extension to use for the spectrum or None to use the first that has
        the WCS form (i.e. 'CD1_1'/'CDELT1' and 'CRVAL_1' present).
    :type fluxext: int or None
    :param errext: Fits extension to use for the errors or None to skip errors.
    :type errext: int of None
    :param hdrext:
        The fits extension to use to look for the header keywords - if None,
        this will be the same as `fluxext`, or ignored if `fluxext` is None.
    :type hdrext: int or None
    :param errtype: Form of error data if present: 'err','ierr','var', or 'ivar'
    :type errtype: str

    :returns: An :class:`astropysics.spec.Spectrum` object

    :except ValueError: If the fits file is missing necessary keywords
    :except IOError: If `fluxext` is None and no HDUs have the correct form

    """
    from ..spec import Spectrum

    with open_with_pyfits(fn) as f:
        if fluxext is None:
            for i,hdu in enumerate(f):
                if 'CRVAL1' in hdu.header and ('CD1_1' in hdu.header or 'CDELT1' in hdu.header):
                    fluxext = hdrext = i
                    break
            else:
                raise IOError('Could not locate an HDU with CRVAL1 and CD1_1/CDELT1')

        hdr = f[fluxext if hdrext is None else hdrext].header

        if 'CTYPE1' not in hdr:
            from warnings import warn
            warn('No CTYPE1 keyword, so uncertain if this is a linear spectrum')
        else:
            if not hdr['CTYPE1'] == 'LINEAR':
                raise ValueError('Spectrum coordinates must be linear')

        if 'CRVAL1' not in hdr:
            raise ValueError('missing header keyword CRVAL1')

        if 'CD1_1' in hdr:
            dispersion = hdr['CD1_1']
        elif 'CDELT1' in hdr:
            dispersion = hdr['CDELT1']
        else:
            raise ValueError('missing header keyword CD1_1 or CDELT1')

        flux = f[fluxext].data
        #sometimes this is a 1xN array instead of just N:
        if len(flux.shape)==2 and flux.shape[0]==1:
            flux = flux[0]
        x = np.arange(flux.size)*dispersion+hdr['CRVAL1']

        if errext is not None:
            err = f[errext].data
        else:
            err = None

        if errtype == 'err':
            fobj = Spectrum(x,flux,err=err)
        elif errtype == 'ierr':
            fobj = Spectrum(x,flux,err=1/err)
        elif errtype == 'var':
            fobj = Spectrum(x,flux,ivar=1/err)
        elif errtype == 'ivar':
            fobj = Spectrum(x,flux,ivar=err)
        else:
            raise ValueError('Unrecognized errtype %s'%errtype)

        fobj.name = fn
        return fobj


def load_deimos_spectrum(fn,plot=False,extraction='horne',retdata=False,smoothing=None):
    """
    extraction type can 'horne' or 'boxcar'

    if smoothing is positive, it is gaussian sigmas, if negative,
    boxcar pixels

    returns Spectrum object with ivar, [bdata,rdata]
    """
    import pyfits
    from ..spec import Spectrum

    if 'spec1d' not in fn or 'fits' not in fn:
        raise ValueError('loaded file must be a 1d spectrum from DEEP/spec2d pipeline')

    if hasattr(pyfits,'NP_pyfits'):
        try:
            from pyfits.NP_pyfits import _extensionNameCaseSensitive as pfcasesen
        except ImportError:
            from pyfits.NP_pyfits import EXTENSION_NAME_CASE_SENSITIVE as pfcasesen
    else:
        try:
            from pyfits.core import _extensionNameCaseSensitive as pfcasesen
        except ImportError:
            from pyfits.core import EXTENSION_NAME_CASE_SENSITIVE as pfcasesen
    if extraction == 'horne':
        extname = 'Horne' if pfcasesen else 'HORNE'
    elif extraction == 'boxcar':
        extname = 'Bxspf'if pfcasesen else 'BXSPF'
    else:
        raise ValueError('unrecgnized extraction type %s'%extraction)

    with open_with_pyfits(fn) as f:
        extd = dict([(f[i].header['EXTNAME'],i) for i in range(1,len(f))])
        if extname+'-B' in extd and extname+'-R' in extd:
            bi,ri = extd[extname+'-B'],extd[extname+'-R']
            bd,rd=f[bi].data,f[ri].data
            x=np.concatenate((bd.LAMBDA[0],rd.LAMBDA[0]))
            flux=np.concatenate((bd.SPEC[0],rd.SPEC[0]))
            ivar=np.concatenate((bd.IVAR[0],rd.IVAR[0]))
            sky=np.concatenate((bd.SKYSPEC[0],rd.SKYSPEC[0]))

            changei = len(bd.LAMBDA[0])
        elif extname+'-B' in extd or extname+'-R' in extd:
            if extname+'-B' in extd:
                d = f[extname+'-B'].data
            else:
                d = f[extname+'-R'].data

            x = d.LAMBDA[0]
            flux = d.SPEC[0]
            ivar = d.IVAR[0]
            sky = d.SKYSPEC[0]

            changei = len(d.LAMBDA[0])-1
        else:
            raise KeyError('neither side found for extraction type '+extname )

        fobj = Spectrum(x,flux,ivar=ivar)
        fobj.sky = sky

        if smoothing:
            if smoothing < 0:
                fobj.smooth(smoothing*-1,filtertype='boxcar')
            else:
                fobj.smooth(smoothing,replace=True)

        if plot:
            from matplotlib import pyplot as plt
            if plot != 'noclf':
                plt.clf()
            plt.plot(fobj.x[:changei],fobj.flux[:changei],'-b')
            plt.plot(fobj.x[changei:],fobj.flux[changei:],'-r')

        fobj.name = fn
        if retdata:
            return fobj,bd,rd
        else:
            return fobj

def load_all_deimos_spectra(dir='.',pattern='spec1d*',extraction='horne',
                            smoothing=None,verbose=True):
    """
    loads all deimos spectra found in the specified directory that
    match the requested pattern and returns a list of the file names
    and the Spectrum objects.

    extraction and smoothing are the same as for load_deimos_spectrum

    verbose indicates if information should be printed

    returns dictionary mapping file names to Spectrum objects
    """
    from glob import glob
    from os.path import join

    fns = glob(join(dir,pattern))
    fns.sort()
    specs = []

    fnrem=[]
    for i,fn in enumerate(fns):
        if verbose:
            print 'Loading spectrum',fn
        try:
            s = load_deimos_spectrum(fn,False,extraction,False,smoothing)
            specs.append(s)
        except Exception,e:
            if verbose:
                print 'Exception loading spectrum',fn,'skipping...'
            fnrem.append(i)
    for i in reversed(fnrem):
        del fns[i]

    return dict(zip(fns,specs))

def _load__old_spylot_spectrum(s,bandi):
    from ..spec import Spectrum
    x=s.getCurrentXAxis()
    f=s.getCurrentData(bandi=bandi)
    if s.isContinuumSubtracted():
        e=s.getContinuumError()
    else:
        e=s.getWindowedRMSError(bandi=bandi)
    return Spectrum(x,f,e)


#clean up namespace
del get_config
