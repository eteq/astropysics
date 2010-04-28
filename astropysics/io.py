#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""

==
io
==

The ``io`` module contains classes and functions for loading and saving data in
various relevant formats used in astronomy, as well as convinience functions
retrieval for built-in data.

.. todo:: examples


Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.io
   :parts: 1

Module API
----------

"""

from __future__ import division,with_statement
import numpy as np

try:
    import pyfits
except ImportError:
    from warnings import warn
    warn('pyfits not found - all FITS-related IO will not work')
    
    
#<-----------------------Internal to package----------------------------------->
    
def _get_package_data(dataname):
    """
    Use this function to load data files distributed with astropysics in the 
    astropysics/data directory
    
    dataname is the file name of a file in the data directory, and a string
    with the contents of the file will be returned
    """
    from . import __name__ as rootname
    from . import __file__ as rootfile
    from pkgutil import get_loader
    from os.path import dirname
    path = dirname(rootfile)+'/data/'+dataname
    return get_loader(rootname).get_data(path)

#<-----------------------General IO utilities---------------------------------->

def loadtxt_text_fields(fn,fieldline=0,typedelim=':',asrecarray=True,**kwargs):
    """
    this uses numpy.loadtxt to load a text file into a numpy record 
    array where the field names are inferred from a commented text line.
    
    the format for the field line is:
    
    #field1:typecode1[:unit1] field2:typecode2[:unit2] (... )
    
    with the character in place of the : optionally selected with the 
    typedelim keyword.  Any comments in the line will be removed. 
    type codes are the same as those used in numpy.dtype.  if units are provided
    they are a numerical factor to multiply the column by - if ommitted, it is
    assumed to be 1
    
    fieldline tells which line of the file to use to find the line with
    field information
    
    asrecarray converts the array to a record array before returning, allowing
    attribute-style access
    
    kwargs are passed into numpy.loadtxt
    """
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
        factors = {}
        for fi in fields:
            t = tuple(fi.split(typedelim))
            dtype.append(t[:2])
            if len(t)>2:
                factors[t[0]] = np.array(t[2],dtype=t[1])
        dtype = np.dtype(dtype)
    except TypeError: #figure out where the problem was
        for fi in fields:
            t = tuple(fi.split(typedelim))
            try:
                dtype(t[1])
            except TypeError:
                raise TypeError('dtype code {1} invalid for field {0}'.format(t))
            #if something else went wrong, re-raise
            raise
    arr = np.loadtxt(fn,dtype=dtype,**kwargs)
    for n,f in factors.iteritems():
        arr[n]*=f
    
    
    
    if asrecarray:
        return arr.view(np.recarray)
    return arr


class FixedColumnData(object):
    """
    Parses a data file composed of lines with a fixed set of columns
    with the same number of bytes in each
    
    """
    
    
    def __init__(self,skiprows=0,oneindexed=True,commentchars='#'):
        """
        :param skiprows: The number of rows to skip initially
        :type skiprows: int
        :param oneindexed: 
            If True, the column numbers will be treated such that the first
            column is column #"1" - otherwise, the first column is #"0".
        :type oneindexed: bool
        :param commentchars: 
            The comment character - any lines that start with these characters
            will be ignored.
        :type commentchars: string (interpreted as sequence of characters)
        """ 
        self.cols = {}
        self.skiprows = skiprows
        self.commentchars = list(commentchars)
        self.oneindexed = oneindexed
        
    def _overlapcheck(self,lower,upper,exclude=None):
        for n,(l,u,f) in self.cols.iteritems():
            if n != exclude and lower <= u and upper >= l:
                raise ValueError('input range %i-%i overlaps on column %s'%(lower,upper,n))
            
        
    
    def addColumn(self,name,lower,upper,format=None):
        """
        :param name: Name of the column
        :type name: string
        :param lower: 
            lowest character index for this column *including* this one.
        :type lower: int
        :param upper:
            highest character index for this column *including* this one.
        :type upper: int
        :param format: 
            Format for this column, or None to infer from the contents.
        :type format: input to :class:`numpy.dtype` or None
        
        """
        if format is not None: #check that its a valid format specifier
            format = np.dtype(format)
        self._overlapcheck(lower,upper,name)
        self.cols[name] = (lower,upper,format)
        
    def addColumnsFromFile(self,fn,linestart=None,sep=None,useskiprows=True,
                                maxcols=None):
        """
        Adds columns parsed from a data file itself.
        
        The parsed file is expected to have lines that begins with the
        `linestart` argument and the rest should be able to be split into three
        or four columns in order name,lower,upper,format .
        
        :param fn: File name of the file to parse
        :type fn: string
        :param linestart: String that indicates the line is a column line.
        :type linestart: 
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
                if linestart is None:
                    linesplit = ls.split(sep) if sep else ls.split()
                    linesplit[1] = int(linesplit[1])
                    linesplit[2] = int(linesplit[2])
                    if len(linesplit)>=4:
                        linesplit[3] = maybe_convert_string(linesplit)
                    self.addColumn(*linesplit)
                    icols += 1
                elif ls.startswith(linestart):
                    lsr = ls.replace(linestart,'')
                    linesplit = (lsr.split(sep) if sep else lsr.split())
                    linesplit[1] = int(linesplit[1])
                    linesplit[2] = int(linesplit[2])
                    if len(linesplit)>=4:
                        linesplit[3] = maybe_convert_string(linesplit)
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
        
    def parseFile(self,fn):
        """
        Parse a file that follows this object's format.
        
        :param fn: File name of the file to parse
        :type fn: string
        
        :returns: 
            A tuple (recarr,masks) where recarr is the data and masks is a 
            dictionary mapping column names to masks.  The masks are True if
            the value is valid, and False if not.
        :rtype: a :class:`numpy.core.records.recarray` and a dict
        """
        lists = dict([(k,list()) for k in self.cols])
        addi = -int(self.oneindexed)
        
        with open(fn) as f:
            for i,row in enumerate(f):
                rs = row.strip()
                validrow = i >= self.skiprows and rs != '' and rs[0] not in self.commentchars
                if validrow:
                    for n,(l,u,f) in self.cols.iteritems():
                        li,ui = l+addi,u+addi+1
                        lists[n].append(row[li:ui])
        
        sorti = np.argsort([l for l,u,f in self.cols.values()])
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
                alist.append(np.array(lst))
        
        recarr = np.rec.fromarrays(alist,names=','.join(sortnames))
        
        return recarr,dict(masks)
    
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
                    
                    
                
         
def loadtxt_fixed_column_fields(fn,fncol=None,skiprows=0,comments='#',columnlinestart='#',
                                columnsep=None,maxcols=None,oneindexed=True):
    """
    Loads a fixed column data file using :class:`FixedColumnData`. 
    
    If `fncol` is None, it will assumed to be the same as `fn`.
    
    For details on arguments controlling how columns are inferred, See
    :meth:`FixedColumnData.addColumnsFromFile`. For other arguments see
    :class:`FixedColumnData`.
    
    See :meth:`FixedColumnData.parseFile` for return type
    """
    if fncol is None:
        fncol = fn
    fcp = FixedColumnData(skiprows,oneindexed,comments)
    fcp.addColumnsFromFile(fncol,columnlinestart,columnsep,True,maxcols)
    return fcp.parseFile(fn)

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
try:
    import vo.table
except ImportError:
    from warnings import warn
    warn('vo.table not found - VOTable processing limited to VOTableReader class')

class VOTableReader(object):
    """
    This class represents a VOTable.  Currently, it is read-only, and will 
    probably not be enhanced due to the existence of Michael Droettboom's vo
    package
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

def load_wcs_spectrum(fn,fluxext=1,errext=None,hdrext=None,errtype='err'):
    """
    Loads a spectrum from a fits file with WCS keywords CD1_1 and CRVAL_1
    
    fluxext specifies the extension to use for the flux data, while errext
    specifies err source (or None for no errors) - errtype gives the form
    of the error data - either 'err','ierr','var', or 'ivar'
    
    hdrext specifies which extension to use to look for the header keywords - 
    by default this is the same as the flux extension
    """
    import pyfits
    from .spec import Spectrum
    
    f=pyfits.open(fn)
    try:
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
    finally:
        f.close()
        

def load_deimos_spectrum(fn,plot=False,extraction='horne',retdata=False,smoothing=None):
    """
    extraction type can 'horne' or 'boxcar'
    
    if smoothing is positive, it is gaussian sigmas, if negative, 
    boxcar pixels
    
    returns Spectrum object with ivar, [bdata,rdata]
    """
    import pyfits
    from .spec import Spectrum
    if 'spec1d' not in fn or 'fits' not in fn:
        raise ValueError('loaded file must be a 1d spectrum from DEEP/spec2d pipeline')
    
    if extraction == 'horne':
        extname = 'Horne' if pyfits.NP_pyfits._extensionNameCaseSensitive else 'HORNE'
    elif extraction == 'boxcar':
        extname = 'Bxspf'if pyfits.NP_pyfits._extensionNameCaseSensitive else 'BXSPF'
    else:
        raise ValueError('unrecgnized extraction type %s'%extraction)
    
    f=pyfits.open(fn)
    try:
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
    finally:
        f.close()
        
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
    from .spec import Spectrum
    x=s.getCurrentXAxis()
    f=s.getCurrentData(bandi=bandi)
    if s.isContinuumSubtracted():
        e=s.getContinuumError()
    else:
        e=s.getWindowedRMSError(bandi=bandi)
    return Spectrum(x,f,e)
