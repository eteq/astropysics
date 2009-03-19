#Copyright (c) 2009 Erik Tollerud (etolleru@uci.edu) 

"""
This module contains functions and classes for various forms of data/file
input and output As well as library retrieval for built-in data.
"""

from __future__ import division
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

#<------------------------VOTable classes-------------------------------------->
class VOTable(object):
    """
    This class represents a VOTable.  Currently, it is read-only, although 
    later re-designs may change this
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
        return [str(f[1]) for f in fs]
    
    def getTableFieldExtras(self,table=0):
        nm = self._tableToName(table)
        fs = self._fields[nm]
        return [f[2] for f in fs]
    
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
        
        extrad={}
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
        
        extrad={}
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
            
        
#class VOTableElem(object):
#    pass

#class VOTable(VOTableElem):
#    description=''
#    infos=[]
#    groups=[]
#    params=[]
#    resources=[]
    
#class Resource(VOTableElem):
#    description=''
#    infos=[]
#    groups=[]
#    params=[]
#    resources=[]
#    links=[]
#    tables=[]
    
#class Table(VOTableElem):
#    description=''
#    infos=[]
#    groups=[]
#    params=[]
#    links=[]
#    fields=[]
#    data=None
    
#class Field(VOTableElem):
#    description=''
#    values=None
#    links=[]
    
#class Data(VOTableElem):
#    pass

#class Group(VOTableElem):
#    description=''
#    fieldrefs=[]
#    params=[]
#    paramrefs=[]
#    groups=[]
    
#class Info(VOTableElem):
#    description=''
#    values=None
#    links=[]

#class Values(VOTableElem):
#    min=None
#    max=None
#    options=[]

#<--------------------------Spectrum loaders----------------------------------->

def load_deimos_spectrum(fn,plot=True,extraction='horne',retdata=False,smoothing=None):
    """
    extraction type can 'horne' or 'boxcar'
    
    returns Spectrum object with ivar, [bdata,rdata]
    """
    import pyfits
    from .spec import Spectrum
    if 'spec1d' not in fn or 'fits' not in fn:
        raise ValueError('loaded file must be a 1d spectrum from DEEP/spec2d pipeline')
    
    if extraction == 'horne':
        extname = 'HORNE'
    elif extraction == 'boxcar':
        extname = 'BXSPF'
    else:
        raise ValueError('unrecgnized extraction type %s'%extraction)
    
    f=pyfits.open(fn)
    try:
        extd = dict([(f[i].header['EXTNAME'],i) for i in range(1,len(f))])
        bi,ri = extd[extname+'-B'],extd[extname+'-R']
        bd,rd=f[bi].data,f[ri].data
        x=np.concatenate((bd.LAMBDA[0],rd.LAMBDA[0]))
        flux=np.concatenate((bd.SPEC[0],rd.SPEC[0]))
        ivar=np.concatenate((bd.IVAR[0],rd.IVAR[0]))
        sky=np.concatenate((bd.SKYSPEC[0],rd.SKYSPEC[0]))
        
        changei = len(bd.LAMBDA[0])
        
        fobj = Spectrum(x,flux,ivar)
        fobj.sky = sky
        
        if smoothing:
            fobj.smooth(smoothing,replace=True)
        
        if plot:
#            from operator import isMappingType
#            if not isMappingType(plot):
#                plot={}
            from matplotlib import pyplot as plt
            plt.plot(fobj.x[:changei],fobj.flux[:changei],'-b')
            plt.plot(fobj.x[changei:],fobj.flux[changei:],'-r')
        
        if retdata:
            return fobj,bd,rd
        else:
            return fobj
    finally:
        f.close()
        
def load_deimos_templates(fns):
    """
    This will generate a dictionary of Spectrum objects from the specified  
    templates 
    fn can be a single file, a list of files, or a unix-style pattern
    """
    import pyfits
    from .spec import Spectrum
    from operator import isSequenceType
    from warnings import warn
    
    if isinstance(fns,basestring):
        from glob import glob
        fns = glob(fns)
        
    if not isSequenceType(fns):
        raise ValueError('improper filename format')
    
    tempd={}
    xd={}
    for fn in fns:
        f = pyfits.open(fn)
        try:
            h=f[0].header
            d=f[0].data
            if len(d.shape) == 1:
                d = (d,)
            if len(d) == 1:
                if h.has_key('NAME0'):
                    k=h['NAME0'].strip()
                elif h.has_key('NAME'):
                    k=h['NAME'].strip()
                elif h.has_key('OBJECT'):
                    k=h['OBJECT'].strip()
                else:
                    k = fn.strip()
                    warn('could not infer spectrum name - using filename %s as key name'%k)
                    
                xd[k] = 10**(h['COEFF0']+h['COEFF1']*np.arange(d[0].shape[-1]))
                tempd[k] = d[0].copy()
            else:
                x = 10**(h['COEFF0']+h['COEFF1']*np.arange(d.shape[-1]))
                for i,spec in enumerate(d):
                    if  h.has_key('NAME%i'%i):
                        k=h['NAME%i'%i].strip()
                    elif h.has_key('OBJECT') and h.has_key('EIGEN%i'%i):
                        k = '%s%i'%(h['OBJECT'].strip(),h['EIGEN%i'%i])
                    elif h.has_key('NAME'):
                        k = '%s-%i'%(h['NAME'].strip(),i)
                    elif h.has_key('OBJECT'):
                        k = '%s-%i'%(h['OBJECT'].strip(),i)
                    else:
                        k = '%s-%i'%(fn.strip(),i)
                        warn('could not infer spectrum name - using filename %s as key name'%k)
                        
                    xd[k] = x
                    tempd[k]=spec.copy() #TODO: see if its better to not copy here
                    
        finally:
            f.close()
            
    return dict(((k,Spectrum(xd[k],tempd[k])) for k in tempd))
    
def load_spylot_spectrum(s,bandi):
    from .spec import Spectrum
    x=s.getCurrentXAxis()
    f=s.getCurrentData(bandi=bandi)
    if s.isContinuumSubtracted():
        e=s.getContinuumError()
    else:
        e=s.getWindowedRMSError(bandi=bandi)
    return Spectrum(x,f,e)
