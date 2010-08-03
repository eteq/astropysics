#Copyright 2010 Erik Tollerud
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

============================================================
publication -- tools for preparing astronomy papers in LaTeX 
============================================================

The :mod:`publication` module contains classes and funtions to assist in
preparing papers or other publications with an emphasis on common astronomy
journals. The tools here are for use with `LaTeX
<http://www.latex-project.org/>`_ as the actual authoring tools.

.. todo:: examples

Classes and Inheritance Structure
---------------------------------

.. inheritance-diagram:: astropysics.publication
   :parts: 1

Module API
----------

"""

from __future__ import division,with_statement

#<-------------------------------OLD--------------------->

class TeXDoc(object):
    def __init__(self,filename):
        self.lines=ls=[]
        with open(filename) as f:
            self.lines=[l for l in f]
            
    def save(self,fn):
        with open(fn,'w') as f:
            for l in self.lines:
                f.write(l)
                
    def stripComments(self):
        comms=self._stripComms(self.lines)
        return comms,len(comms)
        
    def _stripComms(self,ls):
        comms,itorem=[],[]
        for i,l in enumerate(ls):
            if l.startswith('%'):
                itorem.append(i)
                comms.append(l)
            else:
                lastc=''
                for j,c in enumerate(l):
                    if c=='%' and lastc!='\\':
                        self.lines[i]=l[:j]
                        comms.append(l[j:])
                        break
                    lastc=c
        itorem.sort(reverse=True)
        for i in itorem:
            del ls[i]
                
        return comms
    
    def getDocClass(self):
        """
        returns ind,cls,sty
        """
        ind,sty,cls=None,None,None
        import re
        rex=re.compile(r'.*\\documentclass(\[(.*)\])?{(.*)}.*')
        for i,l in enumerate(self.lines):
            if not l.lstrip().startswith('%'):
                m=rex.match(l)
                if m:
                    if ind is not None:
                        raise ValueError(r'multiple \documentclass entries found')
                    ind = i
                    sty = m.group(2)
                    cls = m.group(3)
        return ind,cls,sty
    
    def setDocClass(self,newclass,newstyle=None,commentold=True):
        i=self.getDocClass()[0]
        posttag=('[%s]'%newstyle if newstyle else '')+'{%s}'%newclass
        if commentold:
            self.lines[i]='%'+self.lines[i]
            self.lines.insert(i,r'\documentclass'+posttag+'\n')
        else:
            self.lines[i]=r'\documentclass'+posttag+'\n'
        
    def printDoc(self):
        for l in self.lines:
            print l
            
    def generateFigureMapping(self,applymapping=False):
        from re import compile
        d={}
        regex=compile(r'.*(\includegraphics|\plot(one|two|fiddle)){(.+?)}.*')
        
        figletters=['a','b','c','d','e','f','g','h','i','j','k','l','m']
        infig=0
        figcount=0
        fignames=[]
        templines=self.lines[:]
        self._stripComms(templines)
        for i,l in enumerate(templines):
            flatl=l.replace(' ','')
            if r'\begin{figure' in flatl:
                infig=1
            elif r'\end{figure' in flatl:
                infig=-1
            
            if infig == 1:
                m=regex.match(l)
                while m:
                    fin=m.group(3)
                    fignames.append(fin.strip())
                    l=l.replace('{'+fin+'}','')
                    m=regex.match(l)
            elif infig == -1:
                figcount+=1
                if len(fignames) > 1:
                    for j,n in enumerate(fignames):
                        d[n]='f%i%s'%(figcount,figletters[j])
                else:
                    d[fignames[0]]='f%i'%figcount
                fignames=[]
                infig=0
                
        for k,v in d.items():
            if '.' in k:
                newk=k.split('.')[0]
                del d[k]
                d[newk]=v
        if applymapping:
            self.applyFigureMapping(d)
        return d
    
    def applyFigureMapping(self,mapping):
        from re import compile
        regex=compile(r'.*(\includegraphics|\plot(one|two|fiddle))({.+})?.*')
        for k,v in mapping.iteritems():
            for i,l in enumerate(self.lines):
                if k in l:
                    m=regex.match(l)
                    if m and ('{%s}'%k in m.group(3) or '{%s.'%k in m.group(3)):
                        self.lines[i]=l.replace(k,v)
                        
    def findTag(self,tagname,repl=None):
        """
        looks for a line of the form \tagname[{r1}{r2}{r3}]
        if repl is not None, the curly braces will be replaced (either as a sequence or string)
        """
        if type(repl) is str:
            repl=(repl,)
        
        from re import compile
        regex=compile(r'.*\\'+tagname+r'({.*})*.*')
        d={}
        if repl is None:
            for i,l in enumerate(self.lines):
                m = regex.match(l)
                
                if m is not None:
                    ti = l.index('\\'+tagname)
                    if '%' not in l[:ti]:
                        lg=m.group(1)
                        if lg is not None:
                            d[i]=m.group(1)
        else:
            replstr='{%s}'%('}{'.join(repl))
            for i,l in enumerate(self.lines):
                m = regex.match(l)
                if m is not None:
                    lg=m.group(1)
                    if lg is not None:
                        lrep=lg.split('}{')
                        if len(lrep) != len(repl):
                            raise ValueError('Tag found, but wrong number of elements in line %i :\n%s'%(i,l))
                        d[i]=lg
                    
            for i,s1 in d.iteritems():
                self.lines[i]=self.lines[i].replace(s1,replstr)
            
        return d
                    
                
    
def prep_for_apj_pub(fnbase=None,targetdir='./pubApJ/',replbib=False):
    from os import sep,mkdir,chmod
    from os.path import exists,isdir
    from shutil import copy,rmtree
    from glob import glob
    import tarfile
    
    if fnbase is None:
        from glob import glob
        texfiles=glob('*.tex')
        if len(texfiles)!=1:
            raise IOError('No unique tex file in current directory!')
        fnbase=texfiles[0].replace('.tex','')
        
    if exists(targetdir) and isdir(targetdir):
        ri=raw_input('Target directory '+targetdir+' exists - overwrite([y]/n)?')
        if 'n' in ri:
            return
        rmtree(targetdir)
        
    elif exists(targetdir):
        raise IOError('targetdir is not a directory')
    mkdir(targetdir)
        
    if not targetdir.endswith(sep):
        targetdir = targetdir+sep
    
    texfn=fnbase+'.tex'
    
    td=TeXDoc(texfn)
    print 'stripped',td.stripComments()[1],'comments'
    fm=td.generateFigureMapping(True)
    td.setDocClass('aastex','manuscript')
    td.save(targetdir+'ms.tex')
    
    for k,v in fm.iteritems():
        copy(k+'.eps',targetdir+v+'.eps')
    
    
    if replbib:
        raise NotImplementedError
    else:
        #copy over bibliography files and rename
        bblfn=fnbase+'.bbl'
        
        if not exists(texfn):
            raise IOError('.tex file not found')
        
        if exists(bblfn):
            copy(bblfn,targetdir+'ms.bbl')
        else:
            raise IOError('no .bbl file found')
        #replace \bibliography{} entry with ms name
        td.findTag('bibliography',repl=('ms',''))
        
    f=tarfile.open(targetdir+'ms.tar.gz','w|gz')
    try:
        fnl=glob(targetdir+'*')
        for fn in fnl:
            f.add(fn)
    finally:
        f.close()
    chmod(targetdir+'ms.tar.gz',33188)
        
def prep_for_arxiv_pub(fnbase=None,targetdir='./pubArXiv/',fnoutbase='ms',replbib=False):
    from os import sep,mkdir,chmod,system
    from os.path import exists,isdir
    from shutil import copy,rmtree
    from glob import glob
    import tarfile
    from warnings import warn
    from subprocess import Popen,PIPE
    
    if fnbase is None:
        from glob import glob
        texfiles=glob('*.tex')
        if len(texfiles)!=1:
            raise IOError('No unique tex file in current directory!')
        fnbase=texfiles[0].replace('.tex','')
    if exists(targetdir) and isdir(targetdir):
        ri=raw_input('Target directory '+targetdir+' exists - overwrite([y]/n)?')
        if 'n' in ri:
            return
        rmtree(targetdir)
        
    elif exists(targetdir):
        raise IOError('%s is not a directory'%targetdir)
    mkdir(targetdir)
    if not targetdir.endswith(sep):
        targetdir = targetdir+sep
    
    texfn=fnbase+'.tex'
    
    td=TeXDoc(texfn)
    print 'stripped',td.stripComments()[1],'comments'
    fm=td.generateFigureMapping(True)
    td.setDocClass('emulateapj')
    
    for k,v in fm.iteritems():
        copy(k+'.eps',targetdir+v+'.eps')
    
    
    if replbib:
        raise NotImplementedError
    else:
        #copy over bibliography files and rename
        bblfn=fnbase+'.bbl'
        
        if not exists(texfn):
            raise IOError('.tex file not found')
        
        if exists(bblfn):
            copy(bblfn,targetdir+fnoutbase+'.bbl')
        else:
            raise IOError('no .bbl file found')
        
        print 'td ',td.findTag('bibliography',repl=(fnoutbase,''))
        
    td.save(targetdir+fnoutbase+'.tex')
    
    f=tarfile.open(targetdir+fnoutbase+'.tar.gz','w|gz')
    try:
        fnl=glob(targetdir+'*')
        for fn in fnl:
            f.add(fn)
    finally:
        f.close()
    chmod(targetdir+fnoutbase+'.tar.gz',33188)
    
    retcode1 = system('cd '+targetdir+';latex '+fnoutbase+'.tex')
    retcode2 = system('cd '+targetdir+';latex '+fnoutbase+'.tex')
    retcode3 = system('cd '+targetdir+';latex '+fnoutbase+'.tex')
    if retcode1 != 0 or retcode2 != 0 or retcode3 != 0:
        warn('LaTeX compilation failed!')
    else:
        print 'converting dvi to pdf for',texfn
        system('cd '+targetdir+';dvipdf '+fnoutbase+'.dvi')
        
        ti = td.findTag('title').values()[0][1:-1]
        au = td.findTag('author').values()
        if len(au)>1:
            warn('multiple authors found!')
        else:
            fignums = []
            for fcode in fm.values():
                fcode = fcode[1:]
                
                done = False
                while not done:
                    try:
                        if fcode != '':
                            fignums.append(int(fcode))
                        done = True
                    except ValueError:
                        fcode = fcode[:-1]
            
            proc = Popen('pdfinfo '+targetdir+fnoutbase+'.pdf',stdout=PIPE,shell=True)
            pout = proc.communicate()[0]
            proc.wait()
            
            print '\nTitle:\n',ti
            au = au[0][1:-1] #strip enclosing curly braces
            au = [a.split('\\')[0] for a in au.split(',') if '{' not in a.split('\\')[0] and '}' not in a.split('\\')[0] ]
            print '\nAuthors:\n',', '.join(au)
            
            for l in pout.split('\n'):
                if 'Pages' in l:
                    print '\nPages: ',l.replace('Pages:','').strip(),'Figures:',len(set(fignums))
                    break
            else:
                warn('no page entry found in pdfinfo')
            

    
#<-----------------------NEW BELOW HERE------------------------------------>

class TeXNode(object):
    """
    TODO:DOC
    """
    
class TeXDoc(TeXNode):
    """
    TODO:DOC
    """
    
#issue warning if abstract too long
#strip comments
#.tar.gz up with appropriate name and date            
def prep_for_arxiv_pub(texfile):
    raise NotImplementedError

#Tasks: redo figures into f##[l].eps and move the files
#set class to aastex
#strip comments
#fix any deluxetable* -> deluxetable (and rotate)
#issue warning if abstract too long
#.tar.gz up with appropriate name and date
def prep_for_apj_pub(texfile,authorlast):
    raise NotImplementedError