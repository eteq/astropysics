#!/usr/bin/env python
"""
<author> <date>

<project description>

<DELETE COMMENTS BELOW HERE AFTER READING>
This file is a template for a data analysis project (typically to go with a 
paper or two).  The main usage scenario is to plass functions and classes 
before the "if __name__=='__main__' line.  Plots/figures go there as well,
with the @plotfunc decorator.  Then the command line usage can be used to 
re-generate any/all plots, and interactive work in ipythons supported by doing:

    run <filename.py>
    make_plots('plotname',...)

Or if changes are made to some of the plots, the following can be used:

    import <filename>
    reload(<filename>)
    <filename>.make_plots('plotname',locals(),...)

Additionally, this script can be called at the command line to show or save
figures.  This is especially useful in conjunction with the mainfig argument 
of plotfunc to save all figures intended for use in a paper.


**Example**

Replace 'ADD CLASSES/FUNCTIONS HERE' with:

@plotfunc(mainfig='paper1')
def aplot(x,data1,data2):
    subplot(2,1,1)
    plot(x,data1,label='1')
    plot(x,data2,label='2')
    
    xlabel('My x')
    ylabel('My data')
    legend(loc=0)
    
    subplot(2,1,2)
    scatter(data1,data2)
    xlabel('data 1')
    ylabel('data 2')
    xlim(2,4.5)
    ylim(1,5)
    
and replace 'ADD ON-RUN OPERATIONS HERE' with:

    x = linspace(0,1,100)
    data1 = x*2.5 + 2
    data2 = data1 + randn(100)
    
    
then you could do:

./myscrit.py -m paper1

at the command line, and a directory 'paper1' will be created with files 
'aplot.eps' and 'aplot.pdf'.

"""

from __future__ import division,with_statement

from numpy import *
from numpy.random import *
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt

#modify this list to change the default file formats for saving figures
defaultsavetypes = ['eps','pdf'] 

#<--------------------Plot support functions----------------------------------->
_plotreg = {}
_mainfigs = {True:[]}
    
def plotfunc(figtype='mpl',figsize=None,tweakbbox=None,mainfig=False):
    """
    :param figtype:
        Specifies the type/format of the figure. Can be any of:
        
            * 'mpl'
                Matplotlib figure saving to both eps and pdf.
            * 'mpl<ext1>,<ext2>'
                Matplotlib figure saving each figure to filenames that have the
                requested extension.
            * 'mpltoeps'
                saves as png and then converts png->eps
            * 'mayavi' 
                Mayavi figure, defaults to png
            * 'mayavi<ext1>,<ext2>' 
                Mayavi figure, saving to the specified extension types
            
    :param figsize: 
        Specifies the size of the figure as a (width,height) tuple or None for default.
    :param tweakbbox:
        Specifies the bounding box to be used for eps output as a (x1,y1,x2,y2)
        tuple or None for the default bounding box.
    :param bool mainfig:
        If True, this is a "main" figure - e.g. one that will actually appear in
        a paper and is used by the '-m' command line option. Can also be a
        string, in which case the figures will be saved in a matching name
        directory and can be selected seperately from those that are True.
    """
    def deco(f):
        fname = f.func_name
        if 'mpl' in figtype:
            if figsize is None:
                fsz = (8,8)
            else:
                fsz = figsize
        elif 'mayavi' in figtype:
            if figsize is None:
                fsz = (800,600)
            else:
                fsz = figsize
        else:
            raise ValueError('unrecognized plot type %s'%figtype)
        _plotreg[fname]=(f,fsz,figtype,tweakbbox)
        if mainfig:
            if mainfig not in _mainfigs:
                _mainfigs[mainfig] = []
            _mainfigs[mainfig].append(fname)
        return f
    
    if callable(figtype):
        f=figtype
        figtype='mpl'
        return deco(f)
    else:
        return deco
    
def plot_names():
    """
    Returns the names of all registered plots for this file.  
    Also returns the 'main' plots, as a dictiontionary mapping the main types to
    figure names
    
    :returns: 2-tuple (plotnames,mainnames)
    
    """
    return _plotreg.keys(),_mainfigs.copy()

#To save by default in a different directory, change the "save='.'" line at the
#start of this function to "save='/path/to/dir'"
def make_plots(plots,argd=None,figs=None,save=False,overwrite=True,
               showfigs=True,savetypes=defaultsavetypes):
    """
    Generate plots and maybe save them.
    
    :param plots: A comma-seperated string of plots to show
    :param argd:
        The dictionary to use to find the arguments of the plot function. If
        None, global variables will be used.
    :param figs: 
        A list of figure numbers to use for each of the figures in `plots`. If
        None, the current figure will be used as the starting point.
    :param str save:
        A directory in which to save the plots. If False, they will not be
        saved, if True, the current directory will be used.
    :param bool overwrite:
        If True and plots are being saved, any existing plot will be
        overwritten. Otherwise, the figure will not be saved.
    :parm bool showfigs:
        If True, the figures will be shown (in whatever GUI is the correct one
        for the type).
    :param list savetype: 
        A list of strings specifying the extensions of the file formats for
        saving figures that don't explicitly specify a file format. Change the
        `defaultsavetypes` variable at the top of the script to change the
        default (initially, ['eps','pdf'])
    
    """
    import os,shutil,inspect,subprocess   
    
    if save is True:
        save='.'
    if isinstance(save,basestring):
        if save.endswith(os.sep):
            save = os.sep[:-len(os.sep)]
        if not os.path.exists(save):
            os.mkdir(save)
        save+=os.sep  
    
    if plots is None:
        plots = _plotreg.keys()
    elif isinstance(plots,basestring):
        if '*' in plots:
            import re
            plots = plots.replace('*','.*')
            regex = re.compile(plots)
            plots = [p for p in plot_names()[0] if regex.match(p)]
        else:
            plots = plots.split(',')
        
    if figs is None:
        figs = arange(len(plots))+1
         
    if len(figs) != len(plots):
        raise ValueError("plots don't match figs")
    
    if argd is None:
        argd = globals()
            
    plotd={}
    for p,fnum in zip(plots,figs):
        if p not in _plotreg:
            raise KeyError('plotname %s not found'%p)
        plotd[p] = (_plotreg[p][0],fnum,_plotreg[p][1],_plotreg[p][2],_plotreg[p][3])
        
    isinter = plt.isinteractive()
    plt.ioff()
    anympl = False
    try:
        for fname,(figf,fignum,figsize,ftype,tweakbbox) in plotd.iteritems(): 
            print 'figure',fname
            
            if 'mpl' in ftype:
                anympl = True
                if isinter:
                    plt.ion()
                plt.close(fignum)
                plt.figure(fignum,figsize=figsize)
                plt.ioff()
            elif 'mayavi' in ftype:
                from enthought.mayavi import mlab as M
                if len(M.get_engine().scenes) > 0:
                    M.gcf().parent.close_scene(M.gcf())
                f=M.figure(fignum,size=figsize) 
                #this is buggy
                #f.scene.set_size(figsize) 
                #M.clf()
            else:
                raise ValueError('unrecognized figtype')
                
            argsp =  inspect.getargspec(figf)
            if argsp[2] is not None:
                figf(**argd)
            else:
                args=[]
                for arg in argsp[0]:
                    args.append(argd[arg])
                figf(*args)
            
            if save:
                epsconv = False
                if ftype=='mpltoeps':
                    epsconv = 'png'
                    savefunc = plt.savefig
                    saveexts = ['png']
                elif ftype.startswith('mpl'):
                    savefunc = plt.savefig
                    saveexts = ftype[3:].strip()
                    if saveexts=='':
                        if isinstance(savetypes,basestring):
                            saveexts = [savetypes]
                        else:
                            saveexts = savetypes
                    else:
                        saveexts = saveexts.split(',')
                elif ftype.startswith('mayavi'):
                    savefunc = M.savefig
                    saveexts = ftype[6:].strip()
                    if saveexts=='':
                        saveexts = ['png']
                    else:
                        saveexts = saveexts.split(',')
                        if 'eps' in saveexts:
                            saveexts.remove('eps')
                            if 'png' not in saveexts:
                                saveexts.append('png')
                            epsconv = 'png'
                else:
                    raise ValueError('unrecognized figtype')
                
                epsfns = []
                for saveext in saveexts:
                    fn = '%s%s.%s'%(save,fname,saveext)
                    
                    if not overwrite and os.path.exists(fn):
                        print fn,'exists, skipping (use -o at command line to overwrite)'
                    else:
                        print 'saving',fn
                        savefunc(fn)
                        if saveext == epsconv:
                            epsfns.append((fn,'%s%s.%s'%(save,fname,'eps')))
                        elif saveexts == 'eps':
                            epsfns.append((False,fn))
                
                for epsfnt in epsfns:
                    from os import system
                    fn,epsfn = epsfnt
                    if fn:
                        print 'converting',fn,'to',epsfn
                        retcode = subprocess.call('convert %s %s'%(fn,epsfn),shell=True)
                        if retcode is not 0:
                            raise ValueError('convert failed to convert %s to %s'%(fn,epsfn))
                    
                    if tweakbbox is not None:
                        print 'updating eps bounding box on',epsfn,'to',tweakbbox
                        with open(epsfn+'.tmpwrite','w') as fw:
                            with open(epsfn,'r') as fr:
                                for l in fr:
                                    if 'BoundingBox' in l and 'Page' not in l:
                                        newline = l.split(':')[0]+(': %s %s %s %s\n'%tweakbbox)
                                        fw.write(newline)
                                    else:
                                        fw.write(l)
                        shutil.move(epsfn+'.tmpwrite',epsfn)
                    
    finally:
        if isinter:
            plt.ion()
    if showfigs and anympl:
        plt.draw()


#<-----ADD CLASSES/FUNCTIONS HERE - decorate plot functions with @plotfunc----->

if __name__ == '__main__':
    from optparse import OptionParser
    import sys
    
    #ADD ON-RUN OPERATIONS HERE - i.e. load data files, name variables, etc.
    
    op = OptionParser()
    op.usage = '%prog [options] [plot1 plot2 ...]|[maintype]'
    
    op.add_option('-m','--main-figs',help='Shows all main figures.  If an argument is given, it tells which group of main figures to show.',
                        dest='mainfigs',action='store_true',default=False)
    op.add_option('-s','--save',help='save in the specified location',default=None)
    op.add_option('-o','--overwrite',help='Overwrite existing figure',action='store_true',default=False)
    
    ops,args = op.parse_args()
    
    if ops.mainfigs:
        if len(args)==0:
            figstomake = _mainfigs[True][:]
            if len(figstomake)>0:    
                make_plots(set(figstomake),save=ops.save,overwrite=ops.overwrite)
        else:
            for arg in args:
                if arg not in _mainfigs:
                    print 'Invalid main figure group',arg,'Exiting...'
                    sys.exit(1)
            for arg in args:
                figstomake = _mainfigs[True][:]
                figstomake.extend(_mainfigs[arg])
                if ops.save is None:
                    save = arg
                else:
                    save = ops.save
                if len(figstomake)>0:
                    make_plots(set(figstomake),save=save,overwrite=ops.overwrite)
        
        
    else:
        figstomake = args
        if len(figstomake)>0:    
            make_plots(set(figstomake),save=ops.save,overwrite=ops.overwrite)
        
        
