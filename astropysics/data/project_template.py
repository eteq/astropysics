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

#<--------------------Plot support functions----------------------------------->
_plotreg = {}
_mainfigs = []
    
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
        a paper (used for the command-line call)
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
            _mainfigs.append(fname)
        return f
    
    if callable(figtype):
        f=figtype
        figtype='mpl'
        return deco(f)
    else:
        return deco
    
def plot_names():
    """
    Returns the names of all registered plots for this file, as well as the
    names of just the 'main' plots.
    
    :returns: 2-tuple (plotnames,mainnames)
    
    """
    return _plotreg.keys(),_mainfigs

def make_plots(plots,argd=None,figs=None,save=False,overwrite=True,showfigs=True):
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
    """
    import os,shutil,inspect,subprocess   
    
    if save is True:
        save='.'
    if type(save) is str and not save.endswith(os.sep):
        save+=os.sep        
    
    if plots is None:
        plots = _plotreg.keys()
    elif isinstance(plots,basestring):
        if '*' in plots:
            import re
            plots = plots.replace('*','.*')
            regex = re.compile(plots)
            plots = [p for p in plot_names() if regex.match(p)]
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
                        saveexts = ['eps','pdf']
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
                
                for saveext in saveexts:
                    fn = '%s%s.%s'%(save,fname,saveext)
                    epsfns = None
                    
                    if not overwrite and exists(fn):
                        print fn,'exists, skipping (use -o at command line to overwrite)'
                    else:
                        print 'saving',fn
                        savefunc(fn)
                        if saveext == epsconv:
                            epsfns = (fn,'%s%s.%s'%(save,fname,'eps'))
                
                if epsfns is not None:
                    from os import system
                    fn,epsfn = epsfns
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

def main():
    try: 
        from argparse import ArgumentParser
        p = ArgumentParser()
        
        p.add_argument('-m','--main-figs',help='Shows all main figures even if they are not included in command line',
                            dest='mainfigs',action='store_true',default=False)
        p.add_argument('-s','--save',help='save in the specified location',default=None)
        p.add_argument('-o','--overwrite',help='Overwrite existing figure',action='store_true',default=False)
        
        add_to_parser(p)
        args = p.parse_args()
        on_run(args)
        
        if args.mainfigs:
            figstomake = _mainfigs[:]
            figstomake.extend(args.plots)
        else:
            figstomake = args.plots
        
        save = args.save
        overwrite = args.overwrite
        
    except ImportError: #argparse not present
        from optparse import OptionParser
        
        op = OptionParser()
        op.usage = '%prog [options] [plot1 plot2 ...]'
        
        op.add_option('-m','--main-figs',help='Shows all main figures even if they are not included in command line',
                            dest='mainfigs',action='store_true',default=False)
        op.add_option('-s','--save',help='save in the specified location',default=None)
        op.add_option('-o','--overwrite',help='Overwrite existing figure',action='store_true',default=False)
        
        add_to_parser(op)
        ops,args = op.parse_args()
        on_run((ops,args))
        
        if ops.mainfigs:
            figstomake = _mainfigs[:]
            figstomake.extend(args)
        else:
            figstomake = args
            
        save = ops.save
        overwrite = ops.overwrite
    
    if len(figstomake)>0:    
        make_plots(set(figstomake),save=save,overwrite=overwrite)

#<-------------------------THESE CAN BE CUSTOMIZED----------------------------->
def add_to_parser(p):
    """
    Runs before arguments are parsed, to support adding command-line arguments.
    
    :param p: 
        A :class:`argparse.ArgumentParser` object unless your python version 
        doesn't support it, in which case it is a
        :class:`optparse.OptionParser` object.
    """

def on_run(args):
    """
    Run additional tasks before plotting begins for this script.
    
    :param args: Either an :class:`argparse.Namespace` object or an (ops,args) 
    tuple, depending on whether or not :mod:`argparse` is supported by your
    python version.
    """

#<--ADD CLASSES/FUNCTIONS BELOW HERE - decorate plot functions with @plotfunc-->



#DON'T DELETE THIS:
if __name__ == '__main__':
    main()