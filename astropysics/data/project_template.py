#!/usr/bin/env python
"""
<author> <date>

<project description>

<DELETE COMMENTS BELOW HERE AFTER READING>
This file is a template for a data analysis project (typically to go with a 
paper or two).  The main usage scenario is to plass functions and classes 
before the "if __name__=='__main__' line.  Plots/figures go there as well,
with the @plotfunc decorator.  Then the command line usage can be used to 
re-generate any/all plots, and interactive work in ipython is supported by doing:

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

Replace "ADD CLASSES/FUNCTIONS HERE" with:

@plotfunc(mainfig=True)
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
    
and put the following in the "on_run" function:

    x = linspace(0,1,100)
    data1 = x*2.5 + 2
    data2 = data1 + randn(100)
    
    
then you could do:

./myscript.py -m

at the command line, and figues 'aplot.eps' and 'aplot.pdf' will be saved in 
the current directory.

Also, if you replace the "@plotfunc(mainfig=True)" with 
"@plotfunc(mainfig='paper1')", and call it as:

./myscript.py -m paper1

A new directory "paper1" will be created with the 'aplot.eps' and 'aplot.pdf'
figures.


Finally, keep in mind that this script was written to be used with ipython.
Be sure to start ipython as "ipython -pylab" (or "ipython --pylab" if you have
an old version of ipython). You can then just do 
"run myscript.py <command line args>" and it will load all your variables from
on_run into the ipython namespace.  You can then "import myscript", and 
whenever you want to try to re-make a plot, either because you changed the 
value of some variable or because you changed the plot function, just do:

reload(myscript)
myscript.make_plot('plotname',locals())

And your new plot should appear with your changes and using the value of any
of the plot functions arguments from the same-name variable in ipython.

"""
from __future__ import division,with_statement

from numpy import *
from numpy.random import *
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from collections import defaultdict

#<--------------------Plot support functions----------------------------------->
_plotreg = {}
_mainfigs = defaultdict(list)
    
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
        a paper (if '-m' is used as a command line argument).  If it is a 
        string, this figure will be registered as a figure to be generated when 
        the '-m string' syntax is used as a command line argument.
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
            _mainfigs[mainfig].append(fname)
        return f
    
    if callable(figtype):
        f=figtype
        figtype='mpl'
        return deco(f)
    else:
        return deco
    
def plot_names(printnms=False):
    """
    Returns the names of all registered plots for this file, as well as the
    names of just the 'main' plots.
    
    :params printnms: If True, the plot names will be printed.
    
    :returns: 
        2-tuple (plotnames,mainnamedict) where `mainnamedict` is a dict mapping
        main group names to main figure names (the 'True' key maps to those
        without a group name).
    
    """
    if printnms:
        if len(_mainfigs)>0:
            mainnms = []
        
            print 'Main figure names:'
            for n in _mainfigs[True]:
                print n
                mainnms.append(n)
            for grp in _mainfigs:
                if grp != True:
                    print 'Group "%s":'%grp
                    for n in _mainfigs[grp]:
                        print n
                        mainnms.append(n)
                
            if len(_plotreg)>len(mainnms):
                print '\nOther plot names:'
                for n in _plotreg.keys():
                    if n not in mainnms:
                        print n
        elif len(_plotreg)>0:
            print 'Plot names:'
            for n in _plotreg.keys():
                print n
        else:
            print 'No plots defined.'
        
    return _plotreg.keys(),dict(_mainfigs)

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
    if isinstance(save,basestring) and not save.endswith(os.sep):
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
                    
                    if not overwrite and os.path.exists(fn):
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
    """
    Runs when the function is called from the command line.  Don't alter this 
    function directly - instead, edit :func:`on_run` and :func:`add_to_parser`
    
    :returns: A dictionary of the local variables from :func:`on_run`.
    """
    try: 
        from argparse import ArgumentParser
        p = ArgumentParser()
        
        p.add_argument('-l','--list',help='Show plot names and exit (overrides all other arguments)',action='store_true',default=False)
        p.add_argument('-m','--main-figs',dest='mainfigs',const=True,default=False,nargs='?',
            help='Shows all main figures even if they are not included in the command line.'+
                 ' Can also provide a mainfig group name to create only mainfigs associated with that group.')
        p.add_argument('-s','--save',help='save in the specified location',default=None)
        p.add_argument('-o','--overwrite',help='Overwrite existing figure',action='store_true',default=False)
        p.add_argument('plots',metavar='pltnm',nargs='*',help='Names of plots to show')
        
        add_to_parser(p)
        args = p.parse_args()
        
        if args.list:
            plot_names(True)
            return {}
                
        rvars = on_run(args)
        
        if args.mainfigs:
            if args.mainfigs==True:
                figstomake = []
                for fs in _mainfigs.values():
                    figstomake.extend(fs)
            else:
                if args.mainfigs not in _mainfigs:
                    raise KeyError('Main group %s not found!'%args.mainfigs)
                figstomake = _mainfigs[args.mainfigs][:]
            figstomake.extend(args.plots)
        else:
            figstomake = args.plots
        
        save = args.save
        overwrite = args.overwrite
        
    except ImportError: #argparse not present
        from optparse import OptionParser
        
        op = OptionParser()
        op.usage = '%prog [options] [plot1 plot2 ...]'
        
        op.add_option('-l','--list',help='Show plot names and exit (overrides all other arguments)',action='store_true',default=False)
        op.add_option('-m','--main-figs',dest='mainfigs',action='store_true',default=False,
            help='Shows all main figures even if they are not included in the command line.'+
                 ' Can also provide a mainfig group name to create only mainfigs associated with that group.')
        op.add_option('-s','--save',help='save in the specified location',default=None)
        op.add_option('-o','--overwrite',help='Overwrite existing figure',action='store_true',default=False)
        
        add_to_parser(op)
        ops,args = op.parse_args()
        
        if ops.list:
            plot_names(True)
            return {}
        
        rvars = on_run((ops,args))
        
        if ops.mainfigs:
            try:
                mainfigrp = args.pop(0)
                if mainfigrp not in _mainfigs:
                    raise KeyError('Main group %s not found!'%mainfigrp)
                figstomake = _mainfigs[mainfigrp][:]
            except IndexError:
                figstomake = []
                for fs in _mainfigs.values():
                    figstomake.extend(fs)
                
            figstomake.extend(args)
        else:
            figstomake = args
            
        save = ops.save
        overwrite = ops.overwrite
    
    if len(figstomake)>0:
        indict = rvars.copy()
        indict.update(globals())
        make_plots(set(figstomake),indict,save=save,overwrite=overwrite)
        
    return rvars

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
    
    Any local variables in this function will end up in *global* scope when 
    this script is run, so this is where you should define the variables that
    the plots are going to use.
    
    :param args: Either an :class:`argparse.Namespace` object or an (ops,args) 
    tuple, depending on whether or not :mod:`argparse` is supported by your
    python version.
    
    :returns: A dictionary of all the local variables created in this function.
    """
    return locals()

#<----ADD YOUR OWN CODE BELOW HERE - decorate plot functions with @plotfunc---->


#DON'T DELETE THIS:
if __name__ == '__main__':
    globals().update(main())
