#!/usr/bin/env python
from __future__ import division,with_statement

import sys,os,difflib
from astropysics import publication

def test_tex_file(fn=None,htmldiff=None):
    if fn is None:
        testdir = os.path.split(os.path.realpath(__file__))[0]
        fn = testdir + os.sep + 'test_pub.tex'
    with open(fn) as f:
        orig = f.read().split('\n')
    tfstr = publication.TeXFile(fn)().split('\n')
    
    diffres = list(difflib.unified_diff(orig,tfstr))
    if htmldiff:
        hd = difflib.HtmlDiff()
        htmlres = hd.make_file(orig,tfstr)
        with open(htmldiff,'w') as f:
            f.write(htmlres)
    else:
        for l in diffres:
            print l
    assert len(diffres)==0,'TeXFile string does not match original file '+os.path.split(fn)[1]
        
if __name__ == '__main__':
    if len(sys.argv) > 1:
        test_tex_file(htmldiff=sys.argv[1])
    else: 
        test_tex_file()
