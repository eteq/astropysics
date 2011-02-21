
if __name__=='__main__':
    from nose import run,config
    from nose.plugins import doctests,testid
    import os,sofatests
    
    print 'Running nose tests'
    cfg = config.Config(verbosity=2)
    run(config=cfg,plugins=[testid.TestId(),doctests.Doctest()])
    
    print 'Completed nose tests, running SOFA-based tests'
    if not os.path.abspath(os.curdir).endswith('tests'):
        os.chdir(os.path.join(os.curdir,'tests'))
    sofatests.main(compile=True,epoch=None)

