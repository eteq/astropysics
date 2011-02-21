
if __name__=='__main__':
    from nose import main,config
    from nose.plugins import doctests,testid
    
    cfg = config.Config(verbosity=2)
    main(config=cfg,plugins=[testid.TestId(),doctests.Doctest()])
    

