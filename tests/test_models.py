#!/usr/bin/env python
from __future__ import division,with_statement
import numpy as np
from astropysics import models

def test_1d_model_funcs():    
    for modname in models.list_models(baseclass=models.FunctionModel1D):
        cls = models.get_model_class(modname)
        # ignore datacentric models because they are data-dependant
        if not issubclass(cls,models.DatacentricModel1D):
            #test both ways of generating the model instance
            if models.get_model_class(modname).isVarnumModel():
                #just try it w/ 3 variable arguments w/ default values
                mod = cls(3,)
                mod = models.get_model_instance(modname,nvarparams=3)
            else:
                mod = cls()
                mod = models.get_model_instance(modname)
                
            if mod.rangehint is None:
                x1 = 0
                x2 = 01
            else:
                x1,x2 = mod.rangehint
               
            avgval = mod((x1+x2)/2) #single value output
            rangevals = mod(np.linspace(x1,x2,12)[1:-1]) #array output
            
            assert np.all(np.isfinite(avgval)),'Non-finite value encountered for model '+modname
            assert np.all(np.isfinite(rangevals)),'Non-finite value encountered for model '+modname
            
        
if __name__ == '__main__':
    import nose
    nose.main()
