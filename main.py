#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""





if __name__ == '__main__':
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass


from platform import system
    
if system() != 'Darwin' and system() != 'Windows':   
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'


import numpy as np
from numpy.random import random_sample as rs
from tiktak import tiktak
print('Hi!')


from residuals import mdl_resid
from calibration_params import calibration_params

if __name__ == '__main__':
    
    fix_values = True
    if fix_values:
        xfix = {'u_shift_mar': 1.691189,
                'sigma_psi': 0.4584274,
                'sigma_psi_mult': 3.79805726,
                'util_alp': 0.98122074,
                'util_kap': 0.92348031}
    else:
        xfix = None
        
    
        
    lb, ub, xdef, keys, translator = calibration_params(xfix=xfix)
    
    print('calibration adjusts {}'.format(keys))
    
    print('')
    print('')
    print('running tic tac...')
    print('')
    print('')
    
    
    
    #Tik Tak Optimization
    param=tiktak(N=500,N_st=20,skip_local=True,skip_global=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
    #Now Re do the computation with graphs!
    out, mdl = mdl_resid(param[1],return_format=['distance','model'],
                         verbose=True,draw=True)
    
    
   
        

