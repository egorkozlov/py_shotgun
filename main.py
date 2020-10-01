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
from estimates import get_point
print('Hi!')


from residuals import mdl_resid
from calibration_params import calibration_params

if __name__ == '__main__':
    
    
    #from numba import config
    #config.NUMBA_NUM_THREADS = 2

    
    fix_values = False
    if fix_values:
        
        keys_fix = ['sigma_psi','sigma_psi_init','mu_psi_init',
                    'u_shift_mar','util_alp','util_kap']
        
        x_high, _ = get_point(True,read_wisdom=False)
        xfix = {k : x_high[k] for k in keys_fix}
        
        
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
    param=tiktak(xfix=xfix,N=6000,N_st=100,skip_local=False,skip_global=True,
                             resume_global=False,resume_local=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
   
        

