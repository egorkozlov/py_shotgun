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
    
    
    #from numba import config
    #config.NUMBA_NUM_THREADS = 2

    
    fix_values = True
    if fix_values:
        xfix = {'sigma_psi': 0.2798538779149558,
                'sigma_psi_init': 1.6916827480353398,
                'u_shift_mar': 1.1207829780316174,
                'util_alp': 0.33532441265560753,
                'util_kap': 0.8617581929517057, 
                'couple_rts': 2.80540188195949
                }
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
    param=tiktak(xfix=xfix,N=6000,N_st=250,skip_local=False,skip_global=False,
                             resume_global=False,resume_local=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
    #Now Re do the computation with graphs!
    out, mdl = mdl_resid(param[1],return_format=['distance','model'],
                         verbose=True,draw=True)
    
    
   
        

