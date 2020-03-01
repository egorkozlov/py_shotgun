#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27
 
@author: Egor Kozlov
"""


#import warnings
#warnings.filterwarnings("error")
 

from platform import system
     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
 
 
from residuals import mdl_resid
from targets import target_values

print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
     
    x = {'u_shift_mar': 2.3355210455808297,
         'sigma_psi': 0.6053952015313737,
         'sigma_psi_mult': 3.658469721691887,
         'pmeet': 0.36236559249094785,
         'pmeet_t': -0.010441196032291739,
         'util_alp': 1.724517409748223,
         'util_kap': 1.3477223782840697,
         'preg_20': 0.04812530882453013,
         'preg_30': 0.0900267785184747,
         'sm_shift': -2.812699250993756,
         'util_qbar': 0.5*3.9922428798997034}
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))