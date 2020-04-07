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
    
    x = {'u_shift_mar': 1.1236372092858131,
     'sigma_psi': 0.19039234194260712,
     'sigma_psi_mult': 6.125171460954793,
     'pmeet_0': 0.1602430213239973,
     'pmeet_t': 0.054536185185160437,
     'pmeet_t2': -0.0013660983035121872,
     'util_alp': 1.063925000236826,
     'util_kap': 7.096103099395396,
     'preg_a0': 0.16472997866861458,
     'preg_at': -0.05333624261834616,
     'util_qbar': 1.786295966229198,
     'disutil_marry_sm_mal_coef': 12.41904957062656}

    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))