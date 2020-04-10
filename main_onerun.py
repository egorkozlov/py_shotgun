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
    
    x = {'u_shift_mar': 2.7230034052298824,
         'sigma_psi': 0.17572487206371926,
         'sigma_psi_mult': 4.100699121313946,
         'pmeet_0': 0.5921205737287893,
         'pmeet_t': -0.014010138916286219,
         'pmeet_t2': -0.00139459803708393,
         'util_alp': 1.1170675899324138,
         'util_kap': 6.928722658429808,
         'preg_a0': 0.061387059514552235,
         'preg_at': -0.021028784919175134,
         'preg_at2': 0.0,
         'util_qbar': 2.570238658000373,
         'disutil_marry_sm_mal_coef': 6.62134049339284}
    
    



    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))