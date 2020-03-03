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
     
    x = {'u_shift_mar': 2.5044352358268593,
         'sigma_psi': 0.5222876006269004,
         'sigma_psi_mult': 3.3001050099184672,
         'pmeet': 0.2867917284023725,
         'pmeet_t': -0.014152767736645006,
         'util_alp': 0.635329927428283,
         'util_kap': 3.948202011966941,
         'preg_20': 0.0368466345061881,
         'preg_30': 0.0611339918872953,
         'util_qbar': 1.5720864281376485,
         'disutil_marry_sm_mal_coef': 4.656512556378265}
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))