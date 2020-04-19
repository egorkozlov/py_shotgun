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
    

    x = {'u_shift_mar': 1.309432986233957,
         'sigma_psi': 0.22100216636978237,
         'sigma_psi_mult': 2.528042541325081,
         'pmeet_0': 0.1470379479378908,
         'pmeet_t': 0.09943956197088111,
         'pmeet_t2': -0.005415529671440785,
         'util_alp': 2.489731956618132,
         'util_kap': 8.697384793595958,
         'util_lam': 0.6000410107481078,
         'preg_a0': 0.0925469387893579,
         'preg_at': -0.03388403739204715,
         'preg_at2': -0.005269467250514311,
         'util_qbar': 3.661546058586023,
         'disutil_marry_sm_mal_coef': 12.646588370071647,
         'disutil_shotgun_coef': 2.0855699509471}


    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))