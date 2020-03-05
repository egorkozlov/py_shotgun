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
     
    x = {'u_shift_mar': 1.2375741262803743,
         'sigma_psi': 0.09346593854048857,
         'sigma_psi_mult': 2.471258656852517,
         'pmeet': 0.17036689841688138,
         'pmeet_t': 0.0,
         'util_alp': 1.1316659121221422,
         'util_kap': 6.9221292953673315,
         'preg_20': 0.08473982500019177,
         'preg_30': 0.01,
         'util_qbar': 1.4656511339567124,
         'disutil_marry_sm_mal_coef': 7.487682392559771}


    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))