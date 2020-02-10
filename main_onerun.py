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
 
print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
     
    x = {'u_shift_mar': 0.78012885,
         'sigma_psi': 0.22466639,
         'sigma_psi_mult': 1.34982032,
         'pmeet': 0.17557838,
         'pmeet_t': -0.00900235,
         'util_alp': 0.47896648,
         'util_kap': 1.72740958,
         'preg_20': 0.2184985,
         'preg_30': 0.33067173}
    
    out, mdl, agents, res = mdl_resid(x=x,return_format=['distance','models','agents','scaled residuals'],
                                      load_from='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))