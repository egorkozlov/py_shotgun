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
     
    x = {'u_shift_mar': 0.725863356303563,
         'sigma_psi': 0.13579410386649068,
         'sigma_psi_mult': 4.270601078298901,
         'pmeet': 0.30889908284645357,
         'pmeet_t': -0.0038829067012848967,
         'util_alp': 0.6534190912803465,
         'util_kap': 1.9136130954048896,
         'preg_20': 0.03915073600027051,
         'preg_30': 0.10173830697979773,
         'targets':'high education'}
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,return_format=['distance','models','agents','scaled residuals'],
                                      load_from='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))