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
     
    x = {'u_shift_mar': 1.513283479062665,
     'sigma_psi': 0.5843113309708895,
     'sigma_psi_mult': 3.8226850007122897,
     'pmeet': 0.3271653553125473,
     'pmeet_t': -0.008108358997352225,
     'util_alp': 0.7457467041652874,
     'util_kap': 1.3223533546970248,
     'preg_20': 0.03980279614565633,
     'preg_30': 0.11395486976467334,
     'poutsm': 0.125,
     'util_qbar': 2.8312530407001595}

    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))