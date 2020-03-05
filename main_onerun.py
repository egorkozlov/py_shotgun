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
     
    x = {'u_shift_mar': 1.4622468191770523,
         'sigma_psi': 0.13341343588771484,
         'sigma_psi_mult': 2.7114359061341653,
         'pmeet': 0.16330434124234866,
         'pmeet_t': 0.0,
         'util_alp': 1.5035902442038718,
         'util_kap': 7.862043282658021,
         'preg_20': 0.10034665087728878,
         'preg_30': 0.01,
         'util_qbar': 1.7158729365532353,
         'disutil_marry_sm_mal_coef': 8.047422930867116}

    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))