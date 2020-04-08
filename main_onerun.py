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
    
    x = {'u_shift_mar': 1.2153751773920431,
         'sigma_psi': 0.35129329404362025,
         'sigma_psi_mult': 2.0296007012130284,
         'pmeet_0': 0.9597288580956033,
         'pmeet_t': -0.04040761001542115,
         'pmeet_t2': -0.005449583132069378,
         'util_alp': 1.558283511068975,
         'util_kap': 9.62860689159476,
         'preg_a0': 0.012275448802295209,
         'preg_at': -0.02739241487331874,
         'util_qbar': 2.9271875945093018,
         'disutil_marry_sm_mal_coef': 11.391239869986203}


    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))