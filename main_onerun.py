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
    
    x = {'u_shift_mar': 2.454083991840396,
         'sigma_psi': 0.16666295932923875,
         'sigma_psi_mult': 4.2578420743166365,
         'pmeet_0': 0.512603449144626,
         'pmeet_t': -0.00392071562690105,
         'pmeet_t2': -0.002523247366351184,
         'util_alp': 1.0610200316606404,
         'util_kap': 8.104853182716534,
         'preg_a0': 0.10791736157187538,
         'preg_at': -0.0319145736956649,
         'util_qbar': 3.2532324599246265,
         'disutil_marry_sm_mal_coef': 7.441147983430074}



    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      load_from='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))