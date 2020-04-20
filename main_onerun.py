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
    

    x = {'sigma_psi': 0.2910544737092481,
         'sigma_psi_mult': 3.539469008535118,
         'pmeet_0': 0.26244808909605194,
         'pmeet_t': 0.1328625646023624,
         'pmeet_t2': -0.002715921755452132,
         'preg_a0': 0.3855181445875888,
         'preg_at': 0.16529108344731835,
         'preg_at2': -0.00404874287820581,
         'u_shift_mar': 1.283362837172244,
         'util_alp': 2.693149863108261,
         'util_kap': 12.286338934623597,
         'util_qbar': 0.6043845729601405,
         'sm_shift': -1.991160314738063,
         'disutil_marry_sm_mal_coef': 14.727637231587778,
         'disutil_shotgun_coef': 3.7946494078610344}



    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))