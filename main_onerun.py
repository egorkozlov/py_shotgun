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
    
    x = {'u_shift_mar': 3.491556927367782,
         'sigma_psi': 0.23125088678174552,
         'sigma_psi_mult': 4.190695479590094,
         'pmeet_0': 0.6909755543788532,
         'pmeet_t': -0.03813741282849738,
         'pmeet_t2': -0.0007523304237593556,
         'util_alp': 1.6618397255153359,
         'util_kap': 7.228897081332338,
         'preg_a0': 0.21927455798827009,
         'preg_at': -0.051015156827834196,
         'preg_at2': -0.004000349359161717,
         'util_qbar': 2.440773537269905,
         'disutil_marry_sm_mal_coef': 10.211476319923747}

    
    



    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))