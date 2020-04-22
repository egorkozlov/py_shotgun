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
    

    x = {'sigma_psi': 0.39334993551757474,
         'sigma_psi_mult': 6.461534693482925,
         'pmeet_0': 0.3988891928888862,
         'preg_a0': 0.06072186009040158,
         'preg_at': -0.044917971687432,
         'preg_at2': -0.00997429561655438,
         'u_shift_mar': 1.8866425125831014,
         'util_alp': 3.6302729099123443,
         'util_kap': 14.857375158019218,
         'util_qbar': 4.994817208261213,
         'disutil_marry_sm_mal_coef': 5.619415004495844,
         'disutil_shotgun_coef': 0.672932495972583}





    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))