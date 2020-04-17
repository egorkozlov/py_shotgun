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
    

    x = {'u_shift_mar': 1.210628343658612,
         'sigma_psi': 0.26639709030019953,
         'sigma_psi_mult': 1.8248012429189604,
         'pmeet_0': 0.844167184740344,
         'pmeet_t': -0.04230335372769154,
         'pmeet_t2': -0.0040029283292030565,
         'util_alp': 2.421656979989644,
         'util_kap': 7.410375667727777,
         'util_lam': 0.5459630550173283,
         'preg_a0': 0.08337096159543468,
         'preg_at': -0.006940437283707093,
         'preg_at2': -0.008392522495184568,
         'util_qbar': 1.6547172104823515,
         'disutil_marry_sm_mal_coef': 14.365005988625324,
         'disutil_shotgun_coef': 1.890394571999703}

    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))