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
     
    x = {'u_shift_mar': 2.0415312646276234,
         'sigma_psi': 0.3424730326985349,
         'sigma_psi_mult': 3.3858306398842473,
         'pmeet': 0.4247289616489994,
         'pmeet_t': -0.004604640269509339,
         'util_alp': 1.6320658827348236,
         'util_kap': 1.4259194141270437,
         'preg_20': 0.058445449688699594,
         'preg_30': 0.04351737584566795,
         'sm_shift': 0.09515539741538825,
         'util_qbar': 2.0898404131498056,
         'm_zf': 2.6016884179397683,
         'm_zf0': 2.057811638685381}
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))