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
    
    x = {'u_shift_mar': 1.125586011451983,
         'sigma_psi': 0.2836129982791315,
         'sigma_psi_mult': 4.472677709835214,
         'pmeet_0': 0.3531691856924769,
         'pmeet_t': -0.002789200924016383,
         'pmeet_t2': -0.00015158750562126004,
         'util_alp': 1.6068900506337755,
         'util_kap': 11.438222382035043,
         'util_lam': 0.7,
         'util_xi': 1.75,
         'preg_a0': 0.1378677816718491,
         'preg_at': -0.0205690435756406,
         'preg_at2': -0.009291712345759809,
         'util_qbar': 1.1779594356137741,
         'couple_rts': 0.17659154224613074,
         'disutil_marry_sm_mal_coef': 11.481636871149822}


    
    



    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))