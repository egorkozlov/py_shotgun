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
    
    
    
    
    x = {'sigma_psi': 0.29657734335571323,
         'sigma_psi_mult': 3.729950449694295,
         'pmeet_0': 0.4588845109118548,
         'pmeet_t': 0.028520315381868666,
         'pmeet_t2': -0.0027715221438747405,
         'preg_a0': 0.10568137068532876,
         'preg_at': -0.012412628483134985,
         'preg_at2': -0.005036760005378196,
         'u_shift_mar': 2.130767784054183,
         'util_alp': 0.7388571167758351,
         'util_kap': 0.7655243824821856,
         'util_qbar': 1.4951823800236654,
         'disutil_marry_sm_mal_coef': 9.998061340856006,
         'disutil_shotgun_coef': 0.5298630898532276}






    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))