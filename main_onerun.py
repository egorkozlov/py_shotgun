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
    
    
    
    
    x = {'sigma_psi': 0.3158003793955531,
         'sigma_psi_mult': 4.947443937771622,
         'pmeet_0': 0.4118969594789825,
         'pmeet_t': 0.02597917916494849,
         'pmeet_t2': -0.0014532220037343849,
         'preg_a0': 0.08470713450911965,
         'preg_at': -0.0021214797571710853,
         'preg_at2': -0.0013603726333735015,
         'u_shift_mar': 2.000120246942508,
         'util_alp': 0.6122298960907041,
         'util_kap': 0.7671589794930209,
         'util_qbar': 1.7407912648838693,
         'disutil_marry_sm_mal_coef': 11.898399107945592,
         'disutil_shotgun_coef': 0.5124129679921092}








    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))