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
    
    
    # high educ: baseline
    
    
    for npsi in range(7,41):
        
        x = {'sigma_psi': 0.2476722655689144,
             'sigma_psi_mult': 5.125801752198881,
             'pmeet_21': 0.13090071390433766,
             'pmeet_28': 0.38392851771344155,
             'pmeet_35': 0.40800438661571153,
             'preg_21': 0.07677173244887601,
             'preg_28': 0.06061469851763078,
             'preg_35': 0.01825557056243586,
             'u_shift_mar': 1.7329041070973545,
             'util_alp': 0.6182672481649074,
             'util_kap': 0.8081836080864513,
             'util_qbar': 0.5163163798943308,
             'disutil_marry_sm_mal_coef': 14.446603934890161,
             'disutil_shotgun_coef': 0.4904309002252879,
             'taste_shock_mult': 4.116448914683272,
             'high education':True,
             'n_psi':npsi}
        
        
        tar = target_values('high education')        
        out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                          return_format=['distance','models','agents','scaled residuals','moments'],
                                          #load_from='mdl.pkl',
                                          verbose=False,draw=True,cs_moments=False,
                                          moments_repeat=2)
        
        mdl[0].time_statistics()
        print('Done. Residual in point x0 is {}'.format(out))
        del((out, mdl, agents, res, mom))