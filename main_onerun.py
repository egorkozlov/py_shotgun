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
    
    
    
    
    x = {'sigma_psi': 0.2574534718663512,
         'sigma_psi_mult': 4.843521651707444,
         'pmeet_0': 0.4118969594789825,
         'pmeet_t': 0.02597917916494849,
         'pmeet_t2': -0.0014532220037343849,
         'preg_a0': 0.08470713450911965,
         'preg_at': -0.0021214797571710853,
         'preg_at2': -0.0013603726333735015,
         'u_shift_mar': 2.141116670945274,
         'util_alp': 0.7548334306376878,
         'util_kap': 0.7468838375304478,
         'util_qbar': 1.4138765424769943,
         'disutil_marry_sm_mal_coef': 13.993615199797048,
         'disutil_shotgun_coef': 0.8626166830980367}


    






    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))
    
    from fit_plot import make_fit_plots
    make_fit_plots(mdl[0].setup,agents.compute_moments(),target_values('high education'))