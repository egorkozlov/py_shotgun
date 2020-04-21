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
    

    x = {'sigma_psi': 0.36748173935405126,
         'sigma_psi_mult': 6.235752108631493,
         'pmeet_0': 0.3471712860271766,
         'preg_a0': 0.08712286357271183,
         'preg_at': -0.03072066904819082,
         'preg_at2': -0.008558694160374839,
         'u_shift_mar': 1.1888100917273723,
         'util_alp': 0.9489361772874477,
         'util_kap': 4.314946848831633,
         'util_qbar': 0.3044538970955245,
         'disutil_marry_sm_mal_coef': 0.7000150539157969,
         'disutil_shotgun_coef': 3.8177626705415153}





    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))