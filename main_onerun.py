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
    
    
    
    
    x = {'sigma_psi': 0.15708192700193488,
         'sigma_psi_mult': 5.305288050025261,
         'pmeet_0': 0.37323350562593977,
         'pmeet_t': 0.021980702978686817,
         'pmeet_t2': -0.0013690241411399573,
         'preg_a0': 0.05073239657355087,
         'preg_at': 0.0029740693848109845,
         'preg_at2': -0.0006947723476003372,
         'u_shift_mar': 1.6483875285794014,
         'util_alp': 0.6790146096417299,
         'util_kap': 0.7644279343790064,
         'util_qbar': 0.08356166794865397,
         'disutil_marry_sm_mal_coef': 4.779135319968562,
         'disutil_shotgun_coef': 1.682061613383758}





    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=None,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))