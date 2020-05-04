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
    
    
    '''
    x = {'sigma_psi': 0.4042887868032187,
         'sigma_psi_mult': 5.159831214740526,
         'pmeet_21': 0.09462288438881994,
         'pmeet_28': 0.3497836550541352,
         'pmeet_35': 0.45958239075216556,
         'preg_21': 0.020184941653201283,
         'preg_28': 0.040679132216882644,
         'preg_35': -0.00613815415792901,
         'u_shift_mar': 1.9391514168600832,
         'util_alp': 0.6740601543545901,
         'util_kap': 0.8245547047683262,
         'util_qbar': 0.494312244622116,
         'disutil_marry_sm_mal_coef': 14.575313199536971,
         'disutil_shotgun_coef': 0.5747914283838933}
    '''
    x = {'sigma_psi': 0.3992869835377583,
         'sigma_psi_mult': 5.072185498196008,
         'pmeet_21': 0.0859842703527165,
         'pmeet_28': 0.3700582376952367,
         'pmeet_35': 0.47137306696712405,
         'preg_21': 0.02232032891602287,
         'preg_28': 0.04157025892605673,
         'preg_35': 0.028035339707370655,
         'u_shift_mar': 1.836551924933839,
         'util_alp': 0.6271358861497329,
         'util_kap': 0.8341827838900215,
         'util_qbar': 0.0,
         'disutil_marry_sm_mal_coef': 15.500831380402422,
         'disutil_shotgun_coef': 0.6284153641357494}


    

    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #load_from='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))
    
    from fit_plot import make_fit_plots
    make_fit_plots(agents,target_values('high education'))