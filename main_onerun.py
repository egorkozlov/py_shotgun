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
    x = {'sigma_psi': 0.4203961804901682,
         'sigma_psi_mult': 4.774783102380805,
         'pmeet_21': 0.1150734338507825,
         'pmeet_28': 0.3658702881000681,
         'pmeet_35': 0.49328894668466444,
         'preg_21': 0.026909130608297693,
         'preg_28': 0.03678288557173975,
         'preg_35': -0.010004147713862654,
         'u_shift_mar': 2.037465051515264,
         'util_alp': 0.7158378198840009,
         'util_kap': 0.826411106925244,
         'util_qbar': 0.7521355119903652,
         'disutil_marry_sm_mal_coef': 14.595453417926215,
         'disutil_shotgun_coef': 0.33579823302000383}
    '''
    x = {'sigma_psi': 0.37621958,
         'sigma_psi_mult': 5.09055414,
         'pmeet_21': 0.11759438,
         'pmeet_28': 0.36227985,
         'pmeet_35': 0.47577083,
         'preg_21': 0.04601108,
         'preg_28': 0.03504144,
         'preg_35': 0.0610367,
         'u_shift_mar': 1.94371434,
         'util_alp': 0.68140333,
         'util_kap': 0.82638585,
         'util_qbar': 0.67972182,
         'disutil_marry_sm_mal_coef': 13.89400622,
         'disutil_shotgun_coef': 0.49459884}


    

    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))
    
    from fit_plot import make_fit_plots
    make_fit_plots(agents,target_values('high education'))