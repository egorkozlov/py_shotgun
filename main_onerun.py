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
         'pmeet_21': 0.06037336,
         'pmeet_28': 0.35412571,
         'pmeet_35': 0.50546231,
         'preg_21': -0.00638973,
         'preg_28': 0.0835086,
         'preg_35': 0.04009042,
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
    
    from fit_plot import make_fit_plots
    make_fit_plots(mdl[0].setup,agents.compute_moments(),target_values('high education'))