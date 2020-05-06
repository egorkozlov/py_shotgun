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
    
    out, mdl, agents, res = mdl_resid(x=None,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      save_to='mdl.pkl',
                                      verbose=True,draw=True,cs_moments=False,
                                      moments_repeat=2)
    
    mdl[0].time_statistics()
                         
    
    print('Done. Residual in point x0 is {}'.format(out))
    
    
    
    
    from fit_plot import make_fit_plots
    make_fit_plots(agents,target_values('high education'))
    
    from crosssection import CrossSection
    from simulations import Agents
    mom_list_cs = []
    mom_list_pl = []
    n_repeat = 3
    import numpy as np
    np.random.seed(13)
    for _ in range(n_repeat):
        cs = CrossSection(mdl,fix_seed=False,verbose=False,N_total=30000)
        mom_list_cs = mom_list_cs + [cs.compute_moments().copy()]
        del(cs)
        
        pl = Agents(mdl,N=15000,T=22,verbose=False,fix_seed=False)
        mom_list_pl = mom_list_pl + [pl.compute_moments().copy()]
        del(pl)
        
        
    moments_cs = dict()
    moments_pl = dict()
    keys = mom_list_cs[0].keys()
    keys_pl = mom_list_pl[0].keys()
    
    for key in keys:
        if key not in keys_pl: continue
        moments_cs[key] = np.sum([m[key] for m in mom_list_cs])/n_repeat
        moments_pl[key] = np.sum([m[key] for m in mom_list_pl])/n_repeat
    
    moments_comp = {key: (moments_cs[key],moments_pl[key]) for key in moments_cs}
    