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
         'taste_shock_mult': 4.116448914683272}
    

    x = {'sigma_psi': 0.22225792188972474,
         'sigma_psi_mult': 7.461120766010872,
         'pmeet_21': 0.2841209803418755,
         'pmeet_28': 0.5972782973618359,
         'pmeet_35': 0.8495146505392482,
         'preg_21': 0.03351967307448346,
         'preg_28': 0.02658760503060691,
         'preg_35': 0.023934310550573212,
         'u_shift_mar': 1.1813684776362563,
         'util_alp': 0.41580007423909804,
         'util_kap': 0.9018094315778211,
         'util_qbar': 0.9224730817253374,
         'disutil_marry_sm_mal_coef': 9.191322692652617,
         'disutil_shotgun_coef': 0.23934745280098635,
         'taste_shock_mult': 5.224732545765163}


    

    tar = target_values('high education')
    
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True,cs_moments=False,
                                      moments_repeat=2)
    
    mdl[0].time_statistics()
                         
    
    print('Done. Residual in point x0 is {}'.format(out))
    
    
    #from simulations import Agents
    #moments_aux = Agents( mdl, N=10000, T=18, female=False, verbose=False).aux_moments()
    from fit_plot import make_fit_plots
    make_fit_plots(agents,target_values('high education'))#,moments_aux=moments_aux)
    
    
    '''
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
    '''