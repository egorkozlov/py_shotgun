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
from targets import target_values, all_targets

print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
    
    high_e = False
    
    if high_e:
        x = {'sigma_psi': 0.13140335512127152,
             'sigma_psi_mult': 7.318494317792161,
             'pmeet_21': 0.31499861958165964,
             'pmeet_28': 0.60067071370404,
             'pmeet_35': 0.7472621428312332,
             'preg_21': 0.040420105606407006,
             'preg_28': 0.02296063935559003,
             'preg_35': 0.0698790264219899,
             'u_shift_mar': 0.9606447172924117,
             'util_alp': 0.3468964692327379,
             'util_kap': 0.9378088317753515,
             'util_qbar': 0.9082817629134391,
             'disutil_marry_sm_mal_coef': 8.44506597458977,
             'disutil_shotgun_coef': 0.25841666575918587,
             'taste_shock_mult': 5.476599495906876,
             'high education':True}
        
        targ_mode = 'high education'
    else:
        x = {'sigma_psi': 0.480991547005569,
             'sigma_psi_mult': 4.600215309475741,
             'pmeet_21': 0.48886846739080425,
             'pmeet_28': 0.17451514254902903,
             'pmeet_35': 0.5918506163621335,
             'preg_21': 0.3960300615120441,
             'preg_28': 0.115344941557771,
             'preg_35': 0.5167972159585046,
             'u_shift_mar': 1.9291671160578177,
             'util_alp': 0.47953144662272207,
             'util_kap': 0.8247164477867641,
             'util_qbar': 0.015315764308810709,
             'disutil_marry_sm_mal_coef': 7.70003003525308,
             'disutil_shotgun_coef': 0.4186198614019019,
             'taste_shock_mult': 3.1008647867466848,
             'high education':False}
        targ_mode = 'low education'




    
    tar = target_values(targ_mode)
    
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      #load_from='mdl.pkl',
                                      verbose=True,draw=True,cs_moments=False,
                                      moments_repeat=2)
    
    mdl[0].time_statistics()
                         
    

    
    print('Done. Residual in point x0 is {}'.format(out))
    
    
    #from simulations import Agents
    #moments_aux = Agents( mdl, N=10000, T=18, female=False, verbose=False).aux_moments()
    from fit_plot import make_fit_plots
    make_fit_plots(agents,all_targets(targ_mode))#,moments_aux=moments_aux)
    
    
    mdl[0].mar_graphs()
    
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