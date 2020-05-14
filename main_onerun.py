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
    x = {'sigma_psi': 0.10766216584844603,
         'sigma_psi_mult': 6.093923962547796,
         'pmeet_21': 0.19813694064318865,
         'pmeet_28': 0.3537652569134399,
         'pmeet_35': 0.49174708026674585,
         'preg_21': 0.03209762599843419,
         'preg_28': 0.06586923576148838,
         'preg_35': 0.09454363441239991,
         'u_shift_mar': 1.1177968780093115,
         'util_alp': 0.2825390734412878,
         'util_kap': 0.9707843516044867,
         'util_qbar': 0.8576029389166274,
         'disutil_marry_sm_mal_coef': 15.490767514548592,
         'disutil_shotgun_coef': 0.31299196911695903,
         'taste_shock_mult': 8.0}
    '''
    '''
    x = {'sigma_psi': 0.154819596272499,
         'sigma_psi_mult': 5.369369910975616,
         'pmeet_21': 0.1947578521767194,
         'pmeet_28': 0.35317172001916,
         'pmeet_35': 0.4462870780336654,
         'preg_21': 0.03348310308060995,
         'preg_28': 0.06187733817261159,
         'preg_35': 0.015229569104669022,
         'u_shift_mar': 1.3651711748701434,
         'util_alp': 0.3796474450996089,
         'util_kap': 0.9134203294944718,
         'util_qbar': 0.8586883551899729,
         'disutil_marry_sm_mal_coef': 11.807589876890876,
         'disutil_shotgun_coef': 0.3178880464362836,
         'taste_shock_mult': 7.999210546887254}
    '''
    '''
    x = {'sigma_psi': 0.12236224527121549,
         'sigma_psi_mult': 5.704201308752313,
         'pmeet_21': 0.20423533383772224,
         'pmeet_28': 0.35543855140338027,
         'pmeet_35': 0.4672960637883365,
         'preg_21': 0.02543432257425196,
         'preg_28': 0.06978730091192242,
         'preg_35': 0.04652343564738676,
         'u_shift_mar': 1.1488025339016394,
         'util_alp': 0.29512469844357203,
         'util_kap': 0.983901025229017,
         'util_qbar': 0.8529706428214711,
         'disutil_marry_sm_mal_coef': 13.082089056032371,
         'disutil_shotgun_coef': 0.3153246960829357,
         'taste_shock_mult': 8.0}
    '''
        
    x = {'sigma_psi': 0.11505749117016362,
         'sigma_psi_mult': 5.852340399857861,
         'pmeet_21': 0.24984766129215077,
         'pmeet_28': 0.31577591457938914,
         'pmeet_35': 0.47673262106126735,
         'preg_21': 0.015410747143281972,
         'preg_28': 0.06345376608392728,
         'preg_35': -0.04610072125697056,
         'u_shift_mar': 1.408902651023607,
         'util_alp': 0.35196894037684356,
         'util_kap': 0.8544636657426437,
         'util_qbar': 0.9007017830322368,
         'disutil_marry_sm_mal_coef': 15.192983709511992,
         'disutil_shotgun_coef': 0.25682426989088775,
         'taste_shock_mult': 14.015173971372116}




    

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