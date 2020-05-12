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
    x = {'sigma_psi': 0.5814121121665838,
         'sigma_psi_mult': 5.166848399191245,
         'pmeet_21': 0.16381625559794244,
         'pmeet_28': 0.3691365690289941,
         'pmeet_35': 0.3879313316537848,
         'preg_21': 0.02794620197864167,
         'preg_28': 0.030777480850298217,
         'preg_35': 0.009246598699316738,
         'u_shift_mar': 1.8014493130394245,
         'util_alp': 0.6057866999755354,
         'util_kap': 0.8433330370075438,
         'disutil_marry_sm_mal_coef': 11.010141916862054,
         'disutil_shotgun_coef': 0.8316044174306312,
         'taste_shock_mult':1.0}
    '''
    '''
    x = {'sigma_psi': 0.4902891581150292,
         'sigma_psi_mult': 5.702606626187299,
         'pmeet_21': 0.15015478921321068,
         'pmeet_28': 0.3710088343850946,
         'pmeet_35': 0.44012274619061426,
         'preg_21': 0.10903649039754515,
         'preg_28': 0.030979962029004545,
         'preg_35': 0.029828705099284203,
         'u_shift_mar': 1.2,
         'util_alp': 0.1462819544897464,
         'util_kap': 1.100359301455061,
         'util_qbar': 0.9546130214227615,
         'disutil_marry_sm_mal_coef': 8.404348024770918,
         'disutil_shotgun_coef': 0.8347214347081617,
         'taste_shock_mult': 4.930701825497131}
    '''
    x = {'sigma_psi': 0.26138747738791257,
         'sigma_psi_mult': 5.179254639435019,
         'pmeet_21': 0.17986335029641237,
         'pmeet_28': 0.32408751651872386,
         'pmeet_35': 0.36426615699680626,
         'preg_21': 0.06251572908756534,
         'preg_28': 0.031170308112898497,
         'preg_35': 0.10674870611243295,
         'u_shift_mar': 1.2732222822605368,
         'util_alp': 0.35015283870308955,
         'util_kap': 0.9462115202985493,
         'util_qbar': 0.8449866159744631,
         'disutil_marry_sm_mal_coef': 8.321891582811833,
         'disutil_shotgun_coef': 0.23939689770978267,
         'taste_shock_mult': 4.7870173423387525}


    

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