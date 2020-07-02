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
from estimates import get_point
 
if __name__ == '__main__':
    
    high_e = True
    x, targ_mode = get_point(high_e)
    x = {'sigma_psi': 0.1342648840919196,
         'sigma_psi_mult': 3.9213102635901445,
         'pmeet_21': 0.19462609220769422,
         'pmeet_28': 0.46904254700807024,
         'pmeet_35': 0.4853799346774549,
         'preg_21': 0.07039396372104786,
         'preg_28': 0.03603579972872427,
         'preg_35': 0.06930057037990973,
         'u_shift_mar': 1.31679256124115,
         'util_alp': 0.4172232117671171,
         'util_kap': 0.8335867276157953,
         'util_qbar': 0.909496552645544,
         'disutil_marry_sm_mal_coef': 28.469082232141393,
         'disutil_shotgun_coef': 0.10979562467284859,
         'abortion_costs_mult': 0.6355790686264702,
         'p_abortion_access': 0.5485976858824111,
         'u_lost_divorce_mult': 4.157201180190372}





    
    tar = target_values(targ_mode)
    
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      #load_from='mdl.pkl',
                                      verbose=True,draw=True,cs_moments=False,
                                      moments_save_name = 'college move social stigma',
                                      moments_repeat=2)
    
    mdl[0].time_statistics()
                         
    

    
    print('Done. Residual in point x0 is {}'.format(out))
    
    from targets import all_targets
    tar_fk = all_targets('low education')
    mom_fk = {key: tar_fk[key][0] for key in tar_fk}
    
    
    #from simulations import Agents
    #moments_aux = Agents( mdl, N=10000, T=18, female=False, verbose=False).aux_moments()
    
    from fit_plot import FitPlots
    fp = FitPlots(targ_mode=targ_mode,
                   base='college move social stigma.pkl',
                   compare=None,
                   base_name='Model',
                   compare_name='Data',
                   #graphs_title_add="Experiment: Removing Subsistence Constraint",
                   moments_aux=None)#,moments_aux=moments_aux)
    
    '''
    from fit_plot import FitPlots
    fp = FitPlots(targ_mode=targ_mode,
                   base='college no taste shocks.pkl',
                   compare='college sigma zero.pkl',
                   base_name='no shocks',
                   compare_name='sigma zero',
                   #graphs_title_add="Experiment: Removing Subsistence Constraint",
                   moments_aux=None)#,moments_aux=moments_aux)
    
    mdl[0].mar_graphs()
    '''
    
