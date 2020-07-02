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
    x = {'sigma_psi': 0.32141972647150646,
         'sigma_psi_mult': 4.89470848401967,
         'pmeet_21': 0.15003551504477958,
         'pmeet_28': 0.364733559917808,
         'pmeet_35': 0.45932227232680173,
         'preg_21': 0.11705221172526127,
         'preg_28': 0.06031520195756093,
         'preg_35': 0.13847403463524244,
         'u_shift_mar': 1.4609834826218564,
         'util_alp': 0.3983042535166076,
         'util_kap': 0.8601277198595394,
         'util_qbar': 0.9201475740242766,
         'disutil_marry_sm_mal_coef': 17.10336823543867,
         'disutil_shotgun_coef': 0.5285761051437314,
         'abortion_costs_mult': 20.16630464635698,
         'p_abortion_access': 0.5626163254818022,
         'u_lost_divorce_mult':2.0}




    
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
                   compare='college no taste shocks.pkl',
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
    
