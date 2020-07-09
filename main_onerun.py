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
    x = {'sigma_psi': 0.8903998377617001,
         'sigma_psi_mult': 5.008454375933519,
         'pmeet_21': 0.26461668508226976,
         'pmeet_28': 0.3808517813188633,
         'pmeet_35': 0.3537668845748406,
         'preg_21': 0.03301912411330818,
         'preg_28': 0.04337421298029398,
         'preg_35': 0.06499192886306442,
         'u_shift_mar': 1.4329636375369272,
         'util_alp': 0.4980641699631754,
         'util_kap': 0.8967664551452008,
         'util_qbar': 1.0457400305750946,
         'disutil_marry_sm_mal_coef': 14.444200806643657,
         'disutil_shotgun_coef': 2.117334990378265,
         'abortion_costs_mult': 8.430713206729187,
         'p_abortion_access': 1.0,
         'u_lost_divorce_mult': 2.690742151483817}




    
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
    
