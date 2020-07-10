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
    x = {'sigma_psi': 0.10056612648397834,
         'sigma_psi_init': 0.26233143570897033,
         'pmeet_21': 0.1335781709881927,
         'pmeet_28': 0.36676421412123034,
         'pmeet_35': 0.6577934697353877,
         'preg_21': 0.03405237241832415,
         'preg_28': 0.055095768985394006,
         'preg_35': 0.04333002738196415,
         'u_shift_mar': 1.3054616179142995,
         'util_alp': 0.5055633884552083,
         'util_kap': 0.7587624279095282,
         'util_qbar': 0.0,
         'disutil_marry_sm_mal': 27.07777234023702,
         'disutil_shotgun': 5.2127711833143024,
         'abortion_costs': 0.0,
         'p_abortion_access': 1.0,
         'u_lost_divorce': 6.643996484294591}





    
    tar = target_values(targ_mode)
    
    
    this_name = 'real child support'
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      #load_from='mdl.pkl',
                                      verbose=True,draw=True,cs_moments=False,
                                      moments_save_name = this_name,
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
                   compare='college 0 0 child support.pkl',
                   base=this_name+'.pkl',
                   base_name=this_name,
                   compare_name='No child support',
                   #graphs_title_add="Experiment: Removing Subsistence Constraint",
                   moments_aux=None) #,moments_aux=moments_aux)
    
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
    
