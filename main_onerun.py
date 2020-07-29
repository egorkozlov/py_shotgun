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
    
    x = {'sigma_psi': 0.17547979922272952,
         'sigma_psi_init': 0.2935764995685123,
         'pmeet_21': 0.1375957154357594,
         'pmeet_28': 0.3577852762275396,
         'pmeet_35': 0.7015588540881759,
         'preg_21': 0.044896015921624136,
         'preg_28': 0.05190356560351397,
         'preg_35': 0.039391507653039434,
         'u_shift_mar': 1.0763146779858002,
         'util_alp': 0.10000192022650033,
         'util_kap': 0.9050768407157648,
         'util_qbar': 1.8065820813416464,
         'disutil_marry_sm_mal': 40.54954577056158,
         'disutil_shotgun': 7.734445729198316,
         'abortion_costs': 0.7995715430427783,
         'p_abortion_access': 0.9970209421052573,
         'u_lost_divorce': 10.00156200327671}






    
    tar = target_values(targ_mode)
    
    
    this_name = 'baseline'
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
                   compare=None,
                   base=this_name+'.pkl',
                   compare_name='data',
                   base_name='baseline',#this_name,
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
    
