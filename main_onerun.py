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
    
    
    
    x = {'sigma_psi': 0.1853324578327308,
         'sigma_psi_init': 0.9399973974581646,
         'pmeet_21': 0.23057044853525627, 
         'pmeet_28': 0.6302582022794753,
         'pmeet_35': 1.0,
         'preg_21': 0.004632574153187789,
         'preg_28': 0.030915854690938757,
         'preg_35': -0.03835926491296945,
         'u_shift_mar': 1.1634888476823635,
         'util_alp': 0.21921319884901105,
         'util_kap': 0.8304799413261448,
         'util_qbar': 1.5124778292649947,
         'disutil_marry_sm_mal': 48.78355862063415,
         'disutil_shotgun': 4.611347004670638,
         'abortion_costs': 39.894265116705,
         'p_abortion_access': 0.990470993900194,
         'u_lost_divorce': 8.788992431309284}



    
    tar = target_values(targ_mode)
    
    
    this_name = 'full child support'
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
                   compare='baseline.pkl',
                   base=this_name+'.pkl',
                   compare_name='baseline',
                   base_name=this_name,
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
    
