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
    
    
    
    x = {'sigma_psi': 0.23571313928571477,
         'sigma_psi_init': 2.123636664878597,
         'pmeet_21': 0.2873647689783723,
         'pmeet_28': 0.7553774259051128,
         'pmeet_35': 1.0,
         'preg_21': 0.017461033039670462,
         'preg_28': 0.024901095174551184,
         'preg_35': -0.015475818457370151, 
         'u_shift_mar': 1.2306521824428305,
         'util_alp': 0.3854997739795115,
         'util_kap': 0.8103116344047732,
         'util_qbar': 1.0505842569863009,
         'disutil_marry_sm_mal': 60.0,
         'disutil_shotgun': 5.581480022921534,
         'abortion_costs': 30.940624456940288,
         'p_abortion_access': 0.1799031427804904,
         'u_lost_divorce': 6.673691715831503}



    
    tar = target_values(targ_mode)
    
    
    this_name = 'out mar child support'
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
    
