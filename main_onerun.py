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
 
def main():
    
    high_e = True
    x, targ_mode = get_point(high_e)
    
    x = {'sigma_psi': 0.06883045907429694,
         'sigma_psi_init': 0.08943321614503934,
         'pmeet_21': 0.14245350436785614,
         'pmeet_28': 0.35351576474635815,
         'pmeet_35': 0.5803847682254825,
         'preg_21': 0.020712755871815328,
         'preg_28': 0.06153082169980899,
         'preg_35': 0.03295250149252342,
         'u_shift_mar': 1.3469652448502023,
         'util_alp': 0.5643002665831136,
         'util_kap': 0.7267898392026289,
         'util_qbar': 0.11114171525939373,
         'disutil_marry_sm_mal': 18.850868358463234,
         'disutil_shotgun': 0.011217472786367675,
         'abortion_costs': 4.138337665129658,
         'p_abortion_access': 0.9992984327036426,
         'u_lost_divorce': 5.912838282187521}






    
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
    '''
    from fit_plot import FitPlots
    fp = FitPlots(targ_mode=targ_mode,
                   compare='baseline.pkl',
                   base=this_name+'.pkl',
                   compare_name='baseline',
                   base_name=this_name,
                   #graphs_title_add="Experiment: Removing Subsistence Constraint",
                   moments_aux=None) #,moments_aux=moments_aux)
    '''
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
    
if __name__ == '__main__':
    main()