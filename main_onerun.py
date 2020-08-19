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


    x = {'sigma_psi': 0.3606833154484394, 'sigma_psi_init': 3.031242467682322, 'pmeet_21': 0.2981951262856033, 'pmeet_28': 0.7468405296611746, 'pmeet_35': 0.9700424253585994, 'preg_21': 0.023056886381180183, 'preg_28': 0.03585409029595745, 'preg_35': 0.013599676361164492, 'u_shift_mar': 1.1114058257254988, 'util_alp': 0.2509204297743008, 'util_kap': 0.8802848731701842, 'util_qbar': 1.157339328421708, 'disutil_marry_sm_mal': 76.65446606843315, 'disutil_shotgun': 7.079428268293057, 'abortion_costs': 33.74015195752216, 'p_abortion_access': 0.9974677343975493, 'u_lost_divorce': 8.105777902557344}



    tar = target_values(targ_mode)


    this_name = 'baseline'
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      #save_to='mdl.pkl',
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
                   compare_name='Data',
                   base_name='baseline',
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
    return locals()

if __name__ == '__main__':
    allthings = main()
    globals().update(**allthings)
