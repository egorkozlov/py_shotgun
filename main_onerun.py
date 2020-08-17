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



    x = {'sigma_psi': 0.24000883142940208,
        'sigma_psi_init': 2.270662128160131,
        'pmeet_21': 0.282197067097691,
        'pmeet_28': 0.717301599527013,
        'pmeet_35': 0.9998565249924385,
        'preg_21': 0.027242779157600316,
        'preg_28': 0.034354384984791576,
        'preg_35': 0.024562299247326433,
        'u_shift_mar': 0.8600153678932919,
        'util_alp': 0.1, 'util_kap': 0.9399627350061319,
        'util_qbar': 1.2755593891876527,
        'disutil_marry_sm_mal': 69.13458440788543,
        'disutil_shotgun': 4.870682077608209,
        'abortion_costs': 30.308034843501463,
        'p_abortion_access': 0.8956773350793297,
        'u_lost_divorce': 5.492799438755628}




    tar = target_values(targ_mode)


    this_name = 'baseline'
    out, mdl, agents, res, mom = mdl_resid(x=x,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals','moments'],
                                      save_to='mdl.pkl',
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
                   compare='baseline.pkl', #None,
                   base='no divorce costs.pkl',#this_name+'.pkl',
                   compare_name='baseline', #'Data',
                   base_name='no divorce costs',
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
