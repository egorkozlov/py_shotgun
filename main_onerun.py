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
    
    high_e = False
    x, targ_mode = get_point(high_e)



    x = {'sigma_psi': 0.39953601382428866, 'sigma_psi_init': 3.0725683838557276, 'pmeet_21': 0.432991116061913, 'pmeet_28': 0.26180839230577446, 'pmeet_35': 0.5691992566974007, 'preg_21': 0.170975611759682, 'preg_28': 0.21779042931083295, 'preg_35': 0.17575004499028102, 'u_shift_mar': 1.4998917495906658, 'util_alp': 0.4179831829995812, 'util_kap': 1.3061403769026076, 'util_qbar': 0.003578373701123387, 'disutil_marry_sm_mal': 51.48144648331411, 'disutil_shotgun': 0.15849135146243648, 'abortion_costs': 40.05274285709042, 'p_abortion_access': 0.9759747406920513, 'u_lost_divorce': 0.5704761664715526}



    tar = target_values(targ_mode)


    this_name = 'no college baseline'
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
