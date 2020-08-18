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



    x = {'sigma_psi': 0.25,'sigma_psi_init': 2.332155241562433,'pmeet_21': 0.30325973411684065,'pmeet_28': 0.7329390078838265,'pmeet_35': 1.0,'preg_21': 0.017160684393485684,'preg_28': 0.03733205098574577,'preg_35': -0.0014982571526763744,'u_shift_mar': 1.1123578939396237,'util_alp': 0.3061592839066737,'util_kap': 0.8512902458486851,'util_qbar': 1.1171096340029305,'disutil_marry_sm_mal': 64.12382814306783,'disutil_shotgun': 5.78268701428576,'abortion_costs': 30.926325607734366,'p_abortion_access': 0.9916006844636375,'u_lost_divorce': 6.555383565811324}




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
