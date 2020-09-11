#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27

@author: Egor Kozlov
"""


#import warnings
#warnings.filterwarnings("error")


from platform import system
import os

if system() != 'Darwin' and system() != 'Windows':
    os.environ['QT_QPA_PLATFORM']='offscreen'



from residuals import mdl_resid
from targets import target_values

print('Hi!')


os.environ['MKL_CBWR']='AUTO'

from estimates import get_point
 
def main(read_wisdom=False,erase=False):
    
    high_e = True
    
    if erase:
        try:
            os.remove('az_dist_fem.pkl')
            print('removed')
        except:
            pass

        try:
            os.remove('az_dist_mal.pkl')
            print('removed')
        except:
            pass
    
    x, targ_mode = get_point(high_e,read_wisdom=read_wisdom)
    
        
    
    tar = target_values(targ_mode)


    this_name = 'college baseline'
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
                   base_name='high school baseline',
                   moments_aux=None) #,moments_aux=moments_aux)
    
    mdl[0].mar_graphs()
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
