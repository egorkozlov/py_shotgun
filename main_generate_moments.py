#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:45:21 2020

@author: egorkozlov
"""

from estimates import get_point
from residuals import mdl_resid
from targets import target_values
from fit_plot import FitPlots
from tiktak import filer
   
    
# this is basic function to solve a model with possible fix to parameters

def run(adj_name,fix,educ_name,resume=False,noplot=False):
    # at first this takes point x saved in estimates.py 
    # by calling get_point(high_e). x is a dict with parameter names and values
    # then it applies adjustment fix (a dict) to the point x
    # resulting file name is educ_name + adj_name
    # adj_name -- name of adjustment
    # fix -- dictionary containing parameters of setup.py to change
    # educ_name -- name of education group ("col" or "hs")
    # high_e -- input to get_point function, True for college, False for HS
    
    high_e = (educ_name == 'col' or educ_name == 'college')
    low_e = (educ_name == 'hs' or educ_name == 'high school')
    
    assert (high_e or low_e), 'wrong specifier for education'
    
    x, targ_mode = get_point(high_e,read_wisdom=False)

    tar = target_values(targ_mode)

    print('\n\n\n\n\n\n')
    print('doing {} {}'.format(educ_name,adj_name))

    x_new = x.copy()
    
    fix = fix.copy()
    
    if 'multiply' in fix:
        mult = fix.pop('multiply')
        for m in mult: x_new[m] *= mult[m]
        
    x_new.update(fix)
    
    print(x_new)
    
    name = '{} {}'.format(educ_name,adj_name)
    fname = '{}.pkl'.format(name)

    try:


        if resume:
            try:
                mom = filer(fname,0,0,repeat=False)
                skip = True
            except:
                skip = False
        else:
            skip = False

        if not skip:
            print("computing {}".format(fname))
            out,  mom = mdl_resid(x=x_new,targets=tar,
                                          return_format=['distance','moments'],
                                          verbose=False,draw=False,cs_moments=False,
                                          save_to ='mdl for {}'.format(fname),
                                          moments_save_name = name,
                                          moments_repeat=5,Tsim=42)
            print("file {} saved".format(fname))
        else:
            print("file {} already exists".format(fname))

        
        print('done, doing fit plots')

        try:                 
            if adj_name == 'baseline':
                fp = FitPlots(targ_mode=targ_mode,
                               compare=None,
                               base='{} baseline.pkl'.format(educ_name),
                               compare_name='Data',
                               base_name='baseline',
                               moments_aux=None,
                               noplot=noplot)
            else:
                fp = FitPlots(targ_mode=targ_mode,
                               compare='{} baseline.pkl'.format(educ_name),
                               base=fname,
                               compare_name='baseline',
                               base_name=adj_name,
                               moments_aux=None,
                               noplot=noplot)
        except:
            print('something wrong with fit plots...')

    except KeyboardInterrupt:
        raise(KeyboardInterrupt)
    except BaseException as a:
        print("file {} {}.pkl could not be produced. Exception:".format(name,adj_name))
        print(a)
        print(' ')
    return mom, fp

def adj_list(return_dict=False):
    # this is a list of possible countefactual scenarios with names on them
    
    
    adjustments = [('baseline',{}),
                   ('double social stigma',{'multiply':{'disutil_shotgun':2.0}}),
                   #('plus ten percent social stigma',{'multiply':{'disutil_shotgun':1.1}}),
                   ('infinite social stigma',{'disutil_shotgun':400.0}),
                   ('abortions and no stigma',{'disutil_shotgun':0.0,
                                              'p_abortion_access':1.0,
                                              'abortion_costs':0.0}),
                   ('no social stigma',{'disutil_shotgun':0.0}),
                   ('no child support',{'child_support_awarded_nm':0.0,
                                        'child_support_awarded_div':0.0}),
                   ('full child support',{'child_support_awarded_nm':1.0,
                                        'child_support_awarded_div':1.0}),
                   ('in mar child support',{'child_support_awarded_nm':0.0,
                                        'child_support_awarded_div':1.0}),                   
                   ('double divorce costs',{'multiply':{'u_lost_divorce':2.0}}),
                   ('half divorce costs',{'multiply':{'u_lost_divorce':0.5}}),
                   ('no unplanned pregnancy unanticipated',{'ppreg_sim_mult':0.0}),
                   #('no taxes to couples',{'tax_childless_couples':False,
                   #                        'tax_couples_woth_children':False}),
                   #('no taxes to couples with children',{'tax_couples_woth_children':False}),
                   #('no taxes to single mothers',{'tax_single_mothers':False}),
                   #('no taxes to ones with children',{'tax_single_mothers':False,
                   #                                   'tax_couples_woth_children':False}),
                   #('out mar child support',{'child_support_awarded_nm':1.0,
                   #                     'child_support_awarded_div':0.0}),
                   ('infinite divorce costs',{'u_lost_divorce':200.0}),
                   ('no divorce costs',{'u_lost_divorce':0.0}),
                   #('full access to abortion',{'p_abortion_access':1.0}),
                   ('no abortions',{'p_abortion_access':0.0}),
                   ('costless abortion',{'p_abortion_access':1.0,
                                         'abortion_costs':0.0}),
                   ('no remar penalty',{'disutil_marry_sm_mal':0.0}),
                   ('infinite remar penalty',{'disutil_marry_sm_mal':500.0}),
                   ('no subsistence constraint',{'util_qbar':0.0}),
                   ('no unplanned pregnancy',{'preg_21': 0.0,
                                              'preg_28': 0.0,
                                              'preg_35': 0.0}),
                   ('no skills depreciation',{'z_drift':0.0}),
                   ('no pay gap',{'pay_gap':False}),
                   ('no kids',{'any kids':False}),
                    #('no partners',{'pmeet_21': 0.0,
                    #                'pmeet_28': 0.0,
                    #                'pmeet_35': 0.0}),
                   ('no home production',{'util_kap':0.001})]
    
    
    
    
    
    
    if return_dict: adjustments = dict(adjustments)
    return adjustments
    
    
def generate_counterfactuals(resume=True):
    # this generates many instances of counterfactuals
    
    
    adjustments = adj_list()
    
    # produce baseline
    
    
    '''
    # runs baseline first
    run(*adjustments[0],'col')
    run(*adjustments[0],'hs')
    '''
    
    for adj_name, fix in adjustments:    
        for educ_name in ['col','hs']:
            run(adj_name,fix,educ_name,resume=resume)
    
    
if __name__ == '__main__':
    generate_counterfactuals()
