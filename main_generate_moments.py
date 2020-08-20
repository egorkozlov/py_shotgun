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
    


def generate_counterfactuals(resume=True):
    for name, high_e in zip(['col','hs'],[True,False]):
        x, targ_mode = get_point(high_e)
        
        tar = target_values(targ_mode)
        
        
        adjustments = [('baseline',{}),
                       ('abortions and no stigma',{'disutil_shotgun':0.0,
                                                  'p_abortion_access':1.0,
                                                  'abortion_costs':0.0}),
                       ('no child support',{'child_support_awarded_nm':0.0,
                                            'child_support_awarded_div':0.0}),
                       ('full child support',{'child_support_awarded_nm':1.0,
                                            'child_support_awarded_div':1.0}),
                       ('in mar child support',{'child_support_awarded_nm':0.0,
                                            'child_support_awarded_div':1.0}),
                       ('out mar child support',{'child_support_awarded_nm':1.0,
                                            'child_support_awarded_div':0.0}),
                       ('infinite divorce costs',{'u_lost_divorce':200.0}),
                       ('no divorce costs',{'u_lost_divorce':0.0}),
                       ('full access to abortion',{'p_abortion_access':1.0}),
                       ('no abortions',{'p_abortion_access':0.0}),
                       ('costless abortion',{'p_abortion_access':1.0,
                                             'abortion_costs':0.0}),
                       ('no social stigma',{'disutil_shotgun':0.0}),
                       ('no remar penalty',{'disutil_marry_sm_mal':0.0}),
                       ('no remarriage',{'disutil_marry_sm_mal':500.0}),
                       ('no qbar',{'util_qbar':0.0}),
                       ('no unplanned pregnancy',{'preg_21': 0.0,
                                                  'preg_28': 0.0,
                                                  'preg_35': 0.0}),
                        ('no skills depreciation',{'z_drift':0.0}),
                        ('no pay gap',{'pay_gap':False}),
                        ('no kids',{'any kids':False}),
                        ('no partners',{'pmeet_21': 0.0,
                                        'pmeet_28': 0.0,
                                        'pmeet_35': 0.0}),
                        ('no home production',{'util_kappa':0.001})]
        
        for adj_name, fix in adjustments:
            
            print('\n\n\n\n\n\n')
            print('doing {} {}'.format(name,adj_name))
            
            x_new = x.copy()
            x_new.update(fix)
            
            try:
                
                fname = '{} {}.pkl'.format(name,adj_name)
                
                if resume:
                    try:
                        m = filer(fname,0,0,repeat=False)
                        del(m)
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
                                                  moments_save_name = '{} {}'.format(name,adj_name),
                                                  moments_repeat=5)
                    print("file {} saved".format(fname))
                else:
                    print("file {} already exists".format(fname))
                
           
                print('done, doing fit plots')
                
                try:                 
                    if adj_name == 'baseline':
                        fp = FitPlots(targ_mode=targ_mode,
                                       compare=None,
                                       base='{} baseline.pkl'.format(name),
                                       compare_name='Data',
                                       base_name='baseline',
                                       moments_aux=None)
                    else:
                        fp = FitPlots(targ_mode=targ_mode,
                                       compare='{} baseline.pkl'.format(name),
                                       base=fname,
                                       compare_name='baseline',
                                       base_name=adj_name,
                                       moments_aux=None)
                except:
                    print('something wrong with fit plots...')
                
            except KeyboardInterrupt:
                raise(KeyboardInterrupt)
            except BaseException as a:
                print("file {} {}.pkl could not be produced. Exception:".format(name,adj_name))
                print(a)
                print(' ')
    
    
if __name__ == '__main__':
    generate_counterfactuals()
