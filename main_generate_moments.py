#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:45:21 2020

@author: egorkozlov
"""

from estimates import get_point
from residuals import mdl_resid
from targets import target_values

for name, high_e in zip(['college','high school'],[True,False]):
    x, targ_mode = get_point(high_e)
    
    tar = target_values(targ_mode)
    
    adjustments = [('baseline',{}),
                   ('no social stigma',{'disutil_shotgun_coef':0.0}),
                   ('no remar penalry',{'disutil_marry_sm_mal_coef':0.0}),
                   ('no qbar',{'util_qbar':0.0}),
                   ('no unplanned pregnancy',{'preg_21': 0.0,
                                              'preg_28': 0.0,
                                              'preg_35': 0.0}),
                    ('2O percent child support',{'child_support_share':0.2}),
                    ('no skills depreciation',{'z_drift':0.0}),
                    ('no pay gap',{'pay_gap':False}),
                    ('no kids',{'any kids':False}),
                    ('no partners',{'pmeet_21': 0.0,
                                    'pmeet_28': 0.0,
                                    'pmeet_35': 0.0}),
                    ('no home production',{'util_kappa':0.001})]
    
    for adj_name, fix in adjustments:
        x_new = x.copy().update(fix)
        try:
            out,  mom = mdl_resid(x=x_new,targets=tar,
                                          return_format=['distance','moments'],
                                          verbose=False,draw=False,cs_moments=False,
                                          moments_save_name = '{} {}'.format(name,adj_name),
                                          moments_repeat=5)
            print("file {} {}.pkl saved".format(name,adj_name))
        except KeyboardInterrupt:
            raise(KeyboardInterrupt)
        except BaseException as a:
            print("file {} {}.pkl could not be produced. Exception:".format(name,adj_name))
            print(a)
            print(' ')
    
    