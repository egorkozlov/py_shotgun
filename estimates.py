#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

from tiktak import filer
from calibration_params import calibration_params

def get_point(high_e,read_wisdom=False):

    if high_e:
        
        x = {'sigma_psi': 0.026922821456356934, 'sigma_psi_init': 0.06373681116518592, 'pmeet_21': 0.35569255597322275, 'pmeet_30': 1.4226673263521636, 'pmeet_40': 1.3684862955252783, 'preg_21': 0.008423036981683868, 'preg_28': 0.02771044116527361, 'preg_35': 0.025547961093756344, 'u_shift_kid': 0.3153482576154041, 'util_alp': 0.7482613736915342, 'util_kap': 6.917191837443176, 'util_qbar': 17.43025856523062, 'disutil_marry_sm_mal': 5.793717997829951, 'disutil_shotgun': 0.36650554551336484, 'abortion_costs': 2.960674281198284, 'p_abortion_access': 0.9827091238108862, 'u_lost_divorce': 0.9477696161601913, 'u_shift_couple': -0.16854069608836048,
             'high education': True}
        
        if read_wisdom:
            try:
                print('read wisdom from file!')
                o = filer('wisdom.pkl',0,0,repeat=False)[0]
                x = calibration_params()[-1](o[1])
                print(x)
                print('saved distance is {}'.format(o[0]))
                x.update({'high education': True})
            except:
                print('failed to read from wisdom file')
            
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.04071232742309963, 'sigma_psi_init': 0.0976736663373226, 'pmeet_21': 0.904865576861321, 'pmeet_30': 0.43036332267363914, 'pmeet_40': 1.2, 'preg_21': 0.03864821742604077, 'preg_28': 0.1401149063211062, 'preg_35': 0.08023116446696105, 'u_shift_kid': 0.37454369856900793, 'util_alp': 0.5583843839353859, 'util_kap': 6.13910894928199, 'util_qbar': 0.42904988385364523, 'disutil_marry_sm_mal': 2.795264868966648, 'disutil_shotgun': 0.2713254094887785, 'abortion_costs': 1.866347898322101, 'p_abortion_access': 0.8779040533392931, 'u_lost_divorce': 0.6174707512636085, 'u_shift_couple': -0.2758802590994992,
             'high education': False}
        
        targ_mode = 'low education'
        
        if read_wisdom:
            try:
                print('read wisdom from file!')
                o = filer('wisdom.pkl',0,0,repeat=False)[0]
                x = calibration_params()[-1](o[1])
                print(x)
                print('saved distance is {}'.format(o[0]))
                x.update({'high education': False})
            except:
                print('failed to read from wisdom file')
            
    
    return x, targ_mode
