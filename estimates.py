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
        
        x = {'sigma_psi': 0.02297391391828134, 'sigma_psi_init': 0.061389470944663196, 'pmeet_21': 0.35427476693487037, 'pmeet_30': 1.101469983799187, 'pmeet_40': 0.0, 'preg_21': 0.0021972818487067947, 'preg_28': 0.03194765436380882, 'preg_35': 0.02558129855186625, 'u_shift_kid': 0.34592995681350347, 'util_alp': 0.8554745878335215, 'util_kap': 6.749928240029728, 'util_qbar': 16.614101270795228, 'disutil_marry_sm_mal': 3.262093358061867, 'disutil_shotgun': 0.5084201219022069, 'abortion_costs': 2.251099749734383, 'p_abortion_access': 0.8949286724553163, 'u_lost_divorce': 1.0648120485130006, 'u_shift_couple': -0.16004255555405467,
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
        x = {'sigma_psi': 0.04119975516565719, 'sigma_psi_init': 0.07184509981781, 'pmeet_21': 0.7300641341551373, 'pmeet_30': 0.38552526708748397, 'pmeet_40': 1.4132304041226518, 'preg_21': 0.1029100967053943, 'preg_28': 0.11241132276639117, 'preg_35': 0.11203564468462099, 'u_shift_mar': 0.338428482678413, 'util_alp': 0.5195282434982275, 'util_kap': 7.152398760885778, 'util_qbar': 0.0, 'disutil_marry_sm_mal': 3.18966037249299, 'disutil_shotgun': 0.3647670950676456, 'abortion_costs': 0.2962878054482049, 'p_abortion_access': 0.6662167114665236, 'u_lost_divorce': 0.5275074834332285, 'mu_psi_init': -0.24342175587968384,
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
