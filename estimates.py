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
        
        x = {'sigma_psi': 0.025092055726001986, 'sigma_psi_init': 0.057524117123893596, 'pmeet_21': 0.37414286769758254, 'pmeet_30': 1.1932756253573578, 'pmeet_40': 1.2, 'preg_21': 0.016846303459086942, 'preg_28': 0.026031363571378165, 'preg_35': 0.028147668865726228, 'u_shift_kid': 0.25811605249933434, 'util_alp': 0.6202419382122233, 'util_kap': 8.328279006156304, 'util_qbar': 19.58211188499201, 'disutil_marry_sm_mal': 11.036758925366986, 'disutil_shotgun': 0.49670652381536956, 'abortion_costs': 2.2609160794229095, 'p_abortion_access': 0.9891922518917426, 'u_lost_divorce': 1.0689656664994374, 'u_shift_couple': -0.16143081649350327,
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
        x = {'sigma_psi': 0.05016131011798489, 'sigma_psi_init': 0.06012613381321384, 'pmeet_21': 0.6668575809221471, 'pmeet_30': 0.3370659273435091, 'pmeet_40': 1.1932072216686185, 'preg_21': 0.14464471221184672, 'preg_28': 0.08409002905918418, 'preg_35': 0.09758353505187373, 'u_shift_kid': 0.3244994388255927, 'util_alp': 0.4871696303910731, 'util_kap': 7.407624026588406, 'util_qbar': 0.0, 'disutil_marry_sm_mal': 3.573678892532197, 'disutil_shotgun': 0.6086735804092326, 'abortion_costs': 0.0, 'p_abortion_access': 0.693696849009367, 'u_lost_divorce': 0.8217865600943532, 'u_shift_couple': -0.24510193475666364,
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
