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
        
        x = {'sigma_psi': 0.02959122856845663, 'sigma_psi_init': 0.0645817678385954, 'pmeet_21': 0.4003805778323115, 'pmeet_30': 1.2967939688421732, 'pmeet_40': 1.1772265098112105, 'preg_21': 0.01210926599130474, 'preg_28': 0.026879685039007825, 'preg_35': 0.028107396649081767, 'u_shift_kid': 0.32950464792871087, 'util_alp': 0.7920911344385484, 'util_kap': 6.976799354216199, 'util_qbar': 19.664047254785043, 'disutil_marry_sm_mal': 8.012136894740857, 'disutil_shotgun': 0.4049001885442925, 'abortion_costs': 3.22160688092188, 'p_abortion_access': 1.0, 'u_lost_divorce': 1.138783218504859, 'u_shift_couple': -0.1748334763349187,
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
        x = {'sigma_psi': 0.02481015097391902, 'sigma_psi_init': 0.08871073033936372, 'pmeet_21': 1.1287440759318557, 'pmeet_30': 0.4661001967264672, 'pmeet_40': 1.2, 'preg_21': 0.21599765897582557, 'preg_28': 0.1319052282520419, 'preg_35': 0.3051665590352528, 'u_shift_kid': 0.3124032931085869, 'util_alp': 0.6832818585750131, 'util_kap': 8.324765874087191, 'util_qbar': 0.0, 'disutil_marry_sm_mal': 1.9080112991697709, 'disutil_shotgun': 0.0, 'abortion_costs': 0.0, 'p_abortion_access': 0.5025036452746912, 'u_lost_divorce': 0.21258980489241242, 'u_shift_couple': -0.23626889664541162,
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
