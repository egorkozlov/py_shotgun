#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.1616829565907847, 'sigma_psi_init': 0.4925868917949559, 'pmeet_21': 0.19233181471138022, 'pmeet_30': 0.33949129640650677, 'pmeet_40': 0.7915036603428206, 'preg_21': 0.05071833739024931, 'preg_28': 0.03937540253500338, 'preg_35': 0.11636925499797958, 'u_shift_mar': 1.3565176613394831, 'util_alp': 0.010086848139632504, 'util_kap': 0.8217332522677137, 'util_qbar': 2.336393248618432, 'disutil_marry_sm_mal': 52.35330678635207, 'disutil_shotgun': 4.08263314583455, 'abortion_costs': 0.0, 'p_abortion_access': 0.7335225030510213, 'u_lost_divorce': 8.634602818759316, 'couple_rts': 0.4489943331461631,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.07262028069004121, 'sigma_psi_init': 0.5975384122500378, 'pmeet_21': 0.32747839443931864, 'pmeet_30': 0.040400190976707655, 'pmeet_40': 0.6472740842394364, 'preg_21': 0.1617772402053484, 'preg_28': 0.13779244582019745, 'preg_35': 0.11488115937628847, 'u_shift_mar': 2.0, 'util_alp': 0.01, 'util_kap': 0.7371085628065419, 'util_qbar': 1.6315860400437356, 'disutil_marry_sm_mal': 64.65098573395561, 'disutil_shotgun': 2.7831665371545045, 'abortion_costs': 0.0, 'p_abortion_access': 0.6288965038601756, 'u_lost_divorce': 0.0, 'couple_rts': 0.029883097408165055,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode