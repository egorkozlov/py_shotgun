#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.1477726124226512, 'sigma_psi_init': 0.12383104968331363, 'pmeet_21': 0.10747866149044277, 'pmeet_30': 0.21280921124682906, 'pmeet_40': 0.34735213888811073, 'preg_21': 0.01483992873884013, 'preg_28': 0.043396851654720174, 'preg_35': 0.0024733399057522906, 'u_shift_mar': 1.4368759259133626, 'util_alp': 0.3002084657313606, 'util_kap': 0.8288657954190843, 'util_qbar': 1.9303252338152943, 'disutil_marry_sm_mal': 46.66344149694794, 'disutil_shotgun': 4.361021788613633, 'abortion_costs': 62.28228453253508, 'p_abortion_access': 0.9223803654650251, 'u_lost_divorce': 6.5404973607586685,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.1995078090554912, 'sigma_psi_init': 0.9, 'pmeet_21': 0.34295958151401773, 'pmeet_30': 0.024091666504720194, 'pmeet_40': 0.8551404126131842, 'preg_21': 0.13755613313369197, 'preg_28': 0.11912197905077387, 'preg_35': -0.04998999555370965, 'u_shift_mar': 2.0, 'util_alp': 0.08890551692982442, 'util_kap': 0.7419226543585731, 'util_qbar': 1.3893245736440478, 'disutil_marry_sm_mal': 70.95135376615653, 'disutil_shotgun': 3.5184521159271522, 'abortion_costs': 3.004395052551469, 'p_abortion_access': 0.00026736771206510443, 'u_lost_divorce': 0.0370961835703107,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode