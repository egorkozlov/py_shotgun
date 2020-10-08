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
        
        x = {'sigma_psi': 0.024447077354810887, 'sigma_psi_init': 0.026287447678799637, 'pmeet_21': 0.16181842291258391, 'pmeet_30': 0.4844756482942769, 'pmeet_40': 0.6278134903824727, 'preg_21': 0.008825649168577587, 'preg_28': 0.029699717484135257, 'preg_35': 0.01697979902262614, 'u_shift_kid': 0.29959929396641477, 'util_alp': 0.7103357592468402, 'util_kap': 7.232876352397827, 'util_qbar': 40.0, 'disutil_marry_sm_mal': 2.109356551078698, 'disutil_shotgun': 0.812424825412403, 'abortion_costs': 1.0847666887625578, 'p_abortion_access': 1.0, 'u_lost_divorce': 0.930254920866958, 'u_shift_couple': -0.055377440502576314,
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
        x = {'sigma_psi': 0.04660016232690957, 'sigma_psi_init': 0.06011470653209264, 'pmeet_21': 0.4870963730334307, 'pmeet_30': 0.22882085040071745, 'pmeet_40': 1.2755288412134194, 'preg_21': -0.09067995907595618, 'preg_28': 0.09795613506755489, 'preg_35': -0.01955824358876551, 'u_shift_kid': 0.32832176081503095, 'util_alp': 0.5231414770939987, 'util_kap': 5.489896027104729, 'util_qbar': 1.3521599036304641, 'disutil_marry_sm_mal': 1.3820898067762082, 'disutil_shotgun': 0.0, 'abortion_costs': 0.0, 'p_abortion_access': 0.6888998913368529, 'u_lost_divorce': 0.6717570527862875, 'u_shift_couple': -0.16838304061498055,
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
