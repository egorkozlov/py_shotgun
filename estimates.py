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
        
        x = {'sigma_psi': 0.02271227468311225, 'sigma_psi_init': 0.058638977435951305, 'pmeet_21': 0.3775858999749023, 'pmeet_30': 1.1114312962137847, 'pmeet_40': 0.31014313380175096, 'preg_21': 0.007431431021009238, 'preg_28': 0.02892174895848225, 'preg_35': 0.024209123328699594, 'u_shift_kid': 0.31771592374588375, 'util_alp': 0.771834843542828, 'util_kap': 6.902234564557551, 'util_qbar': 9.473335136133066, 'disutil_marry_sm_mal': 3.2472913078044545, 'disutil_shotgun': 0.4890157219311151, 'abortion_costs': 1.5886358570521772, 'p_abortion_access': 0.8174304666969183, 'u_lost_divorce': 0.9539818957506494, 'u_shift_couple': -0.15110699130997146,
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
        x = {'sigma_psi': 0.1, 'sigma_psi_init': 0.06116524742294942, 'pmeet_21': 0.8113766193192298, 'pmeet_30': 0.5810854411953447, 'pmeet_40': 1.2888197016633065, 'preg_21': 0.10931640666793829, 'preg_28': 0.05149518424513294, 'preg_35': 0.10182921767475042, 'u_shift_kid': 0.4031345423311973, 'util_alp': 0.6494460777054719, 'util_kap': 6.043394416570878, 'util_qbar': 0.0, 'disutil_marry_sm_mal': 5.768217858740022, 'disutil_shotgun': 1.3338166895953723, 'abortion_costs': 0.0, 'p_abortion_access': 0.7538834914491163, 'u_lost_divorce': 2.008607758320293, 'u_shift_couple': -0.3431072199441293,
             'high education': False}
        
        targ_mode = 'low education'
        
        if read_wisdom:
            try:
                print('read wisdom from file!')
                o = filer('/projects/p31069/py_shotgun/wisdom.pkl',0,0,repeat=False)[0]
                x = calibration_params()[-1](o[1])
                print(x)
                print('saved distance is {}'.format(o[0]))
                x.update({'high education': False})
            except:
                print('failed to read from wisdom file')
            
    
    return x, targ_mode
