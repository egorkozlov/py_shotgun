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
        
        x = {'sigma_psi': 0.6656407364967266, 'sigma_psi_init': 1.2788831549969952, 'pmeet_21': 0.32769455723380525, 'pmeet_30': 0.858815496396281, 'pmeet_40': 0.49584447542713284, 'preg_21': 0.015656403414957637, 'preg_28': 0.03488107732234185, 'preg_35': 0.03338858066901612, 'u_shift_mar': 0.3035027068097137, 'util_alp': 0.4451218545526679, 'util_kap': 11.281001265137967, 'util_qbar': 34.38536826299227, 'disutil_marry_sm_mal': 20.923836572884063, 'disutil_shotgun': 5.336872131318279, 'abortion_costs': 6.265297679220528, 'p_abortion_access': 0.8640897302237546, 'u_lost_divorce': 14.947099697801866, 'mu_psi_init': -1.4448913321100947,
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
        x = {'pmeet_21': 0.8098504080749165, 'pmeet_30': 0.29793588617668143, 'pmeet_40': 0.9999510121934609, 'preg_21': 0.1450434235061589, 'preg_28': 0.057965809812286634, 'preg_35': 0.5041573799847998, 'util_qbar': 1.2263397073452056, 'disutil_marry_sm_mal': 20.593804749350948, 'disutil_shotgun': 6.3198680053625305, 'abortion_costs': 0.024332694220526502, 'p_abortion_access': 0.9980595290220098, 'u_lost_divorce': 10.134014640329196, 'sigma_psi': 0.6656407364967266, 'sigma_psi_init': 1.2788831549969952, 'u_shift_mar': 0.3035027068097137, 'util_alp': 0.4451218545526679, 'util_kap': 11.281001265137967, 'mu_psi_init': -1.4448913321100947,
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
