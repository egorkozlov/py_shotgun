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
        
        x = {'sigma_psi': 0.10340355168700123, 'sigma_psi_init': 0.28017722631138803, 'pmeet_21': 0.26993007994477947, 'pmeet_30': 1.0490562026238832, 'pmeet_40': 1.0342929278856385, 'preg_21': 0.02874577432003439, 'preg_28': 0.024417898585832182, 'preg_35': 0.022279437447616027, 'u_shift_kid': 0.35605266302313615, 'util_alp': 0.7025642892133058, 'util_kap': 8.125707413107408, 'util_qbar': 12.745842110837797, 'disutil_marry_sm_mal': 10.134832930047057, 'disutil_shotgun': 0.8510908959996006, 'abortion_costs': 3.4476511198411384, 'p_abortion_access': 0.8926511246646444, 'u_lost_divorce': 2.4279083763158056, 'u_shift_couple': -0.3261368368243413,
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
        x = {'sigma_psi': 0.186842639244787, 'sigma_psi_init': 0.4128915613753357, 'pmeet_21': 0.6431225690807463, 'pmeet_30': 0.3368238760736062, 'pmeet_40': 1.110476816634265, 'preg_21': 0.0898662454344269, 'preg_28': 0.12145396060498559, 'preg_35': 0.08559982472302477, 'u_shift_kid': 0.3866008674665175, 'util_alp': 0.5473655316173959, 'util_kap': 7.879383388684942, 'util_qbar': 2.1745937042903325, 'disutil_marry_sm_mal': 9.107458565157478, 'disutil_shotgun': 1.1608982210461134, 'abortion_costs': 1.7460179105235598, 'p_abortion_access': 0.9354312883811831, 'u_lost_divorce': 2.252471467060377, 'u_shift_couple': -0.5081838945110324,
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
