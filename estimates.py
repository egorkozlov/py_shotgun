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
        
        x = {'sigma_psi': 0.6542948852950291, 'sigma_psi_init': 1.2686939847399528, 'pmeet_21': 0.34270461302365784, 'pmeet_30': 0.895549399258434, 'pmeet_40': 0.5174108737115869, 'preg_21': 0.01292270977905803, 'preg_28': 0.0344102613105532, 'preg_35': 0.03221323578182611, 'u_shift_mar': 0.28499330541232826, 'util_alp': 0.3735536978974196, 'util_kap': 11.637017584856121, 'util_qbar': 37.63218891154509, 'disutil_marry_sm_mal': 20.36893014147807, 'disutil_shotgun': 5.705192721941784, 'abortion_costs': 6.472136197200477, 'p_abortion_access': 0.8752395298709243, 'u_lost_divorce': 14.902577701560604, 'mu_psi_init': -1.43931830126317,
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
        x = {'pmeet_21': 0.728698619306192, 'pmeet_30': 0.23310107046750475, 'pmeet_40': 0.9629955398495929, 'preg_21': 0.03687725755909954, 'preg_28': 0.0743925695472168, 'preg_35': 0.19133138664351215, 'util_qbar': 25.181454069571686, 'disutil_marry_sm_mal': 14.040986098063753, 'disutil_shotgun': 5.8766374053947485, 'abortion_costs': 4.317887395330901, 'p_abortion_access': 0.6662455812243002, 'u_lost_divorce': 5.439790568560208, 'sigma_psi': 0.3869773077864603, 'sigma_psi_init': 0.9938164169266216, 'u_shift_mar': 0.23416090232520456, 'util_alp': 0.4080762765646524, 'util_kap': 11.265720476636952, 'mu_psi_init': -0.5566081252903454,
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
