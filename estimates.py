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
        
        x = {'sigma_psi': 0.3869773077864603, 'sigma_psi_init': 0.9938164169266216, 'pmeet_21': 0.26138291062442615, 'pmeet_30': 0.7186472606306563, 'pmeet_40': 0.6593653665181554, 'preg_21': 0.01355186141662465, 'preg_28': 0.026776949709968387, 'preg_35': 0.037440655536951584, 'u_shift_mar': 0.23416090232520456, 'util_alp': 0.4080762765646524, 'util_kap': 11.265720476636952, 'util_qbar': 30.846312906184785, 'disutil_marry_sm_mal': 16.579301983282175, 'disutil_shotgun': 7.775996065960521, 'abortion_costs': 7.10104593516042, 'p_abortion_access': 0.5697660445900381, 'u_lost_divorce': 13.108545716002517, 'mu_psi_init': -0.5566081252903454,
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
