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
        x = {'pmeet_21': 0.45883585688212, 'pmeet_30': 0.18734805268478968, 'pmeet_40': 0.9999913332916613, 'preg_21': 0.017231341009001214, 'preg_28': 0.182654503841859, 'preg_35': 0.2978696650103295, 'util_qbar': 0.6406286963250206, 'disutil_marry_sm_mal': 19.855662599071852, 'disutil_shotgun': 3.024338048927361, 'abortion_costs': 5.761894369635252, 'p_abortion_access': 0.5081043774395521, 'u_lost_divorce': 2.222141109928797, 'sigma_psi': 0.2798538779149558, 'sigma_psi_init': 1.6916827480353398, 'u_shift_mar': 1.1207829780316174, 'util_alp': 0.33532441265560753, 'util_kap': 0.8617581929517057, 'couple_rts': 2.80540188195949,
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
