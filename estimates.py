#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.2798538779149558, 'sigma_psi_init': 1.6916827480353398, 'pmeet_21': 0.22119951722233813, 'pmeet_30': 0.6944929872458302, 'pmeet_40': 0.9401168334472371, 'preg_21': 0.006423586619122713, 'preg_28': 0.027174530214511234, 'preg_35': -0.008892724807002941, 'u_shift_mar': 0.11207829780316174, 'util_alp': 0.33532441265560753, 'util_kap': 0.8617581929517057, 'util_qbar': 1.0366661841345464, 'disutil_marry_sm_mal': 36.86138371293624, 'disutil_shotgun': 8.880815017095989, 'abortion_costs': 29.230222498956632, 'p_abortion_access': 0.15865079162079268, 'u_lost_divorce': 9.747677263576907, 'couple_rts': 2.805401881959498,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'pmeet_21': 0.45883585688212, 'pmeet_30': 0.18734805268478968, 'pmeet_40': 0.9999913332916613, 'preg_21': 0.017231341009001214, 'preg_28': 0.182654503841859, 'preg_35': 0.2978696650103295, 'util_qbar': 0.6406286963250206, 'disutil_marry_sm_mal': 19.855662599071852, 'disutil_shotgun': 3.024338048927361, 'abortion_costs': 5.761894369635252, 'p_abortion_access': 0.5081043774395521, 'u_lost_divorce': 2.222141109928797, 'sigma_psi': 0.2798538779149558, 'sigma_psi_init': 1.6916827480353398, 'u_shift_mar': 1.1207829780316174, 'util_alp': 0.33532441265560753, 'util_kap': 0.8617581929517057, 'couple_rts': 2.80540188195949,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode