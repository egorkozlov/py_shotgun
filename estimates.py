#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.2798538779149558, 'sigma_psi_init': 1.6916827480353398, 'pmeet_21': 0.22119951722233813, 'pmeet_30': 0.6944929872458302, 'pmeet_40': 0.9401168334472371, 'preg_21': 0.006423586619122713, 'preg_28': 0.027174530214511234, 'preg_35': -0.008892724807002941, 'u_shift_mar': 1.1207829780316174, 'util_alp': 0.33532441265560753, 'util_kap': 0.8617581929517057, 'util_qbar': 1.0366661841345464, 'disutil_marry_sm_mal': 36.86138371293624, 'disutil_shotgun': 8.880815017095989, 'abortion_costs': 29.230222498956632, 'p_abortion_access': 0.15865079162079268, 'u_lost_divorce': 9.747677263576907, 'couple_rts': 2.805401881959498,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.07262028069004121, 'sigma_psi_init': 0.5975384122500378, 'pmeet_21': 0.32747839443931864, 'pmeet_30': 0.040400190976707655, 'pmeet_40': 0.6472740842394364, 'preg_21': 0.1617772402053484, 'preg_28': 0.13779244582019745, 'preg_35': 0.11488115937628847, 'u_shift_mar': 2.0, 'util_alp': 0.01, 'util_kap': 0.7371085628065419, 'util_qbar': 1.6315860400437356, 'disutil_marry_sm_mal': 64.65098573395561, 'disutil_shotgun': 2.7831665371545045, 'abortion_costs': 0.0, 'p_abortion_access': 0.6288965038601756, 'u_lost_divorce': 0.0, 'couple_rts': 0.029883097408165055,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode