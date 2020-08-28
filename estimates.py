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
        x = {'sigma_psi': 0.8783208275619998, 'sigma_psi_init': 4.841913371090404, 'pmeet_21': 0.4429662745957216, 'pmeet_30': 0.09976862644045448, 'pmeet_40': 0.6556590675095522, 'preg_21': 0.09259183950519333, 'preg_28': 0.07428802631882996, 'preg_35': -0.049794628003451157, 'u_shift_mar': 1.7959460915784782, 'util_alp': 0.4593118057932805, 'util_kap': 1.1639238152564353, 'util_qbar': 0.6935697519247619, 'disutil_marry_sm_mal': 79.68148741525432, 'disutil_shotgun': 5.315001707974356, 'abortion_costs': 5.516265535524452, 'p_abortion_access': 0.18163829572926798, 'u_lost_divorce': 1.184990222239299,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode