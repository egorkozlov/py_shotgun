#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.10031363321851126, 'sigma_psi_init': 0.1427541210832547, 'pmeet_21': 0.11256508666615017, 'pmeet_30': 0.2097190890621539, 'pmeet_40': 0.34536994626423106, 'preg_21': 0.01383288422378677, 'preg_28': 0.04441030757493903, 'preg_35': 0.005721727876573733, 'u_shift_mar': 1.3872323315378159, 'util_alp': 0.35785582638744484, 'util_kap': 0.7571415322805706, 'util_qbar': 1.6972149130997232, 'disutil_marry_sm_mal': 43.335194296938695, 'disutil_shotgun': 4.173717574930719, 'abortion_costs': 37.53003145645212, 'p_abortion_access': 0.22711552946941432, 'u_lost_divorce': 5.296394363464528,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.8783208275619998, 'sigma_psi_init': 4.841913371090404, 'pmeet_21': 0.4429662745957216, 'pmeet_30': 0.09976862644045448, 'pmeet_40': 0.6556590675095522, 'preg_21': 0.09259183950519333, 'preg_28': 0.07428802631882996, 'preg_35': -0.049794628003451157, 'u_shift_mar': 1.7959460915784782, 'util_alp': 0.4593118057932805, 'util_kap': 1.1639238152564353, 'util_qbar': 0.6935697519247619, 'disutil_marry_sm_mal': 79.68148741525432, 'disutil_shotgun': 5.315001707974356, 'abortion_costs': 5.516265535524452, 'p_abortion_access': 0.18163829572926798, 'u_lost_divorce': 1.184990222239299,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode