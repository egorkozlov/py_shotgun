#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.3606833154484394, 'sigma_psi_init': 3.031242467682322, 'pmeet_21': 0.2981951262856033, 'pmeet_28': 0.7468405296611746, 'pmeet_35': 0.9700424253585994, 'preg_21': 0.023056886381180183, 'preg_28': 0.03585409029595745, 'preg_35': 0.013599676361164492, 'u_shift_mar': 1.1114058257254988, 'util_alp': 0.2509204297743008, 'util_kap': 0.8802848731701842, 'util_qbar': 1.157339328421708, 'disutil_marry_sm_mal': 76.65446606843315, 'disutil_shotgun': 7.079428268293057, 'abortion_costs': 33.74015195752216, 'p_abortion_access': 0.9974677343975493, 'u_lost_divorce': 8.105777902557344,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.39953601382428866, 'sigma_psi_init': 3.0725683838557276, 'pmeet_21': 0.432991116061913, 'pmeet_28': 0.26180839230577446, 'pmeet_35': 0.5691992566974007, 'preg_21': 0.170975611759682, 'preg_28': 0.21779042931083295, 'preg_35': 0.17575004499028102, 'u_shift_mar': 1.4998917495906658, 'util_alp': 0.4179831829995812, 'util_kap': 1.3061403769026076, 'util_qbar': 0.003578373701123387, 'disutil_marry_sm_mal': 51.48144648331411, 'disutil_shotgun': 0.15849135146243648, 'abortion_costs': 40.05274285709042, 'p_abortion_access': 0.9759747406920513, 'u_lost_divorce': 0.5704761664715526,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode