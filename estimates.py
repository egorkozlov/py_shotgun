#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.3468795175077058, 'sigma_psi_init': 3.0036894143144446, 'pmeet_21': 0.32466290299554007, 'pmeet_30': 0.7503841331742629, 'pmeet_40': 0.9707518006465453, 'preg_21': 0.023486509188282753, 'preg_28': 0.033696244677538734, 'preg_35': 0.015412565551107149, 'u_shift_mar': 1.1368287856062929, 'util_alp': 0.2665588052372776, 'util_kap': 0.8499122688236413, 'util_qbar': 1.1664837957597296, 'disutil_marry_sm_mal': 73.11952049466217, 'disutil_shotgun': 7.162567149411041, 'abortion_costs': 32.57253104022974, 'p_abortion_access': 0.9206031442215913, 'u_lost_divorce': 8.124964732039244,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.6673685297096009, 'sigma_psi_init': 3.6205110483031033, 'pmeet_21': 0.4309613549586911, 'pmeet_30': 0.19312711918546246, 'pmeet_40': 0.6835834989287874, 'preg_21': 0.06377808161612711, 'preg_28': 0.19535538035750696, 'preg_35': 0.12337894101313086, 'u_shift_mar': 1.6703044676180692, 'util_alp': 0.476084153505717, 'util_kap': 1.3146941483772354, 'util_qbar': 0.112008848477948, 'disutil_marry_sm_mal': 62.91356489051875, 'disutil_shotgun': 2.7095830973085238, 'abortion_costs': 49.911811812631406, 'p_abortion_access': 0.9963075177010862, 'u_lost_divorce': 2.44039174374095,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode