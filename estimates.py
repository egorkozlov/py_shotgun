#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.38461165046941515,
             'sigma_psi_mult': 7.022292622575774,
             'pmeet_21': 0.23394357876540253,
             'pmeet_28': 0.6735306610248678,
             'pmeet_35': 0.717212780520049,
             'preg_21': 0.05584337164444875,
             'preg_28': 0.012522069404703729,
             'preg_35': 0.058071081606685906,
             'u_shift_mar': 0.960273879120309,
             'util_alp': 0.3780867792811883,
             'util_kap': 0.9037789352416816,
             'util_qbar': 0.8677125627698352,
             'disutil_marry_sm_mal_coef': 14.814083168430468,
             'disutil_shotgun_coef': 0.25066659188222656,
             'taste_shock_mult': 9.984701891898439,
             'high education': True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.5282645781191165,
             'sigma_psi_mult': 4.891223345649808,
             'pmeet_21': 0.48748331898178726,
             'pmeet_28': 0.1675112412415896,
             'pmeet_35': 0.4901448272323501,
             'preg_21': 0.39984814046827155,
             'preg_28': 0.08209777096634693,
             'preg_35': 0.5214583896915517,
             'u_shift_mar': 1.8867538482146373,
             'util_alp': 0.45560209855695966,
             'util_kap': 0.8513099021480224,
             'util_qbar': 0.30468765776205914,
             'disutil_marry_sm_mal_coef': 10.818160954316083,
             'disutil_shotgun_coef': 0.43228962765165996,
             'taste_shock_mult': 2.2965454818692232,
             'high education': False}
        
        targ_mode = 'low education'
    
    return x, targ_mode