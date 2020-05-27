#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 19:43:08 2020

@author: egorkozlov
"""

def get_point(high_e):

    if high_e:
        
        x = {'sigma_psi': 0.13140335512127152,
             'sigma_psi_mult': 7.318494317792161,
             'pmeet_21': 0.31499861958165964,
             'pmeet_28': 0.60067071370404,
             'pmeet_35': 0.7472621428312332,
             'preg_21': 0.040420105606407006,
             'preg_28': 0.02296063935559003,
             'preg_35': 0.0698790264219899,
             'u_shift_mar': 0.9606447172924117,
             'util_alp': 0.3468964692327379,
             'util_kap': 0.9378088317753515,
             'util_qbar': 0.9082817629134391,
             'disutil_marry_sm_mal_coef': 8.44506597458977,
             'disutil_shotgun_coef': 0.25841666575918587,
             'taste_shock_mult': 5.476599495906876,
             'high education':True}
    
        targ_mode = 'high education'
        
    else:
        x = {'sigma_psi': 0.480991547005569,
             'sigma_psi_mult': 4.600215309475741,
             'pmeet_21': 0.48886846739080425,
             'pmeet_28': 0.17451514254902903,
             'pmeet_35': 0.5918506163621335,
             'preg_21': 0.3960300615120441,
             'preg_28': 0.115344941557771,
             'preg_35': 0.5167972159585046,
             'u_shift_mar': 1.9291671160578177,
             'util_alp': 0.47953144662272207,
             'util_kap': 0.8247164477867641,
             'util_qbar': 0.015315764308810709,
             'disutil_marry_sm_mal_coef': 7.70003003525308,
             'disutil_shotgun_coef': 0.4186198614019019,
             'taste_shock_mult': 3.1008647867466848,
             'high education':False}
        
        targ_mode = 'low education'
    
    return x, targ_mode