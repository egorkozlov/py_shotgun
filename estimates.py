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
        
        x = {'sigma_psi': 0.025092055726001986, 'sigma_psi_init': 0.057524117123893596, 'pmeet_21': 0.37414286769758254, 'pmeet_30': 1.1932756253573578, 'pmeet_40': 1.2, 'preg_21': 0.016846303459086942, 'preg_28': 0.026031363571378165, 'preg_35': 0.028147668865726228, 'u_shift_kid': 0.25811605249933434, 'util_alp': 0.6202419382122233, 'util_kap': 8.328279006156304, 'util_qbar': 19.58211188499201, 'disutil_marry_sm_mal': 11.036758925366986, 'disutil_shotgun': 0.49670652381536956, 'abortion_costs': 2.2609160794229095, 'p_abortion_access': 0.9891922518917426, 'u_lost_divorce': 1.0689656664994374, 'u_shift_couple': -0.16143081649350327,
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
        x = {'sigma_psi': 0.05247477257301193, 'sigma_psi_init': 0.058638974847773986, 'pmeet_21': 0.6383223607570659, 'pmeet_30': 0.30766246819404597, 'pmeet_40': 1.1587500618700641, 'preg_21': 0.168868268810343, 'preg_28': 0.07600966216926533, 'preg_35': 0.11428422406786634, 'u_shift_kid': 0.31694962728508075, 'util_alp': 0.4693379117706382, 'util_kap': 7.6785444482668765, 'util_qbar': 1.0222307940469122, 'disutil_marry_sm_mal': 3.933857467181001, 'disutil_shotgun': 0.6025347639597007, 'abortion_costs': 0.13283860874257264, 'p_abortion_access': 0.7321046104040729, 'u_lost_divorce': 0.8144152428939953, 'u_shift_couple': -0.24325764797700283,
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
