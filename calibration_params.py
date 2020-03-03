#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 13:11:16 2020

@author: egorkozlov
"""

import numpy as np
from collections import OrderedDict

def calibration_params(xin=None,xfix=None):
    # NOTE: xfix overwrites xin
    
    
    # format is name: (lb,ub,xinit)
    # I am not sure if this should be ordered or not but let's do ordered
    # just in case...
    params = OrderedDict(
              u_shift_mar=(0.0,3.0,1.0),
              sigma_psi=(0.0,1.0,0.1),
              sigma_psi_mult=(1.0,8.0,3.0),
              pmeet=(0.0,1.0,0.4),
              pmeet_t=(-0.1,0.0,-0.05),
              util_alp=(0.05,4.0,1.5),
              util_kap=(0.1,8.0,1.5),
              preg_20=(0.01,0.5,0.1),
              preg_30=(0.01,0.5,0.1),
              util_qbar=(0.0,8.0,0.1),
              poutsm = (1/6,2/3,1/3)
                        )
             
    
    # update params is some dict with values was supplied
    if xin is not None:
        assert type(xin) is dict
        for key in xin:
            if type(xin[key]) is tuple:
                params[key] = xin[key]
            else:
                assert key in params
                v_old = params[key]
                params[key] = (v_old[0],v_old[1],xin[key])
                
    if xfix is not None:
        assert type(xfix) is dict
        for key in xfix:
            assert key in params
            if xin is not None and key in xin and xin[key] != xfix[key]:
                print('Warning: xfix overwrites xin')
            xval = xfix[key]
            params[key] = (xval,xval,xval)
                           
    
    keys_fixed, x_fixed = list(), list()
    
    keys, lb, ub, x0 = list(), list(), list(), list()
    
    for x in params:    
        lb_here = params[x][0]
        ub_here = params[x][1]
        x_here = params[x][2]
        
        if np.abs(ub_here - lb_here) < 1e-4: # if ub and lb are equal
            assert lb_here <= x_here <= ub_here
            keys_fixed.append(x)
            x_fixed.append(x_here)
        else:        
            keys.append(x)
            lb.append(params[x][0])
            ub.append(params[x][1])
            x0.append(params[x][2])
                      
            
    lb, ub, x_here, x_fixed = (np.array(x) for x in (lb,ub,x_here,x_fixed))
    
    def translator(x):
        # in case x has a weird dimensionality
        try:
            x.squeeze()
        except:
            pass
        
        assert len(keys) == len(x), 'Wrong lenght of x!'
        d_var = dict(zip(keys,x))
        if len(keys_fixed) > 0:
            d_fixed = dict(zip(keys_fixed,x_fixed))
            d_var.update(d_fixed)
        return d_var
    
    return lb, ub, x0, keys, translator
    
