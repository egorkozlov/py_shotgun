#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27
 
@author: Egor Kozlov
"""


#import warnings
#warnings.filterwarnings("error")
 

from platform import system
     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
 
 
from residuals import mdl_resid
from targets import target_values

print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
     
    x = {'u_shift_mar': 2.3355210455808297,
         'sigma_psi': 0.6053952015313737,
         'sigma_psi_mult': 3.658469721691887,
         'pmeet': 0.36236559249094785,
         'pmeet_t': -0.010441196032291739,
         'util_alp': 1.724517409748223,
         'util_kap': 1.3477223782840697,
         'preg_20': 0.04812530882453013,
         'preg_30': 0.0900267785184747,
         'sm_shift': -2.812699250993756,
         'util_qbar': 3.9922428798997034}
    
    q = x['util_qbar']
    sm = x['sm_shift']
    request = {'util_qbar':[0.0*q,0.25*q,0.5*q,0.75*q,q,1.25*q,1.5*q,1.75*q,2*q],
               'sm_shift':[0.1*sm,0.25*sm,0.5*sm,0.75*sm,sm,1.25*sm,1.5*sm,1.75*sm,2*sm],
               'preg_mult':[0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0],
               'poutsm':[0.05,0.1,0.2,1/3,1/2,2/3,0.99],
               'z_drift':[-0.4,-0.3,-0.2,-0.15,-0.1,-0.05,-0.025,0.0]
               }
    
    inp = list()
    change = list()
    
    for par, vlist in request.items():
        for val in vlist:
            change.append((par,val))
            xin = x.copy()
            xin.update({par:val})
            inp.append(('moments',xin))
            
    from p_client import compute_for_values
    
    
    print('generated {} points'.format(len(change)))
    result = compute_for_values(inp)
    
    
    from tiktak import filer
    filer('sens_results.py',{'init':x,'input':change,'output':result},True)
    
    
    
            