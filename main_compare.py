#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 21:54:21 2020

@author: egorkozlov
"""

import dill

from platform import system
     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
 
 
from residuals import mdl_resid
 
print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
     
    x_base = {'u_shift_mar': 0.725863356303563,
         'sigma_psi': 0.13579410386649068,
         'sigma_psi_mult': 4.270601078298901,
         'pmeet': 0.30889908284645357,
         'pmeet_t': -0.0038829067012848967,
         'util_alp': 0.6534190912803465,
         'util_kap': 1.9136130954048896,
         'preg_20': 0.03915073600027051,
         'preg_30': 0.10173830697979773,
         'high education': True}
    
    
    
    x_comp = {'u_shift_mar': 0.78012885,
         'sigma_psi': 0.22466639,
         'sigma_psi_mult': 1.34982032,
         'pmeet': 0.17557838,
         'pmeet_t': -0.00900235,
         'util_alp': 0.47896648,
         'util_kap': 1.72740958,
         'preg_20': 0.2184985,
         'preg_30': 0.33067173,
         'high education': False}
    
    # form xlist
    xlist = [x_base]
    he = 'high education'
    replace_fields = [[he],
                      ['pmeet','pmeet_t'],
                      ['pmeet','pmeet_t','preg_20','preg_30'],
                      ['sigma_psi','sigma_psi_mult'],
                      ['u_shift_mar','util_kap','util_alp'],
                      ['sigma_psi','sigma_psi_mult','u_shift_mar','util_kap','util_alp']
                     ]
    for rf in replace_fields:
        xnew = x_base.copy()
        for f in rf:
            xnew[f] = x_comp[f]
        xlist.append(xnew)
    
    
    replace_fields = ['base'] + replace_fields + ['compare']
    xlist.append(x_comp)
    
    assert False
    agents_list = list()
    for x, rf in zip(xlist,replace_fields):        
        print(rf)
        print(x)
        out, agents = mdl_resid(x=x,return_format=['distance','agents'],
                                verbose=False,draw=False)
        agents_list.append(agents)
                         
        print('Done. Residual in point x0 is {}'.format(out))
        
    dill.dump((agents_list,replace_fields,xlist),open('comparison.pkl','wb+'))            
        