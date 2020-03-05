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
     
    x_base = {'u_shift_mar': 1.4622468191770523,
             'sigma_psi': 0.13341343588771484,
             'sigma_psi_mult': 2.7114359061341653,
             'pmeet': 0.16330434124234866,
             'pmeet_t': 0.0,
             'util_alp': 1.5035902442038718,
             'util_kap': 7.862043282658021,
             'preg_20': 0.10034665087728878,
             'preg_30': 0.01,
             'util_qbar': 1.7158729365532353,
             'disutil_marry_sm_mal_coef': 8.047422930867116,
             'high education': True}
    
    
    
    x_comp = {'u_shift_mar': 1.4622468191770523,
             'sigma_psi': 1.97568527e-01,
             'sigma_psi_mult': 2.60091410,
             'pmeet': 3.33322285e-01,
             'pmeet_t': -3.73953602e-04,
             'util_alp': 1.5035902442038718,
             'util_kap': 7.862043282658021,
             'preg_20': 3.07441811e-01,
             'preg_30': 1.07009581e-01,
             'util_qbar': 3.44970307e+00,
             'disutil_marry_sm_mal_coef': 9.89430509e+00,
             'high education': False}
    
    # form xlist
    xlist = [x_base]
    
    he = 'high education'
    
    replace_fields = [[he],
                      ['pmeet','pmeet_t'],
                      ['preg_20','preg_30'],                      
                      ['sigma_psi','sigma_psi_mult'],
                      ['util_qbar'],
                      ['disutil_marry_sm_mal_coef'],
                      ['pmeet','pmeet_t','preg_20','preg_30'],
                      ['util_qbar','disutil_marry_sm_mal_coef'],
                      ['util_qbar',he]
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
        