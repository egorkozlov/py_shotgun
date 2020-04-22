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
    

    x = {'sigma_psi': 0.15917598496253035,
         'sigma_psi_mult': 5.498633340981902,
         'pmeet_0': 0.4079340929161473,
         'preg_a0': 0.02742812287263574,
         'preg_at': 0.012246548651264927,
         'preg_at2': -5.392375114762066e-05,
         'u_shift_mar': 1.7533370629319474,
         'util_alp': 0.6932061300683797,
         'util_kap': 0.710499604679927,
         'util_qbar': 0.0057462460677633245,
         'disutil_marry_sm_mal_coef': 5.4907450877050135,
         'disutil_shotgun_coef': 1.585538253532418}






    
    
    tar = target_values('high education')
    
    out, mdl, agents, res = mdl_resid(x=None,targets=tar,
                                      return_format=['distance','models','agents','scaled residuals'],
                                      #save_to='mdl.pkl',
                                      verbose=True,draw=True)
    
    mdl[0].time_statistics()
    #mdl[0].diagnostics()
                         
    print('Done. Residual in point x0 is {}'.format(out))