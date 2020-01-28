#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27
 
@author: Egor Kozlov
"""
 

from platform import system
     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
 
 
from residuals import mdl_resid
 
print('Hi!')
 
 
 
 
if __name__ == '__main__':
     
    
    out, mdl, agents, res = mdl_resid(return_format=['distance','models','agents','scaled residuals'],
                                      #load_from=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      #save_to=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      verbose=True,draw=True)
                         
    print('Done. Residual in point x0 is {}'.format(out))
     
        
