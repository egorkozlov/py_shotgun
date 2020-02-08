#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""





if __name__ == '__main__':
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass


from platform import system
    
if system() != 'Darwin' and system() != 'Windows':   
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'


import numpy as np
from numpy.random import random_sample as rs
from tiktak import tiktak
print('Hi!')

from residuals import mdl_resid

if __name__ == '__main__':
    
    
    #Build  data moments and pickle them
    
    
    
    #Initialize the file with parameters

    lb = np.array(   [ 0.0,         0.01,        0.5,         0.1,        -0.1,        0.5,         0.5,         0.01,       0.01,       0.00])
    ub = np.array(   [ 3.0,         0.5,         8.0,         1.0,        0.0,         3.0,         3.0,         0.5,        0.5,        4.00])
    xdef = np.array( [ 1.47052128,  0.31739663,  2.2436033 ,  0.2004341 , -0.00240084, 1.68427564,  1.7914976 ,  0.30045406, 0.01858392, 0.01])

    
    
    ##### FIRST LET'S TRY TO RUN THE FUNCTION IN FEW POINTS
    
    #print('Testing the workers...')
    #from p_client import compute_for_values
    #pts = [lb + rs(lb.shape)*(ub-lb) for _ in range(3)]
    #pts = [('compute',x) for x in pts]    
    #outs = compute_for_values(pts,timeout=3600.0)
    #print('Everything worked, output is {}'.format(outs))
    
    
    print('')
    print('')
    print('running tic tac...')
    print('')
    print('')
    
   

    #Tik Tak Optimization
    param=tiktak(400,400,20,lb,ub,mdl_resid,tole=1e-3,nelder=False,refine=False,
                 skip_local=True,skip_global=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
    #Now Re do the computation with graphs!
    out, mdl = mdl_resid(param[1],return_format=['distance','model'],calibration_report=False,
                         verbose=True,draw=True)
    
    
   
        

