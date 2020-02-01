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


    lb = np.array([ 0.0,  1e-4,   0.5,  0.1,  -0.2, 0.0,  0.01, 0.05,  0.05])
    ub = np.array( [ 2.0,  0.5,  10.0,  1.0,   0.0, 1.0,   3.0,  3.0,  0.9])
    x0 = np.array( [0.5,  0.05,   2.0,  0.4, -0.05, 0.8,   0.5,  0.6,  0.3 ])
    
    
    
    ##### FIRST LET'S TRY TO RUN THE FUNCTION IN FEW POINTS
    
    print('Testing the workers...')
    from p_client import compute_for_values
    pts = [lb + rs(lb.shape)*(ub-lb) for _ in range(3)]
    pts = [('compute',x) for x in pts]    
    outs = compute_for_values(pts,timeout=3600.0)
    print('Everything worked, output is {}'.format(outs))
    
    
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
    
    
   
        

