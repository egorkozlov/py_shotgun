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
    
import gc
import dfols


 
from residuals import mdl_resid
from targets import target_values

print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
    

    from calibration_params import calibration_params
    
    lb, ub, x0, keys, translator = calibration_params()
    
    
    weight = 0.9
    
    import numpy as np
    
    
    x0, lb, ub = np.array(x0), np.array(lb), np.array(ub)
    x_rand = lb + np.random.random_sample(len(lb))*(ub-lb)
    x_init = x0*weight + (1-weight)*x_rand
    
    print('starting from {}'.format(x_init))
    
    
    tar = target_values('high education')
    def q(pt):
            #print('computing at point {}'.format(translator(pt)))
            try:
                ans = mdl_resid(translator(pt),return_format=['scaled residuals'])
            except:
                print('During optimization function evaluation failed at {}'.format(pt))
                ans = np.array([1e6])
            finally:
                gc.collect()
                return ans
            
            
    res=dfols.solve(q, x_init, rhobeg = 0.02, rhoend=1e-5, maxfun=120, bounds=(lb,ub),
                    scaling_within_bounds=True, objfun_has_noise=False,
                    npt = len(x_init)+5,
                    user_params={'restarts.use_restarts':True,
                                 'restarts.rhoend_scale':0.5,
                                 'restarts.increase_npt':True})
    
    
    print(res)
    print('Result is {}'.format(translator(res.x)))