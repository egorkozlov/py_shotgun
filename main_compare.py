#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 21:54:21 2020

@author: egorkozlov
"""

import dill

from platform import system
from targets import target_values

     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
 
 
from residuals import mdl_resid
 
print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

from estimates import get_point

 
if __name__ == '__main__':
     
    
    
    x_base, tm = get_point(True,read_wisdom=False) # high education
    
    tar = target_values(tm)
    
    x_comp, _ = get_point(False,read_wisdom=False) # low education
    
    
    pmeets = ['pmeet_21','pmeet_30','pmeet_40']
    ppregs = ['preg_21','preg_28','preg_35']
    abortions = ['abortion_costs','p_abortion_access']
    prefs = ['disutil_marry_sm_mal','disutil_shotgun','u_lost_divorce','util_qbar']
    # form xlist
    xlist = [x_base]
    he = ['high education'],
    replace_fields = [('trend only',he),
                      ('pmeets',pmeets),
                      ('ppregs',ppregs),
                      ('remar',['disutil_marry_sm_mal']),
                      ('stigma',['disutil_shotgun']),            
                      ('div costs',['u_lost_divorce']),
                      ('qbar',['util_qbar']),
                      ('abortions',abortions),
                      ('all probs',pmeets+ppregs),
                      ('all prefs',prefs),
                      ('trend and all prefs',he+prefs),
                      ('trend and all probs',he+pmeets+ppregs),   
                      ('trend and abortions',he+abortions),
                      ('all but trend',prefs+pmeets+ppregs+abortions)
                     ]
    for name, rf in replace_fields:
        xnew = x_base.copy()
        for f in rf:
            xnew[f] = x_comp[f]
        xlist.append(xnew)
    
    
    replace_fields = ['base'] + replace_fields + ['compare']
    xlist.append(x_comp)
    
    
    names = ['college'] + [f[0] for f in replace_fields] + ['high school']
    
    names = ['educ comparison: '+ n for n in names]
    
    for x, name in zip(xlist,names):        
        print(name)
        print(x)
        out,  mom = mdl_resid(x=x,targets=tar,
                                          return_format=['distance','moments'],
                                          verbose=False,draw=False,cs_moments=False,
                                          save_to ='mdl for {}'.format(name),
                                          moments_save_name = name,
                                          moments_repeat=5)
                         
        print('Done. Residual in point x0 is {}'.format(out))
        
        