#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

        
def target_values(mode='high education'):
    targets = dict()
        
        
    if mode=='high education':
        
        targets['never married at 25'] = 0.7526164 
        targets['never married at 30'] = 0.3789449
        targets['never married at 35'] = 0.2119567
        targets['never married at 40'] = 0.1519036

    
        targets['divorced right now at 25'] = 0.0407639
        targets['divorced right now at 30'] = 0.0649885
        targets['divorced right now at 35'] = 0.0908754
        targets['divorced right now at 40'] = 0.1231106

        
        targets['no kids at 25'] = 0.9027082
        targets['no kids at 30'] = 0.5998235
        targets['no kids at 35'] = 0.3383375
        
        targets['no kids at 25 if married'] = 0.7139946
        targets['no kids at 30 if married'] = 0.3992189
        targets['no kids at 35 if married'] = 0.177255
        
        targets['no kids 1 year after marriage'] = 0.8111601
        targets['no kids 2 years after marriage'] = 0.6579388
        targets['no kids 3 years after marriage'] = 0.5090055

        
        targets['mean x share'] = 0.4
        
        targets['k then m at 25'] = 0.2069427
        targets['k then m at 30'] = 0.1156287
        targets['k then m at 35'] = 0.1011743
        
        targets['just k & m at 25'] = 0.0058392
        targets['just k & m at 30'] = 0.0081263
        targets['just k & m at 35'] = 0.0046128
        
        targets['divorced if k then m and one marriage'] = 0.147
        targets['divorced if m then k and one marriage'] = 0.053
        
        
        targets['divorced with kids at 30']      = 0.0251494
        targets['divorced never kids at 30']     = 0.0398391
        targets['never married with kids at 30'] = 0.0443802
        targets['more than one mar at 40']       = 0.1225949
    
        targets['std earnings at 24, female'] = 0.41251325
        targets['std earnings at 30, female'] = 0.42495642
        
        targets['log earnings coef at 25'] = -0.0177554
        targets['log earnings coef at 30'] =  0.0045791
        
        #targets['spouse log coef at 25'] = 0.1853794
        #targets['spouse log coef at 40'] = 0.1394379
        #targets['spouse log coef 1 year after'] = 0.2261541
    
        
    elif mode=='low education':
    
        targets['never married at 25'] = 0.63
        targets['never married at 30'] = 0.43
        targets['never married at 35'] = 0.30
        targets['never married at 40'] = 0.22
    
        targets['divorced right now at 25'] = 0.20
        targets['divorced right now at 30'] = 0.22
        targets['divorced right now at 35'] = 0.25
        targets['divorced right now at 40'] = 0.28
        
        targets['no kids at 25'] = 0.50
        targets['no kids at 30'] = 0.34
        targets['no kids at 35'] = 0.25
        
        targets['no kids at 25 if married'] = 0.24
        targets['no kids at 30 if married'] = 0.15
        targets['no kids at 35 if married'] = 0.11
        
        targets['no kids 1 year after marriage'] = 0.36
        targets['no kids 2 years after marriage'] = 0.30
        targets['no kids 3 years after marriage'] = 0.24
        
        targets['mean x share'] = 0.4
        
        targets['k then m at 25'] = 0.44
        targets['k then m at 30'] = 0.37
        targets['k then m at 35'] = 0.33
        
        targets['just k & m at 25'] = 0.0194
        targets['just k & m at 30'] = 0.0143
        targets['just k & m at 35'] = 0.0082
        
        targets['divorced if k then m and one marriage'] = 0.172
        targets['divorced if m then k and one marriage'] = 0.139
        
        
        
        
    else:
        raise Exception('this mode for targets is not found')
    
    return targets