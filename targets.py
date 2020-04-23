#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

        
def target_values(mode='high education'):
    targets = dict()
        
        
    if mode=='high education':
        
        targets['never married by age, b0'] = (0.3948181,    0.000777)
        targets['never married by age, b1'] = (-0.0507542,  0.0000895)
        targets['never married by age, b2'] = (0.0026662,   0.0000174)

    
        targets['ever kids by age, b0'] = (0.4100186,0.0007724)
        #targets['ever kids by age, b1'] = (0.0497961,0.0000801)
        #targets['ever kids by age, b2'] = (-0.0011349,0.0000177)
        
        #targets['ever kids by years after marriage, b0'] = (0.0306335,0.0029289)
        #targets['ever kids by years after marriage, b1'] = (0.1760231,0.0012216)
        #targets['ever kids by years after marriage, b2'] = (-0.0094863,0.0001077)
        
        
        targets['ever kids 2 years after marriage'] = (0.3420612,0.0022199)
        targets['ever kids 4 years after marriage'] = (0.6122718,0.0023161)
        targets['ever kids 6 years after marriage'] = (0.7739091,0.0020495)
        

        targets['divorced by years after marriage, b0'] = (-.0139552,0.0012282)
        targets['divorced by years after marriage, b1'] = (0.0191589,0.0006924)
        targets['divorced by years after marriage, b2'] = (-0.0006421,0.0000691)
        
        targets['divorced at 30 if one marriage'] = (0.0619966, 0.0005809)
        
        targets['mean x share'] = (0.4,0.001)
        
        targets['k then m by age, b0'] = (0.1199747,0.0009933)
        targets['k then m by age, b1'] = (-0.0100934,0.000355)
        targets['k then m by age, b2'] = (0.0006642,0.0000398)
        
        #targets['share of kids in new marriages, b0'] = (0.058445,0.0010725)
        #targets['share of kids in new marriages, b1'] = (0.0031754,0.0002467)
        
        targets['men ever married at 30, ratio'] =  (1.201104,0.0136505)
        targets['men divorced at 30, ratio'] = (0.5981809,0.0493313)
        targets['men ever kids at 30, ratio'] = (0.9476518,0.0151745)
        
        targets['share of kids in new marriages at 25'] = (0.0456158,0.00030936) 
        targets['share of kids in new marriages at 30'] = (0.0571774,0.00037749)
        targets['share of kids in new marriages at 35'] = (0.0867144,0.00078597)
        
        targets['divorced if k then m and one marriage'] = (0.1474027, (1/5)*0.00243)
        targets['divorced if m then k and one marriage'] = (0.0531433, (1/5)*0.0004851)
        
        
        targets['divorced with kids at 30']      = (0.0251494,0.0007414)
        targets['never married with kids at 30'] = (0.0443802,0.0007902)      
        targets['more than one mar at 40']       = (0.1225949,0.0012733)
        
        targets['in labor force at 30 if kids'] = (0.7394259,(1/5)*0.0035188)
        targets['in labor force at 30 if kids ratio'] = (0.8641976,0.0085405)
    
        
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
