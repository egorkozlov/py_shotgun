#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

        
def target_values(mode='high education'):
    targets = dict()
        
        
    if mode=='high education':
        
        targets['never married at 25'] = 0.75 
        targets['never married at 30'] = 0.38
        targets['never married at 35'] = 0.21
        targets['never married at 40'] = 0.15
    
        targets['divorced right now at 25'] = 0.057
        targets['divorced right now at 30'] = 0.084
        targets['divorced right now at 35'] = 0.11
        targets['divorced right now at 40'] = 0.15
        
        targets['no kids at 25'] = 0.90
        targets['no kids at 30'] = 0.60
        targets['no kids at 35'] = 0.34
        
        targets['no kids at 25 if married'] = 0.71
        targets['no kids at 30 if married'] = 0.39
        targets['no kids at 35 if married'] = 0.17
        
        targets['no kids 1 year after marriage'] = 0.81
        targets['no kids 2 years after marriage'] = 0.66
        targets['no kids 3 years after marriage'] = 0.51
        
        targets['mean x share'] = 0.4
        
        targets['k then m at 25'] = 0.21
        targets['k then m at 30'] = 0.12
        targets['k then m at 35'] = 0.10
        
        targets['just k & m at 25'] = 0.0058
        targets['just k & m at 30'] = 0.0081
        targets['just k & m at 35'] = 0.0046
        
    else:
        raise Exception('this mode for targets is not found')
    
    return targets