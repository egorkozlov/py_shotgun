#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 16:54:45 2020

@author: egorkozlov
"""
import numpy as np

def check_value_functions(self):
    for t in range(self.setup.pars['T']):
        Vt = self.V[t]
        for key in Vt:
            a_diff = np.diff(Vt[key]['V'],axis=0)
            monotonic = (a_diff>0)
            if np.all(monotonic):
                print('at t = {} for {} everything is monotonic'.format(t,key))                
            else:
                print('at t = {} for {} share of {} is not monotnoic'.format(t,key,np.mean(~monotonic)))
                #assert False
                
                
    