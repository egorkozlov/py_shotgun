#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 18:14:17 2020

@author: egorkozlov
"""
import matplotlib.pyplot as plt  
import numpy as np

def v_graphs(mdl):
    
    
    
    # this plots graphs for theta
    plt.figure()
    
    for t in [5,10]:
        x, v = mdl.get_graph_values(t=t,state='Couple, no children',dec=True,fun='thetas',
                                    iassets=4,itheta=slice(None),iexo=(4,3,6))
        plt.plot(x,v,label='t = {}'.format(t))
    plt.plot(x,x,label='def')
    plt.legend()
    plt.title('Actual and renegotiated theta') 
    
    
    
    
    plt.figure()
    
    for t in [0,5,10,15,20,25,30,35]:
        x, v = mdl.get_graph_values(t=t,state='Couple and child',dec=False,fun='s',
                                    iassets=slice(None),itheta=5,iexo=(4,3,6))
        plt.plot(x,v,label='t = {}'.format(t))
        
    plt.legend()
    plt.title('Savings as function of assets') 
    
    
    plt.figure()
    xv = np.zeros(35)
    vv = np.zeros(35,dtype=np.float64)
    
    for t in range(35):
        
        xv[t] = t
        _, vv[t] = mdl.get_graph_values(t=t,state='Couple and child',dec=False,fun='s',
                                    iassets=10,itheta=5,iexo=(4,3,6))
    
    plt.plot(xv,vv)        
    plt.title('Savings over time') 
    
    
    
    
    