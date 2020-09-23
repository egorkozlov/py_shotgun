#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:09:22 2020

@author: egorkozlov
"""

import numpy as np

def fun(x): return np.mean(x)

import matplotlib.pyplot as plt


def plot_sim_graphs(self,Tbeg=4,Tend=19):
    
    psi_nokid = np.zeros((Tend,),dtype=np.float64)
    psi_kid = np.zeros((Tend,),dtype=np.float64)
    
    tval = 21 + np.arange(Tend,dtype=np.float64)
    
    psi_kid[0] = None
    psi_nokid[0] = None
    
    
    nmar = self.state_codes['Couple, no children']
    nmark = self.state_codes['Couple and child']
    
    for t in range(1,Tend):
        # pick people who just transitioned
        
        mar_t = (self.state[:,t] == nmar)
        mark_t = ((self.state[:,t] == nmark))
        
        psi_nokid[t] = fun( self.psi_couple[mar_t,t] ) if np.any(mar_t) else None
        psi_kid[t]   = fun( self.psi_couple[mark_t,t] )   if np.any(mark_t) else None
        
    plt.plot(tval[Tbeg:],psi_nokid[Tbeg:],label='No kids')
    plt.plot(tval[Tbeg:],psi_kid[Tbeg:],label='Kids')
    plt.xlabel('female age')
    plt.ylabel('match quality ($\psi$)')
    plt.title('Mean match quality of couples \n by the presence of the kids')
    plt.grid(True)
    plt.legend()
    plt.savefig('risks_fertility.pdf')
