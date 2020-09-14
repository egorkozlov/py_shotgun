#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:00:01 2020

@author: egorkozlov
"""

from simulations import Agents

class AgentsEst(Agents):
    def __init__(self,*args,**kwargs):
        Agents.__init__(self,*args,nosim=True,**kwargs)
        self.simulate()
        self.get_shares()
        
        
        
        
    def simulate_t(self,t):
        self.anext(t)
        self.iexonext(t) 
        self.statenext(t)
    
    def simulate(self,rep=1):
        for t in range(self.T-1):      
            for _ in range(rep): self.simulate_t(t)
            
            
        
            if self.verbose: self.timer('Simulations, iteration')
        
        
    def get_shares(self):
        self.share_singles = (self.state == self.state_codes['Female, single']).mean(axis=0)
        self.share_km = (self.k_m & (self.nmar == 1)).mean(axis=0)
         
        
if __name__ == '__main__':
    import numpy as np
    try:
        print(len(mdl.V))
    except:
        from dill import load
        mdl = load(open('mdl.pkl','rb+'))
        
    
    
    from targets import all_targets
    tar = all_targets('high education')
    
    share_km_tar = tar['k then m in population at 22'][0]
    share_singles_tar = 1 - tar['ever married at 22'][0]
    
    kf_val = [tar['k then m in population at {}'.format(a)][0] for a in range(22,36)]
    #kf_val = [kf_val[1] - (kf_val[2] - kf_val[1])] + kf_val # interpolate
    
    ss_val = [(tar['divorced and no kids in population at {}'.format(a)][0] +
               tar['never married and no kids in population at {}'.format(a)][0])for a in range(23,36)]
    
    ss_val = [ss_val[0] - (ss_val[1] - ss_val[0])] + ss_val
    
    ss_val = np.array(ss_val)
    kf_val = np.array(kf_val)
    
    
    import dfols
    
    
    pmeet_np = []
    ppreg_np = []
    
    
    for t in range(4):
    
        print('t = {}'.format(t))
        
        def resid(x):
            np.random.seed(12)
            a = AgentsEst(mdl,T=t+2,verbose=False,
                      pmeet_exo=np.array(pmeet_np + [x[0]]),
                      ppreg_exo=np.array(ppreg_np + [x[1]]))
            r = np.array([a.share_singles[-1] - ss_val[t],
                          a.share_km[-1] - kf_val[t]])
            print('x is {}, resid is {}'.format(x,np.sum(r**2)))
            return r
        
        if t == 0:
            xinit = np.array([0.6,0.1])
        else:
            xinit = np.array([pmeet_np[-1],ppreg_np[-1]])
            
        res = dfols.solve(resid,xinit,rhoend=1e-3)
        
        pmeet_np = pmeet_np + [res.x[0]]
        ppreg_np = ppreg_np + [res.x[1]]
        
    
    
    
    
    
    