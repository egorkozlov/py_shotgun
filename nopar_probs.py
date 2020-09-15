#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 12:00:01 2020

@author: egorkozlov
"""

from simulations import Agents

import dfols
import numpy as np
    

class AgentsEst(Agents):
    def __init__(self,*args,**kwargs):
        Agents.__init__(self,*args,nosim=True,**kwargs)
        self.define_targets()
        self.run_sim()
        self.get_shares()
        
        
        
        
    def simulate_t(self,t):
        self.anext(t)
        self.iexonext(t) 
        self.statenext(t)
    
    
    def simulate_npsolve(self,t):
        try:
            pmeet0 = self.pmeet_exo[:t] 
            ppreg0 = self.ppreg_exo[:t]
        except:
            pmeet0 = None
            ppreg0 = None
        
        
        
        # this is horrible but works somehow...
        # the way k_m is defined makes it impossible to simulate with diff't
        # probabilities and the same t...
        
        km_save = self.k_m.copy()
        mk_save = self.m_k.copy()
        nmar_save  = self.nmar.copy()
        kmtrue_save = self.k_m_true.copy()
        
        
        
        s_dta = self.ss_val[t]
        kf_dta = self.kf_val[t]
    
        
        def get_residuals(x):
            pmeet_here = x[0]
            ppreg_here = x[1]
            
            
            self.k_m = km_save.copy()
            self.m_k = mk_save.copy()
            self.km_ture = kmtrue_save.copy()
            self.nmar = nmar_save.copy()
            
            
            self.pmeet_exo = [pmeet_here] if pmeet0 is None else pmeet0 + [pmeet_here]
            self.ppreg_exo = [ppreg_here] if ppreg0 is None else ppreg0 + [ppreg_here]
            
            self.anext(t)
            self.iexonext(t) 
            self.statenext(t)
            
            
            s_mdl = (self.state[:,t+1] == self.state_codes['Female, single']).mean() #((self.i_singles())[:,t+1]).mean()
            
            
            kf_mdl = (self.k_m[:,t+1] & (self.nmar[:,t+1] == 1)).mean() #(self.i_km())[:,t+1].mean()
            
            
            #print('targets are : {}, {}'.format(kf_dta,s_dta))
            #print('model values are : {}, {}'.format(kf_mdl,s_mdl))
            
            resid = np.array([s_mdl - s_dta,kf_mdl-kf_dta])
            print('t is {}, x is {}, resid is {}'.format(t,x,np.sum(resid**2)))
            
            return resid
        
        
        
        if t == 0:
            xinit = np.array([0.1,0.1])
        else:
            xinit = np.array([self.pmeet_exo[t-1],self.ppreg_exo[t-1]])
        
    
        res = dfols.solve(get_residuals,xinit,rhoend=1e-3,bounds=(np.array([0.0,0.0]),np.array([1.0,1.0])))
        print('t = {}'.format(t))
        print(res)
    
        get_residuals(res.x) # this writes the probabilities into the object
        
    
    
    def simulate(self,rep=1):
        # this first estimates probabilities then spits them out
        for t in range(self.T-1):      
            print(t)
            try:
                for _ in range(rep): self.simulate_npsolve(t)
                if self.verbose: self.timer('Simulations, estimation')
            except IndexError:
                break
            
        print('estimation done!')
        ppreg = np.array(self.ppreg_exo)
        pmeet = np.array(self.pmeet_exo)
        
        print(ppreg)
        print(pmeet)
        
        ppreg_int = self.interpolate_last(ppreg,nlast=3)
        pmeet_int = self.interpolate_last(pmeet,nlast=3)
        print(ppreg_int)
        print(pmeet_int)
        
        self.ppreg_exo = ppreg_int
        self.pmeet_exo = pmeet_int
        
        print('simulating with estimated probabilities')
        for t in range(self.T-1):
            self.simulate_t(t)
            if self.verbose: self.timer('Simulations, iteration')
        print('done')
        
        
    def interpolate_last(self,xin,nlast=1):
        nx = xin.size
        xin_int = np.zeros(self.T,dtype=np.float64)
        xin_int[:nx] = xin
        xin_int[nx:] = xin[-nlast:].mean()
        return xin_int
        
    def interpolate_average(self,xin):
        nx = xin.size
        xin_int = np.zeros(self.T,dtype=np.float64)
        xin_int[:nx] = xin
        xin_int[nx:] = xin.mean()
        return xin_int
        
    def interpolate_quadratic(self,xin):
        reg_y = xin
        nx = xin.size
        tval = np.arange(nx,dtype=np.float64)
        
        try:
            pol = np.polyfit(tval,reg_y,2)
        except:
            pol = [0,0,0]
            
        b0 = pol[2]
        b1 = pol[1]
        b2 = pol[0]
        
        tval_int = np.arange(self.T,dtype=np.float64)
        xin_int = np.zeros(self.T,dtype=np.float64)
        xin_int[:nx] = xin
        xin_int[nx:] = np.clip(b0 + b1*tval_int[nx:] + b2*tval_int[nx:],0.0,1.0)
        
        return xin_int
            
        
    def define_targets(self):
        from targets import all_targets
        tar = all_targets('high education')
        
        
        kf_val = [tar['k then m in population at {}'.format(a)][0] for a in range(22,36)]
        #kf_val = [kf_val[1] - (kf_val[2] - kf_val[1])] + kf_val # interpolate
        
        ss_val = [(tar['divorced and no kids in population at {}'.format(a)][0] +
                   tar['never married and no kids in population at {}'.format(a)][0])for a in range(23,36)]
        
        ss_val = [ss_val[0] - (ss_val[1] - ss_val[0])] + ss_val
        
        self.ss_val = np.array(ss_val)
        self.kf_val = np.array(kf_val)
        
        
    def i_singles(self):
        return (self.state == self.state_codes['Female, single'])
    
    def i_km(self):
        return (self.k_m & (self.nmar == 1))
        
        
    def get_shares(self):
        self.share_singles = (self.i_singles()).mean(axis=0)
        self.share_km = (self.i_km()).mean(axis=0)
         
        
if __name__ == '__main__':
    import numpy as np
    try:
        print(len(mdl.V))
    except:
        from dill import load
        mdl = load(open('mdl.pkl','rb+'))
        
        
    np.random.seed(12)

    o = AgentsEst(mdl,T=30,verbose=False)
    ss_val = o.ss_val
    kf_val = o.kf_val
    
    
    
    
    '''
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
    '''
    
    '''
    reference values:
    pmeet_np = [0.22565941308804738,0.3338372741573408,
                0.32564675711745394,0.4380962723617894]
    ppreg_np = [0.02788504431884112,0.018360629740272073,
                 0.02202381141320635,0.02776843206522462]
    '''

    
    
    