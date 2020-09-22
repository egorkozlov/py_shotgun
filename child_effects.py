#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:38:35 2020

@author: egorkozlov
"""

from simulations import Agents

class Agents_FixedSim(Agents):
    # this are simulations where shocks are altered such that everyone 
    # meets partner at period meet_at
    def __init__(self,*args,meet_at=4,has_kids=True,**kwargs):
        Agents.__init__(self,*args,nosim=True,**kwargs)
        # replace shocks
        # shcoks are kids of inverted probabilities
        
        self._shocks_single_meet[:,:meet_at] = 1.0        
        self._shocks_single_preg[:,:meet_at] = 1.0
        
        self._shocks_single_meet[:,meet_at] = 0.0        
        self._shocks_single_preg[:,meet_at] = 1.0 # not pregnant
        
        self._shocks_planned_preg[:,meet_at:] = 0.0 if has_kids else 1.0
        
        
        self.run_sim()
        
        
if __name__ == '__main__':
    import numpy as np
    try:
        print(len(mdl.V))
    except:
        from dill import load
        mdl = load(open('mdl.pkl','rb+'))
        
        
    tmeet = 4
    np.random.seed(12)
    a_k = Agents_FixedSim(mdl,T=40,meet_at=tmeet-1,no_sm=True,has_kids=True,fix_seed=False,verbose=False)
    np.random.seed(12)
    a_nk = Agents_FixedSim(mdl,T=40,meet_at=tmeet-1,no_sm=True,has_kids=False,fix_seed=False,verbose=False)
    
    
    
    n_mark = a_k.state_codes['Couple and child']
    n_marnk = a_k.state_codes['Couple, no children']
    n_single = a_k.state_codes['Female, single']
    n_singlek = a_k.state_codes['Female and child']
    
    
    sta = tmeet+1
    fin = tmeet+9
    
    aval = 21 + np.arange(sta,fin)
    div_k  = np.zeros_like(aval,dtype=np.float64)
    div_nk = np.zeros_like(aval,dtype=np.float64)
    
    
    
    
    for h in range(sta,fin):
        i = h-sta
        comp_k = ((a_k.state[:,h]==4) & (a_nk.state[:,h]==3) & (a_nk.state[:,5]==3))
    
        for a, desc, save in zip([a_k,a_nk],['kid','no kid'],[div_k,div_nk]):
            is_mar = (a.state == n_mark) | (a.state == n_marnk)
            is_mark = (a.state == n_mark)
            is_sm = (a.state == n_singlek)
        
            ever_mar = (np.cumsum(is_mar,axis=1) > 0)
            div_now =  (ever_mar) & ((a.state==n_single) | (a.state==n_singlek))
            ever_div = (np.cumsum(div_now,axis=1) > 0)    
            ever_kid = ( np.cumsum( (a.state == n_mark) | (a.state == n_singlek),axis=1) > 0)
            more_mar = (a.nmar > 1)
            
            save[i] = 100*ever_div[comp_k,19].mean()
            
            print('for desc {} at horizon {} ever div at 40 is {}'.format(desc,h,ever_div[comp_k,19].mean()))
    
    
    adiff = np.mean(div_nk - div_k)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(aval,div_k,'-',label='always conceive')
    plt.plot(aval,div_nk,'--',label='never conceive')
    plt.legend()
    plt.title('share of ever divorced by 40\n by conception success, avg diff = {:02.1f}%'.format(adiff))
    plt.xlabel('age at the first attemped conception')
    plt.ylabel('% divorced by 40')
    plt.grid(True)
    plt.savefig('conception_divorce.pdf')
    
    

