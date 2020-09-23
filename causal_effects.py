#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:38:35 2020

@author: egorkozlov
"""
import numpy as np
from simulations import Agents

class Agents_FixedSim(Agents):
    # this are simulations where shocks are altered such that everyone 
    # meets partner at period meet_at
    def __init__(self,*args,meet_at=4,is_preg=True,**kwargs):
        Agents.__init__(self,*args,nosim=True,**kwargs)
        # replace shocks
        # shcoks are kids of inverted probabilities
        
        self._shocks_single_meet[:,:meet_at] = 1.0        
        self._shocks_single_preg[:,:meet_at] = 1.0
        
        self._shocks_single_meet[:,meet_at] = 0.0        
        self._shocks_single_preg[:,meet_at] = (0.0 if is_preg else 1.0)
        
        
        self.run_sim()
        
def get_cfs(mdl):
    tmeet = 5
    np.random.seed(12)
    a_p = Agents_FixedSim(mdl,T=40,meet_at=tmeet-1,no_sm=True,is_preg=True,fix_seed=False,verbose=False)
    np.random.seed(12)
    a_n = Agents_FixedSim(mdl,T=40,meet_at=tmeet-1,no_sm=True,is_preg=False,fix_seed=False,verbose=False)
    
    
    
    n_mark = a_n.state_codes['Couple and child']
    n_marnk = a_n.state_codes['Couple, no children']
    n_single = a_n.state_codes['Female, single']
    n_singlek = a_n.state_codes['Female and child']

    
    ind_agree_p = (a_p.state[:,tmeet]>=3) # index of compliers
    ind_agree_np = (a_n.state[:,tmeet]>=3)
    ind_comp = ind_agree_p & (~ind_agree_np)
    ind_always = (ind_agree_p) & (ind_agree_np)
    ind_def = (~ind_agree_p) & (ind_agree_np)
    ind_never = (~ind_agree_p) & (~ind_agree_np)
    
    
    
    print('{} compliers, {} always-takers, {} defiers, {} never-takers'.\
          format(ind_comp.sum(), ind_always.sum(), ind_def.sum(), ind_never.sum()))
    
    m_p = a_p.compute_moments()
    m_n = a_n.compute_moments()
    
    
    tc = 29
    age = tc+21
    
    
    header = r'''\hline\hline
& \textbf{Offered} & \textbf{Agreed} & \textbf{Compliers} & \textbf{Always} & \textbf{Never} \\\hline
\textit{Marry if pregnant} & $+,-$ & $+$ & $+$ & $+$ & $-$ \\
\textit{Marry if not pregnant} &  $+,-$ & $+,-$ &  $-$  & $-$ & $-$ \\ \hline'''
    
    def r(x): return 100*x.mean()
    
    prop_row = r'\textit{' + 'share in population' + r'}' + '& 100 & {:02.1f} & {:02.1f} & {:02.1f} & {:02.1f}'\
                .format(r(ind_agree_p),r(ind_comp),r(ind_always),r(ind_never)) + r'\\'
    print(header)
    print(prop_row)

    for scenario, a in zip(['pregnant','not pregnant'],[a_p,a_n]):
        is_mar = (a.state == n_mark) | (a.state == n_marnk)
        is_mark = (a.state == n_mark)
        is_sm = (a.state == n_singlek)

        ever_mar = (np.cumsum(is_mar,axis=1) > 0)
        div_now =  (ever_mar) & ((a.state==n_single) | (a.state==n_singlek))
        ever_div = (np.cumsum(div_now,axis=1) > 0)    
        ever_kid = ( np.cumsum( (a.state == n_mark) | (a.state == n_singlek),axis=1) > 0)
        more_mar = (a.nmar > 1)
        
        def e(x): return 100*x[:,tc].mean()
        def m(x): return 100*x[ind_agree_p,tc].mean()
        def c(x): return 100*x[ind_comp,tc].mean()
        def a(x): return 100*x[ind_always,tc].mean()
        def n(x): return 100*x[ind_never,tc].mean()
        
        
        for var, desc in zip([is_mar,ever_div,div_now,more_mar,is_sm],['married','ever divorced','divorced now','more than one marriage','single mother']):
            name = '{} at {}'.format(desc,age)
        #print('{}, {}: {:02.1f}'.format(scenario,name,e(var)))
        #print('{}, {}, married: {:02.1f}'.format(scenario,name,m(var)))
        #print('{}, {}, always: {:02.1f}'.format(scenario,name,a(var)))
        #print('{}, {}, compliers: {:02.1f}'.format(scenario,name,c(var)))
        #print('{}, {}, never: {:02.1f}'.format(scenario,name,n(var)))
        
            table_row = r'\textit{' + '{}, {}'.format(name,scenario) + r'} &' + ' {:02.1f} & {:02.1f} & {:02.1f} & {:02.1f} & {:02.1f}'.format(e(var),m(var),c(var),a(var),n(var)) + r'\\'
            print(table_row)
        
        
if __name__ == '__main__':
    import numpy as np
    try:
        print(len(mdl.V))
    except:
        from dill import load
        mdl = load(open('mdl.pkl','rb+'))
        
    get_cf(mdl)
    