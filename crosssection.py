#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:13:18 2020

@author: egorkozlov
"""

import numpy as np
from simulations import Agents

class CrossSection:
    def __init__(self,Mlist,*,age_distribution,N_total=15000,verbose=False,
                                     fix_seed=True,verbose_sim=False,**kwargs):
        # age distribution is an array of frequencies. N_total independent 
        # agents are simulated according to this distribution
        
        if fix_seed: np.random.seed(18)
        
        self.verbose = verbose
        self.N_total = N_total
        self.verbose_sim = verbose_sim
        
    
        
        age_distribution = np.array(age_distribution)
        
        age_dist_normailzed = age_distribution/age_distribution.sum()
        age_dist_cumulative = np.cumsum(age_dist_normailzed)
        
        self.age_dist = age_dist_normailzed
        
        a_pick = (np.random.random_sample(N_total)[:,None] > age_dist_cumulative[None,:])
        assert not np.any(a_pick[:,-1])
        t_pick = a_pick.sum(axis=1) # at which age we observe each of N_total individuals
        t_pick = np.sort(t_pick)
        self.t = t_pick
        
        self.Tmax = age_distribution.size
        self.tlist = np.arange(self.Tmax)
        self.num_t = np.array([np.sum(t_pick==t) for t in self.tlist])
        
        def agents_create(n,t):
            if verbose: print('creating {} agents for t = {}'.format(n,t))
            return Agents(Mlist,N=n,T=t+2,fix_seed=False,verbose=verbose_sim,**kwargs)
        
        
        self.agents_t = [agents_create(n,t) for n,t in zip(self.num_t,range(self.Tmax))]
        self.state_codes = self.agents_t[0].state_codes.copy()
        self.female = self.agents_t[0].female
        self.parse()
        c, o, m = self.get_marriage_stats()
        self.marriage_stats = {'Counts':c,'Offers':o,'Marriages':m}
        #del(self.agents_t)
        
    def parse(self):
        agents_base = self.agents_t[-1]
        attr_list = agents_base.__dir__()
        for attr in attr_list:
            thing = getattr(agents_base,attr)
            
            try:
                assert thing.ndim==2
                assert thing.shape[0] == self.num_t[-1]        
                assert not attr[0] == '_'
                #if self.verbose: print('{} was picked'.format(attr))
            except:
                continue
            
            assert not hasattr(self,attr) # runs only once!
            
            thing_cs = np.empty(self.N_total,thing.dtype)
            for t in self.tlist:
                thing_cs[self.t==t] = getattr(self.agents_t[t],attr)[:,t]
                
            setattr(self,attr,thing_cs)
            #if self.verbose: print('{} was parsed'.format(attr))
    
    def get_marriage_stats(self):
        c = dict()
        o = dict()
        m = dict()
        
        for t in self.tlist:
            ci, oi, mi = self.agents_t[t].marriage_stats()
            for d, di in zip((c,o,m),(ci,oi,mi)):
                for key in di:  
                    dkt = di[key][t]
                    try:
                        d[key][t] = dkt
                    except KeyError:
                        d[key] = np.zeros(self.Tmax,dtype=np.int32)
                        assert t == 0
                        d[key][0] = dkt
        
        c.update({'Total':self.num_t})
        return c, o, m
            
            
    def compute_cs_moments(self):
        moments = dict()
        
        n_mark = self.state_codes['Couple and child']
        n_marnk = self.state_codes['Couple, no children']
        n_single = self.state_codes['Female, single'] if self.female else self.state_codes['Male, single']
        n_singlek = self.state_codes['Female and child']
        
        
        have_kid = (self.state == n_mark) | (self.state == n_singlek)
        couple_kid = (self.state == n_mark)
        married = (self.state == n_mark) | (self.state == n_marnk)
        divorced = (~married) & (self.nmar>0)
        
        one_mar = (self.nmar == 1)
        
        age = self.t + 21
        
        for t in range(1,11):
            pick = (self.yaftmar==t) & one_mar & married
            moments['ever kids by years after marriage, {}'.format(t)] = \
                have_kid[self.yaftmar==t].mean() if np.any(pick) else 0.0
         
        for t in range(1,11):
            pick = (self.yaftmar==t) & one_mar
            moments['divorced by years after marriage, {}'.format(t)] = divorced[pick].mean() if np.any(pick) else 0.0
        
        m = self.marriage_stats['Marriages']
        c = self.marriage_stats['Counts']
        haz_m = m['Everyone']/c['Everyone']
        newkids_next_period = (self.planned_preg | self.unplanned_preg)
        childless_this_period = ~have_kid
        
        
        n_newkids = np.zeros(self.Tmax)
        for t in range(n_newkids.size):
            n_newkids[t] = (self.planned_preg | self.unplanned_preg)[age==(21+t)].sum()
        
        n_childless = np.zeros(self.Tmax)
        for t in range(n_newkids.size):
            n_childless[t] = (~have_kid)[age==(21+t)].sum()
            
        print(n_newkids,n_childless)
    
        for t in range(1,16):
            moments['hazard of marriage at {}'.format(21+t)] = haz_m[t] if not np.isnan(haz_m[t]) else 0.0
            
        for t in range(1,16):
            pick = (age==(21+t))
            moments['hazard of new child at {}'.format(21+t)] = n_newkids[t-1]/n_childless[t-1]  if np.any(childless_this_period[pick]) else 0.0
            
            
        k_m_observed = self.k_m & one_mar
        m_k_observed = self.m_k & one_mar
        in_sample = (k_m_observed | m_k_observed)
        
        for t in range(1,16):
            pick = (age==(21+t))
            moments['k then m in population at {}'.format(21+t)] = k_m_observed[pick].mean()
            
        for t in range(1,16):
            pick = (age==(21+t))
            moments['m then k in population at {}'.format(21+t)] = m_k_observed[pick].mean()
            
        for t in range(1,16):
            pick = (age==(21+t)) & in_sample  
            moments['k then m in sample at {}'.format(21+t)] = k_m_observed[pick].mean() if np.any(pick) else 0.0
            
        ls_max = self.labor_supply.max()
        
        moments['mean x share'] = np.mean(self.x[couple_kid]/self.couple_earnings[couple_kid])
        
        pick = one_mar & (age==30)
        moments['divorced at 30 if one marriage'] = divorced[pick].mean() if np.any(pick) else 0.0
        moments['divorced if k then m and one marriage'] = divorced[k_m_observed].mean() if np.any(k_m_observed) else 0.0
        moments['divorced if m then k and one marriage'] = divorced[m_k_observed].mean() if np.any(m_k_observed) else 0.0
        pick = (age==30)
        moments['divorced with kids at 30'] = (divorced & have_kid)[pick].mean() if np.any(pick) else 0.0
        moments['never married with kids at 30'] = ((self.nmar==0) & have_kid)[pick].mean() if     np.any(pick) else 0.0
        pick = (age==40)
        moments['more than one mar at 40'] = (self.nmar>1)[pick].mean() if np.any(pick) else 0.0
        pick = couple_kid & (age==30)
        moments['in labor force at 30 if kids'] = (self.labor_supply > ls_max-1e-4)[pick].mean() if np.any(pick) else 0.0
        
        return moments
            