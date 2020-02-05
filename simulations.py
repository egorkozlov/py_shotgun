#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for simulations
"""

import numpy as np


from mc_tools import mc_simulate, int_prob
from gridvec import VecOnGrid

class Agents:
    
    def __init__(self,Mlist,pswitchlist=None,N=15000,T=None,verbose=True,nosim=False):
            
            
        np.random.seed(18)
  
        # take the stuff from the model and arguments
        # note that this does not induce any copying just creates links
        
        if type(Mlist) is not list:
            Mlist = [Mlist]
            

        
        #Unilateral Divorce
        self.Mlist = Mlist
        self.Vlist = [M.V for M in Mlist]
        self.declist = [M.decisions for M in Mlist]
        self.npol = len(Mlist)
        self.transition = len(self.Mlist)>1
        
            
        if T is None:
            T = self.Mlist[0].setup.pars['Tsim']
            
        self.setup = self.Mlist[0].setup
        self.state_names = self.setup.state_names
        self.N = N
        self.T = T
        self.verbose = verbose
        self.timer = self.Mlist[0].time
        
        # all the randomness is here
        self.shocks_single_iexo = np.random.random_sample((N,T))
        self.shocks_single_meet = np.random.random_sample((N,T))
        self.shocks_couple_iexo = np.random.random_sample((N,T))
        self.shocks_single_a = np.random.random_sample((N,T))
        self.shocks_couple_a = np.random.random_sample((N,T))
        
        
        fem_prob = int_prob(self.setup.exogrid.zf_t[0], sig = self.setup.pars['sig_zf_0'])
        shocks_init = np.random.random_sample((N,))        
        i_fem = np.sum((shocks_init[:,None] > np.cumsum(fem_prob)[None,:]), axis=1)
        iexoinit = i_fem # initial state
        
        
        
        self.shocks_outsm = np.random.random_sample((N,T))
        self.shocks_transition = np.random.random_sample((N,T)) 
        self.shocks_single_preg = np.random.random_sample((N,T))
        # no randomnes past this line please
        
        # initialize assets
        
        self.iassets = np.zeros((N,T),np.int32)
        
        # initialize FLS
        #self.ils=np.ones((N,T),np.float64)
        self.ils_i=np.zeros((N,T),np.int32)#*(len(self.setup.ls_levels)-1)
        
        
        self.ils_i[:,-1] = 5

        # initialize theta
        self.itheta = -np.ones((N,T),np.int32)
        
        # initialize iexo
        self.iexo = np.zeros((N,T),np.int32)
        
        self.c = np.zeros((N,T),np.float32)
        self.x = np.zeros((N,T),np.float32)
        self.s = np.zeros((N,T),np.float32)
        
        self.unplanned_preg = np.zeros((N,T),dtype=np.bool)
        self.planned_preg = np.zeros((N,T),dtype=np.bool)
        self.disagreed = np.zeros((N,T),dtype=np.bool)
        self.met_a_partner = np.zeros((N,T),dtype=np.bool)
        self.agreed = np.zeros((N,T),dtype=np.bool)   
        self.agreed_k = np.zeros((N,T),dtype=np.bool) 
        self.agreed_unplanned = np.zeros((N,T),dtype=np.bool) 
        self.renegotiated = np.zeros((N,T),dtype=np.bool)
        self.just_divorced = np.zeros((N,T),dtype=np.bool)
        self.just_divorced_nk = np.zeros((N,T),dtype=np.bool)
        self.just_divorced_k = np.zeros((N,T),dtype=np.bool)
        self.k_m = np.zeros((N,T),dtype=np.bool)
        self.k_m_true = np.zeros((N,T),dtype=np.bool)
        self.m_k = np.zeros((N,T),dtype=np.bool)
        
        self.iexo[:,0] = iexoinit
        
        
        
        # initialize state
        
        
        
        self.state_codes = dict()
        self.has_theta = list()
        self.has_fls = list()
        for i, name in enumerate(self.setup.state_names):
            self.state_codes[name] = i
            self.has_theta.append((name=='Couple, no children' or name=='Couple and child'))
            self.has_fls.append((name=='Couple, no children' or name=='Couple and child' or
                                     name=='Female and child'))
        
        
        
        self.state = np.zeros((N,T),dtype=np.int32)       
        self.state[:,0] = self.state_codes['Female, single']  # everyone starts as female
        
        self.timer('Simulations, creation',verbose=self.verbose)
        self.ils_def = 0#self.setup.nls - 1
            
        #Create a file with the age of the change foreach person
        
        self.policy_ind = np.zeros((N,T),dtype=np.int8)
        
        if pswitchlist == None:
            pswitchlist = [np.eye(self.npol)]*T
            
        # this simulates "exogenous" transitions of polciy functions
        # policy_ind stands for index of the policies to apply, they are
        # from 0 to (self.npol-1)
        zeros = np.zeros((N,),dtype=np.int8)
        mat_init = pswitchlist[0]
        
        self.policy_ind[:,0] = mc_simulate(zeros,mat_init,shocks=self.shocks_transition[:,0]) # everyone starts with 0
        if self.npol > 1:
            for t in range(T-1):    
                mat = pswitchlist[t+1]
                self.policy_ind[:,t+1] = mc_simulate(self.policy_ind[:,t],mat,shocks=self.shocks_transition[:,t+1])
        else:
            self.policy_ind[:] = 0
            
        if not nosim: self.simulate()
            
        
    def simulate(self):
        
        #Create Variables that stores varibles of interest
        
        
        for t in range(self.T-1):
         
            self.anext(t) 
            self.iexonext(t)            
            self.statenext(t)
            self.timer('Simulations, iteration',verbose=False)
        
        
        #return self.gsavings, self.iexo, self.state,self.gtheta
    
    
    
    def anext(self,t):
        # finds savings (potenitally off-grid)
        
        
        for ipol in range(self.npol):
            for ist, sname in enumerate(self.state_codes):
                
                is_state_any_pol = (self.state[:,t]==ist)  
                is_pol = (self.policy_ind[:,t]==ipol)
                
                is_state = (is_state_any_pol) & (is_pol)
                
                use_theta = self.has_theta[ist]            
                nst = np.sum(is_state)
                
                if nst==0:
                    continue
                
                ind = np.where(is_state)[0]
                
                pol = self.Mlist[ipol].decisions[t][sname]
                    
                
                if not use_theta:
                    
                    # apply for singles
                    
                    
                    
                    anext = pol['s'][self.iassets[ind,t],self.iexo[ind,t]]
                    self.iassets[ind,t+1] = VecOnGrid(self.setup.agrid_s,anext).roll(shocks=self.shocks_single_a[ind,t])
                    self.s[ind,t] = anext
                    self.c[ind,t] = pol['c'][self.iassets[ind,t],self.iexo[ind,t]]
                    self.x[ind,t] = pol['x'][self.iassets[ind,t],self.iexo[ind,t]]
                else:
                    
                    # interpolate in both assets and theta
                    # function apply_2dim is experimental but I checked it at this setup
                    
                    # apply for couples
                                        
                    anext = pol['s'][self.iassets[ind,t],self.iexo[ind,t],self.itheta[ind,t]]
                    self.s[ind,t] = anext
                    self.x[ind,t] = pol['x'][self.iassets[ind,t],self.iexo[ind,t],self.itheta[ind,t]]
                    self.c[ind,t] = pol['c'][self.iassets[ind,t],self.iexo[ind,t],self.itheta[ind,t]]
                    
                    self.iassets[ind,t+1] = VecOnGrid(self.setup.agrid_c,anext).roll(shocks=self.shocks_couple_a[ind,t])
                    
                assert np.all(anext >= 0)
    
    
    def iexonext(self,t):
        
        # let's find out new exogenous state
        
        for ipol in range(self.npol):
            for ist,sname in enumerate(self.state_names):
                is_state_any_pol = (self.state[:,t]==ist)
                is_pol = (self.policy_ind[:,t]==ipol)
                is_state = (is_state_any_pol) & (is_pol)   
                
                nst = np.sum(is_state)
                
                if nst == 0:
                    continue
                
                ind = np.where(is_state)[0]
                sname = self.state_names[ist]
                iexo_now = self.iexo[ind,t].reshape(nst)
                
                
                
                
                if sname == 'Couple, no children' or sname == 'Couple and child' or sname == 'Female and child':
                    
                    nls = self.setup.nls_k if (sname == 'Couple and child')  \
                        else self.setup.nls_nk if (sname == 'Couple, no children') \
                            else self.setup.nls_sk 
                            
                    lvls = self.setup.ls_levels_k if (sname == 'Couple and child') \
                        else self.setup.ls_levels_nk if (sname == 'Couple, no children')  \
                            else self.setup.ls_levels_sk
                    
                    ls_val = self.ils_i[ind,t] 
                    
                    for ils in range(nls):
                        this_ls = (ls_val==ils)                    
                        if not np.any(this_ls): continue
                    
                        cnt = np.sum(this_ls)
                        lvl = lvls[ils]
        
                        
                        if self.verbose: print('At t = {} for {} {} have LS of {}'.format(t,sname,cnt,lvl))
                        
                        
                        mat = self.Mlist[ipol].setup.exo_mats[sname][ils][t]
                        
                        shks = self.shocks_couple_iexo[ind[this_ls],t]
                        
                        iexo_next_this_ls = mc_simulate(iexo_now[this_ls],mat,shocks=shks)
                        self.iexo[ind[this_ls],t+1] = iexo_next_this_ls
                        
                else:
                    mat = self.Mlist[ipol].setup.exo_mats[sname][t]
                    shks = self.shocks_single_iexo[ind,t]                    
                    iexo_next = mc_simulate(iexo_now,mat,shocks=shks) # import + add shocks     
                    self.iexo[ind,t+1] = iexo_next
            
        
    def statenext(self,t):
        
        
        
        
        for ipol in range(self.npol):
            for ist,sname in enumerate(self.state_names):
                is_state_any_pol = (self.state[:,t]==ist)
                is_pol = (self.policy_ind[:,t]==ipol)
                is_state = (is_state_any_pol) & (is_pol) 
                
                
                if not np.any(is_state):
                    continue
                
                
                if self.verbose: print('At t = {} count of {} is {}'.format(t+1,sname,np.sum(is_state)))
                
                ind = np.where(is_state)[0]
                
                nind = ind.size
                
                
                
                if sname == "Female, single":
                    # TODO: this is temporary version, it computes partners for
                    # everyone and after that imposes meet / no meet, this should
                    # not happen.
                    
                    # meet a partner
                    
                    pmeet = self.Mlist[ipol].setup.pars['pmeet_t'][t] # TODO: check timing
                    
                    
                    # divide by 2 subgroups
                    
                    matches = self.Mlist[ipol].decisions[t]['Female, single']
                    
                    
                    ia = self.iassets[ind,t+1] # note that timing is slightly inconsistent  
                    
                    # we use iexo from t and savings from t+1
                    # TODO: fix the seed
                    iznow = self.iexo[ind,t]
                    
                    pmat = matches['Not pregnant']['p'][ia,iznow,:]
                    
                    pmat_cum = pmat.cumsum(axis=1)
                    
                    assert np.all(pmat_cum < 1 + 1e-5) 
                    assert np.all(pmat_cum[:,-1] > 1 - 1e-5)
                    assert np.all(pmat_cum >= 0.0)
                    # there is a rare numerical issue when adding lots of floats
                    # gives imprecise result > 1
                    
                    pmat_cum[:,-1] = 1+1e-5 
                    
                    
                    
                    
                    v = self.shocks_single_iexo[ind,t] #np.random.random_sample(ind.size) # draw uniform dist
                    
                    i_pmat = (v[:,None] > pmat_cum).sum(axis=1)  # index of the position in pmat
                    
                    
                    # these things are the same for pregnant and not pregnant
                    
                    p_preg = matches['Pregnant']['Probability Unplanned'][ia,iznow,i_pmat]
                    # these are individual-specific pregnancy probabilities
                    # for those who are not fertile this is forced to be zero
                        
                    # potential assets position of couple
                    
                    vpreg = self.shocks_single_preg[ind,t]
                    i_preg = (vpreg < p_preg)
                    
                    
                    
                    ic_out = matches['Not pregnant']['iexo'][ia,iznow,i_pmat]*(~i_preg) \
                                + matches['Pregnant']['iexo'][ia,iznow,i_pmat]*(i_preg)
                                
                    ia_out = matches['Not pregnant']['ia'][ia,iznow,i_pmat]*(~i_preg) \
                                + matches['Pregnant']['ia'][ia,iznow,i_pmat]*(i_preg)
                                
                    it_out = matches['Not pregnant']['theta'][ia,iznow,i_pmat]*(~i_preg) + \
                                + matches['Pregnant']['theta'][ia,iznow,i_pmat]*(i_preg)
                    
                    
                    iall, izf, izm, ipsi = self.Mlist[ipol].setup.all_indices(t,ic_out)
                    
                    
                    # compute for everyone
                    
                    
                    vmeet = self.shocks_single_meet[ind,t]
                    i_nomeet =  np.array( vmeet > pmeet )
                    
                    self.met_a_partner[ind,t] = ~i_nomeet
                    self.unplanned_preg[ind,t] = (i_preg) & (~i_nomeet)
                    
                    i_pot_agree = matches['Not pregnant']['Decision'][ia,iznow,i_pmat]*(~i_preg) \
                                    + matches['Pregnant']['Decision'][ia,iznow,i_pmat]*(i_preg)
                                    
                    i_m_preferred = matches['Not pregnant']['Child immediately'][ia,iznow,i_pmat]*(~i_preg) \
                                    + matches['Pregnant']['Child immediately'][ia,iznow,i_pmat]*(i_preg)
                    
                    i_disagree = (~i_pot_agree)
                    
                    
                    i_disagree_or_nomeet = (i_disagree) | (i_nomeet)
                    i_disagree_and_meet = (i_disagree) & ~(i_nomeet)
                    
                    self.disagreed[ind,t] = i_disagree_and_meet
                    
                    become_sm = (i_disagree_and_meet) & (i_preg)
                    become_single = (i_disagree_or_nomeet) & ~(become_sm)   
                    assert np.all(i_disagree_or_nomeet == ((become_sm) | (become_single)))
                    
                    i_agree = ~i_disagree_or_nomeet
    
                    
                    i_agree_mar = (i_agree) & (i_m_preferred)
                    i_agree_coh = (i_agree) & (~i_m_preferred)
                    
                    self.agreed[ind,t] = (i_agree_mar) | (i_agree_coh)
                    self.planned_preg[ind,t] = (i_agree_mar) & ~(i_preg)
                    
                    assert np.all(~i_nomeet[i_agree])
                    
                    
                    nmar, ncoh, ndis, nnom = np.sum(i_agree_mar),np.sum(i_agree_coh),np.sum(i_disagree_and_meet),np.sum(i_nomeet)
                    ntot = sum((nmar, ncoh, ndis, nnom))
                    
                    #if self.verbose: print('{} mar, {} coh,  {} disagreed, {} did not meet ({} total)'.format(nmar,ncoh,ndis,nnom,ntot))
                    #assert np.all(ismar==(i_agree )
                    
                    if np.any(i_agree_mar):
                        
                        self.itheta[ind[i_agree_mar],t+1] = it_out[i_agree_mar]
                        self.iexo[ind[i_agree_mar],t+1] = iall[i_agree_mar]
                        self.state[ind[i_agree_mar],t+1] = self.state_codes['Couple and child']
                        self.iassets[ind[i_agree_mar],t+1] = ia_out[i_agree_mar]
                        
                        self.agreed_k[ind[i_agree_mar],t] = True
                        self.agreed_unplanned[ind[i_agree_mar],t] = i_preg[i_agree_mar]
                        
                        self.k_m[ind[i_agree_mar],t:] = True
                        self.k_m_true[ind[i_agree_mar],t:] = i_preg[i_agree_mar][:,None]
                        
                        
                        # FLS decision
                        #self.ils_i[ind[i_ren],t+1] = 
                        tg = self.Mlist[ipol].setup.v_thetagrid_fine                    
                        fls_policy = self.Mlist[ipol].decisions[t+1]['Couple and child']['fls']
                        
                        self.ils_i[ind[i_agree_mar],t+1] = \
                            fls_policy[self.iassets[ind[i_agree_mar],t+1],self.iexo[ind[i_agree_mar],t+1],self.itheta[ind[i_agree_mar],t+1]]
                        
                        
                    if np.any(i_agree_coh):
                        
                        self.itheta[ind[i_agree_coh],t+1] = it_out[i_agree_coh]
                        self.iexo[ind[i_agree_coh],t+1] = iall[i_agree_coh]
                        self.state[ind[i_agree_coh],t+1] = self.state_codes['Couple, no children']
                        self.iassets[ind[i_agree_coh],t+1] = ia_out[i_agree_coh]
                        
                        # FLS decision
                        tg = self.Mlist[ipol].setup.v_thetagrid_fine
                        #fls_policy = self.V[t+1]['Couple, no children']['fls']
                        fls_policy = self.Mlist[ipol].decisions[t+1]['Couple, no children']['fls']
                        
                        self.ils_i[ind[i_agree_coh],t+1] = \
                            fls_policy[self.iassets[ind[i_agree_coh],t+1],self.iexo[ind[i_agree_coh],t+1],self.itheta[ind[i_agree_coh],t+1]]
                        
                    
                        
                    if np.any(i_disagree_or_nomeet):
                        # do not touch assets
                        fls_policy = self.Mlist[ipol].decisions[t+1]['Female and child']['fls']
                        
                        self.iexo[ind[i_disagree_or_nomeet],t+1] = izf[i_disagree_or_nomeet]
                        #self.state[ind[i_disagree_or_nomeet],t+1] = self.state_codes['Female, single']
                        self.state[ind[become_single],t+1] = self.state_codes['Female, single']
                        self.state[ind[become_sm],t+1] = self.state_codes['Female and child']
                        self.ils_i[ind[become_single],t+1] = self.ils_def
                        self.ils_i[ind[become_sm],t+1] = fls_policy[self.iassets[ind[become_sm],t+1],
                                                                       izf[become_sm]]
                        
                        
                elif sname == "Female and child":
                    
                    pout = self.Mlist[ipol].setup.pars['poutsm_t'][t]
                    i_out_sm = (self.shocks_outsm[ind,t]<pout)
                    self.state[ind[ i_out_sm],t+1] = self.state_codes['Female, single']
                    self.state[ind[~i_out_sm],t+1] = self.state_codes['Female and child']
                    
                    fls_policy = self.Mlist[ipol].decisions[t+1]['Female and child']['fls']
                    
                    self.ils_i[ind[ i_out_sm],t+1] = self.ils_def
                    self.ils_i[ind[~i_out_sm],t+1] = fls_policy[self.iassets[ind[~i_out_sm],t+1],
                                                                       self.iexo[ind[~i_out_sm],t+1]]
                    
                
                
                elif sname == "Couple and child" or sname == "Couple, no children":
                    
                    decision = self.Mlist[ipol].decisions[t][sname]
    
                    haschild = (sname == "Couple and child")
                    
                    # by default keep the same theta and weights
                    
                    self.itheta[ind,t+1] = self.itheta[ind,t]
                    
                    nt = self.Mlist[ipol].setup.ntheta_fine
                    
                    
                    # initiate renegotiation
                    isc = self.iassets[ind,t+1]
                    iall, izf, izm, ipsi = self.Mlist[ipol].setup.all_indices(t+1,self.iexo[ind,t+1])
                    
                    itht = self.itheta[ind,t+1] 
                    agrid =  self.Mlist[ipol].setup.agrid_c                
                    sc = agrid[isc] # needed only for dividing asssets               
                    
                    thts_all = decision['thetas']
                    thts_orig_all = np.broadcast_to(np.arange(nt)[None,None,:],thts_all.shape)
                    
                    
                    thts = thts_all[isc,iall,itht]
                    thts_orig = thts_orig_all[isc,iall,itht]
                    
                    dec = decision['Decision']
                    
                    i_stay = dec[isc,iall] if dec.ndim==2 else dec[isc,iall,itht]
    
                    
                    
                    i_div = ~i_stay    
                    
    
                    i_ren = (i_stay) & (thts_orig != thts)
                    i_renf = (i_stay) & (thts_orig > thts)
                    i_renm = (i_stay) & (thts_orig < thts)
                    i_sq = (i_stay) & (thts_orig == thts)
                        
                    
                    #if self.verbose: print('{} divorce, {} ren-f, {} ren-m, {} sq'.format(np.sum(i_div),np.sum(i_renf),np.sum(i_renm),np.sum(i_sq))                     )
                    
                    self.renegotiated[ind,t] = i_ren
                    
                    zf_grid = self.setup.exo_grids['Female, single'][t]
                    zm_grid = self.setup.exo_grids['Male, single'][t]
                    
                    
                     
                    
                    if np.any(i_div):
                        
                        
                        income_fem = np.exp(zf_grid[izf[i_div]])
                        income_mal = np.exp(zm_grid[izm[i_div]])
                        
                        income_share_fem = income_fem / (income_fem + income_mal)
                        
                        # !!!
                        costs = self.Mlist[ipol].setup.divorce_costs_k if sname == 'Couple and child' else self.Mlist[ipol].setup.divorce_costs_nk
                                   
                        share_f, share_m = costs.shares_if_split(income_share_fem)
                        
                        sf = share_f*sc[i_div]
                        
                        self.just_divorced[ind,t] = i_div
                        if haschild:
                            self.just_divorced_k[ind,t] = i_div
                        else:
                            self.just_divorced_nk[ind,t] = i_div
                        
                        shks = self.shocks_couple_a[ind[i_div],t]
                        self.iassets[ind[i_div],t+1] = VecOnGrid(agrid,sf).roll(shocks=shks)
                        self.itheta[ind[i_div],t+1] = -1
                        self.iexo[ind[i_div],t+1] = izf[i_div]
                        fls_policy = self.Mlist[ipol].decisions[t+1]['Female and child']['fls']
                        
                        if haschild:
                            self.state[ind[i_div],t+1] = self.state_codes['Female and child']
                            self.ils_i[ind[i_div],t+1] = fls_policy[self.iassets[ind[i_div],t+1],
                                                                       izf[i_div]]
                        else:
                            self.state[ind[i_div],t+1] = self.state_codes['Female, single']
                            self.ils_i[ind[i_div],t+1] = self.ils_def
                        
                        
                        
                        
                    if np.any(i_ren):
                        
                        self.itheta[ind[i_ren],t+1] = thts[i_ren]
                        
                        
                        #tg = self.setup.v_thetagrid_fine
                        
                        #Distinguish between marriage and cohabitation
                        if sname == "Couple and child":
                            self.state[ind[i_ren],t+1] = self.state_codes[sname]
                            
                            
                            ipick = (self.iassets[ind[i_ren],t+1],self.iexo[ind[i_ren],t+1],self.itheta[ind[i_ren],t+1])
                            self.ils_i[ind[i_ren],t+1] = self.Mlist[ipol].decisions[t+1][sname]['fls'][ipick]
                        else:
                            i_birth = decision['Give a birth'][isc,iall,thts]
                            i_birth1=i_birth[i_ren]
                            
                            self.planned_preg[ind[i_ren],t] = i_birth1
                            
                            self.m_k[ind[i_ren][i_birth1],t:] = True
                            
                            
                            ipick = (self.iassets[ind[i_ren],t+1],self.iexo[ind[i_ren],t+1],self.itheta[ind[i_ren],t+1])
                            ils_if_k = self.Mlist[ipol].decisions[t+1]["Couple and child"]['fls'][ipick]
                            ils_if_nk = self.Mlist[ipol].decisions[t+1]["Couple, no children"]['fls'][ipick]
                            
                            self.ils_i[ind[i_ren],t+1] = i_birth1*ils_if_k+(1-i_birth1)*ils_if_nk
                            self.state[ind[i_ren],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                          
                                
                           
                        
                    if np.any(i_sq):
                        self.state[ind[i_sq],t+1] = self.state_codes[sname]
                        # do not touch theta as already updated
                        
                        if sname == "Couple and child":
                            self.state[ind[i_sq],t+1] = self.state_codes[sname]
                            
                            ipick = (self.iassets[ind[i_sq],t+1],self.iexo[ind[i_sq],t+1],self.itheta[ind[i_sq],t+1])
                            self.ils_i[ind[i_sq],t+1] = self.Mlist[ipol].decisions[t+1][sname]['fls'][ipick]
                        else:
                            i_birth = decision['Give a birth'][isc,iall,thts]
                            i_birth1=i_birth[i_sq]
                            self.m_k[ind[i_sq][i_birth1],t:] = True                            
                            self.planned_preg[ind[i_sq],t] = i_birth1
                            self.state[ind[i_sq],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                            
                            ipick = (self.iassets[ind[i_sq],t+1],self.iexo[ind[i_sq],t+1],self.itheta[ind[i_sq],t+1])
                            
                            ils_if_k = self.Mlist[ipol].decisions[t+1]["Couple and child"]['fls'][ipick]
                            ils_if_nk = self.Mlist[ipol].decisions[t+1]["Couple, no children"]['fls'][ipick]
                           
                            self.ils_i[ind[i_sq],t+1] = i_birth1*ils_if_k+(1-i_birth1)*ils_if_nk
                            self.state[ind[i_sq],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                
                else:
                    raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))            