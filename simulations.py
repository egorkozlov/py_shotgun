#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for simulations
"""

import numpy as np


from mc_tools import mc_simulate, int_prob
from gridvec import VecOnGrid

class Agents:
    
    def __init__(self,Mlist,pswitchlist=None,female=True,N=15000,T=30,
                 pmeet_exo=None,ppreg_exo=None,
                 no_sm=False,verbose=True,nosim=False,fix_seed=True):
            
            
        if fix_seed: np.random.seed(18)
  
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
        
        
        
        
        
        self.female = female
        self.single_state = 'Female, single' if female else 'Male, single'
        
        
        self.pmeet_exo = pmeet_exo
        self.ppreg_exo = ppreg_exo
        
        # all the randomness is here
        self._shocks_single_iexo = np.random.random_sample((N,T))
        self._shocks_single_meet = np.random.random_sample((N,T))
        self._shocks_couple_iexo = np.random.random_sample((N,T))
        self._shocks_single_a = np.random.random_sample((N,T))
        self._shocks_couple_a = np.random.random_sample((N,T))
        
        
        
        
        if female:
            inc_prob = int_prob(self.setup.exogrid.zf_t[0], sig = self.setup.pars['sig_zf_0'])
        else:
            inc_prob = int_prob(self.setup.exogrid.zm_t[0], sig = self.setup.pars['sig_zm_0'])
            
        _shocks_init = np.random.random_sample((N,))        
        i_inc = np.sum((_shocks_init[:,None] > np.cumsum(inc_prob)[None,:]), axis=1)
        iexoinit = i_inc # initial state
        
        
        
        self._shocks_outsm = np.random.random_sample((N,T))
        self._shocks_transition = np.random.random_sample((N,T)) 
        self._shocks_single_preg = np.random.random_sample((N,T))
        self._shocks_planned_preg = np.random.random_sample((N,T))
        
        self._shocks_child_support_fem = np.random.random_sample((N,T))
        self._shocks_child_support_mal = np.random.random_sample((N,T))
        self._shocks_child_support_award = np.random.random_sample((N,T))
        
        self._shocks_init_sm = np.random.random_sample((N,))

        # no randomnes past this line please
        
        # initialize assets
        
        self.iassets = np.zeros((N,T),np.int32)
        
        # initialize FLS
        #self.ils=np.ones((N,T),np.float64)
        self.ils_i=np.zeros((N,T),np.int8)#*(len(self.setup.ls_levels)-1)
        
        
        self.ils_i[:,-1] = 5

        # initialize theta
        self.itheta = -np.ones((N,T),np.int16)
        
        # initialize iexo
        self.iexo = np.zeros((N,T),np.int16)
        
        self.c = np.zeros((N,T),np.float32)
        self.x = np.zeros((N,T),np.float32)
        self.s = np.zeros((N,T),np.float32)
        
        self.unplanned_preg = np.zeros((N,T),dtype=np.bool)
        self.planned_preg = np.zeros((N,T),dtype=np.bool)
        self.disagreed = np.zeros((N,T),dtype=np.bool)
        self.aborted = np.zeros((N,T),dtype=np.bool)
        self.met_a_partner = np.zeros((N,T),dtype=np.bool)
        self.agreed = np.zeros((N,T),dtype=np.bool)   
        self.agreed_k = np.zeros((N,T),dtype=np.bool) 
        self.agreed_unplanned = np.zeros((N,T),dtype=np.bool) 
        self.renegotiated = np.zeros((N,T),dtype=np.bool)
        self.just_divorced = np.zeros((N,T),dtype=np.bool)
        self.just_divorced_nk = np.zeros((N,T),dtype=np.bool)
        self.just_divorced_k = np.zeros((N,T),dtype=np.bool)
        self.new_child = np.zeros((N,T),dtype=np.bool)
        self.k_m = np.zeros((N,T),dtype=np.bool)
        self.k_m_true = np.zeros((N,T),dtype=np.bool)
        self.m_k = np.zeros((N,T),dtype=np.bool)
        self.has_step = np.zeros((N,T),dtype=np.bool)
        
        self.nmar = np.zeros((N,T),dtype=np.int8)
        self.n_kept = np.zeros((T,),dtype=np.int16)
        self.n_aborted = np.zeros((T,),dtype=np.int16)
        self.share_aborted = np.zeros((T,),dtype=np.float64)
        
        self.ub_hit_single = False
        self.ub_hit_couple = False
        
        self.yaftmar = -np.ones((N,T),dtype=np.int8)
        
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
        
        
        
        self.state = np.zeros((N,T),dtype=np.int8)       
        self.state[:,0] = self.state_codes[self.single_state]  
        
        # initial single mothers:
        i_sm_init = (self._shocks_init_sm <= (self.setup.pars['sm_init'] if not no_sm else 0.0))
        self.state[i_sm_init,0] = (self.state_codes['Female and child'] \
                                        if self.female else \
                                        self.state_codes['Male, single'])
        
                                        
                                                
        
        self.timer('Simulations, creation',verbose=self.verbose)
        self.ils_def = 0  # self.setup.nls - 1
            
        
        
        # replace fls to optimal values
        fls_policy = self.Mlist[0].V[0]['Female and child']['fls'] 
        self.ils_i[i_sm_init,0] = fls_policy[self.iassets[i_sm_init,0],iexoinit[i_sm_init]]
        
        
        
        #Create a file with the age of the change foreach person
        
        self.policy_ind = np.zeros((N,T),dtype=np.int8)
        
        if pswitchlist == None:
            pswitchlist = [np.eye(self.npol)]*T
            
        # this simulates "exogenous" transitions of polciy functions
        # policy_ind stands for index of the policies to apply, they are
        # from 0 to (self.npol-1)
        zeros = np.zeros((N,),dtype=np.int8)
        mat_init = pswitchlist[0]
        
        self.policy_ind[:,0] = mc_simulate(zeros,mat_init,shocks=self._shocks_transition[:,0]) # everyone starts with 0
        if self.npol > 1:
            for t in range(T-1):    
                mat = pswitchlist[t+1]
                self.policy_ind[:,t+1] = mc_simulate(self.policy_ind[:,t],mat,shocks=self._shocks_transition[:,t+1])
        else:
            self.policy_ind[:] = 0
            
            
        if not nosim: self.run_sim()
            
            
    def run_sim(self):        
        self.simulate()
        self.compute_aux()
        counts, offers, marriages = self.marriage_stats()
        self.compute_cse()
    
        if self.verbose and self.ub_hit_single: print('Assets upped bound is reached for singles')
        if self.verbose and self.ub_hit_couple: print('Assets upped bound is reached for couples')
        
    
    def tht_interpolate(self,v,i_other,i_theta):
        it = self.setup.v_thetagrid_fine.i[i_theta]
        itp = self.setup.v_thetagrid_fine.i[i_theta]
        wt = self.setup.v_thetagrid_fine.wthis[i_theta]
        wtp = self.setup.v_thetagrid_fine.wnext[i_theta]
        
        ind_t = i_other + (it,)
        ind_tp = i_other + (itp,)
        return v[ind_t]*wt + v[ind_tp]*wtp
        

           
        
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
                
                #pol = self.Mlist[ipol].decisions[t][sname]
                vpol = self.Mlist[ipol].V[t][sname]
                    
                
                if not use_theta:
                    
                    # apply for singles
                    
                    
                    anext = vpol['s'][self.iassets[ind,t],self.iexo[ind,t]]
                    self.iassets[ind,t+1] = VecOnGrid(self.setup.agrid_s,anext).roll(shocks=self._shocks_single_a[ind,t])
                    self.s[ind,t] = anext
                    self.c[ind,t] = vpol['c'][self.iassets[ind,t],self.iexo[ind,t]]
                    self.x[ind,t] = vpol['x'][self.iassets[ind,t],self.iexo[ind,t]]
                    if np.any(self.iassets[ind,t+1]==self.setup.na-1): self.ub_hit_single=True
                    
                else:
                    
                    # interpolate in both assets and theta
                    # function apply_2dim is experimental but I checked it at this setup
                    
                    # apply for couples
                    
                    def tint(x): return self.tht_interpolate(x,(self.iassets[ind,t],self.iexo[ind,t]),self.itheta[ind,t])
                                    
                    anext = tint(vpol['s'])
                    self.s[ind,t] = anext
                    
                    
                    self.x[ind,t] = tint(vpol['x'])
                    self.c[ind,t] = tint(vpol['c'])
                    
                    self.iassets[ind,t+1] = VecOnGrid(self.setup.agrid_c,anext).roll(shocks=self._shocks_couple_a[ind,t])
                    if np.any(self.iassets[ind,t+1]==self.setup.na-1): self.couple=True
                    
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
                    
                    nls = self.setup.nls[sname]
                            
                    lvls = self.setup.ls_levels[sname]
                    
                    ls_val = self.ils_i[ind,t] 
                    
                    for ils in range(nls):
                        this_ls = (ls_val==ils)                    
                        if not np.any(this_ls): continue
                    
                        cnt = np.sum(this_ls)
                        lvl = lvls[ils]
        
                        
                        if self.verbose: print('At t = {} for {} {} have LS of {}'.format(t,sname,cnt,lvl))
                        
                        
                        mat = self.Mlist[ipol].setup.exo_mats[sname][ils][t]
                        
                        shks = self._shocks_couple_iexo[ind[this_ls],t]
                        
                        iexo_next_this_ls = mc_simulate(iexo_now[this_ls],mat,shocks=shks)
                        self.iexo[ind[this_ls],t+1] = iexo_next_this_ls
                        
                else:
                    mat = self.Mlist[ipol].setup.exo_mats[sname][t]
                    shks = self._shocks_single_iexo[ind,t]                    
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
                
                
                
                if sname == "Female, single" or sname == "Male, single" or sname == "Female and child":
                    # TODO: this is temporary version, it computes partners for
                    # everyone and after that imposes meet / no meet, this should
                    # not happen.
                    
                    # meet a partner
                    
                    
                    ss = self.single_state
                    ss_sm = 'Female and child' if self.female else self.single_state 
                    
                    pcoef = self.Mlist[ipol].setup.pars['pmeet_multiplier_fem']
                    
                    
                    if self.pmeet_exo is None:
                        pmeet = pcoef*self.Mlist[ipol].setup.pars['pmeet_t'][t] # TODO: check timing
                    else:
                        try:
                            pmeet = pcoef*self.pmeet_exo[t]
                        except IndexError:
                            pmeet = pcoef*self.pmeet_exo[-1]
                    
                    p_abortion_access = self.Mlist[ipol].setup.pars['p_abortion_access']
                    
                    # divide by 2 subgroups
                    
                    matches = self.Mlist[ipol].decisions[t][sname]
                    
                    
                    ia = self.iassets[ind,t+1] # note that timing is slightly inconsistent  
                    
                    # we use iexo from t and savings from t+1
                    # TODO: fix the seed
                    iznow = self.iexo[ind,t]
                    
                    
                    if not 'Not pregnant' in matches or (pmeet<1e-5):
                        try:
                            assert pmeet<1e-5
                        except:
                            print(t)
                            print(pmeet)
                            print(sname)
                            print(matches.keys())
                            assert False
                        # propagate yaftmar
                        # !!! this needs to be verified
                        
                        # iexo is set within iexonext
                        izf = self.iexo[ind,t+1]
                        self.yaftmar[ind,t+1] = \
                            (self.yaftmar[ind,t] + 1)*(self.yaftmar[ind,t]>=0)+(-1)*(self.yaftmar[ind,t]<0)
                        
                        if sname != "Female and child":
                            self.ils_i[ind,t+1] = self.ils_def
                            self.state[ind,t+1] = self.state_codes[ss]
                        else:
                            fls_policy = self.Mlist[ipol].V[t+1]['Female and child']['fls']
                            self.state[ind,t+1] = self.state_codes['Female and child']  
                            self.ils_i[ind,t+1] = fls_policy[self.iassets[ind,t+1],izf]
                        continue
                        
                        
                        
                        
                    
                    pmat_np = matches['Not pregnant']['p_mat_extended'][iznow,:] # assuming it is identical !!!
                    pmat_cum_np = pmat_np.cumsum(axis=1)
                    pmat_cum_p = pmat_cum_np
                
                    assert np.all(pmat_cum_np < 1 + 1e-5) 
                    assert np.all(pmat_cum_np[:,-1] > 1 - 1e-5)
                    
                    assert np.all(pmat_cum_p < 1 + 1e-5) 
                    assert np.all(pmat_cum_p[:,-1] > 1 - 1e-5)
                    assert np.all(pmat_cum_p >= 0.0)
                    # there is a rare numerical issue when adding lots of floats
                    # gives imprecise result > 1
                    
                    pmat_cum_np[:,-1] = 1+1e-5 
                    pmat_cum_p[:,-1] = 1+1e-5 
                    
                    
                    
                    
                    v = self._shocks_single_iexo[ind,t] #np.random.random_sample(ind.size) # draw uniform dist
                    
                    i_pmat_np = (v[:,None] > pmat_cum_np).sum(axis=1)  # index of the position in pmat
                    i_pmat_p = (v[:,None] > pmat_cum_p).sum(axis=1)  # index of the position in pmat
                    
                    
                    # these things are the same for pregnant and not pregnant
                    
                    
                    if self.female:
                        fert = self.setup.pars['is fertile'][t]
                    else:
                        fert = self.setup.pars['is fertile'][t-2] if t>=2 else False
                    
                    
                    try:
                        mult = self.setup.pars['ppreg_sim_mult']
                    except:
                        mult = 1.0
                    
                    if sname == 'Female, single':
                        if self.ppreg_exo is None:
                            p_preg = mult*fert*self.setup.upp_precomputed_fem[t][self.iexo[ind,t]]
                        else:
                            p_preg = mult*fert*self.ppreg_exo[t]
                            
                    elif sname == 'Male, single':
                        if self.ppreg_exo is None:
                            p_preg = mult*fert*self.setup.upp_precomputed_mal[t][self.iexo[ind,t]]
                        else:
                            p_preg = mult*fert*self.ppreg_exo[t]
                            #except IndexError:
                            #    p_preg = fert*self.ppreg_exo[-1]
                    elif sname == 'Female and child':
                        p_preg = 1.0 #np.ones_like(ind,dtype=np.float64)[]
                    else:
                        assert False
                        
                        
                    
                    
                    #print('p_preg is {}'.format(p_preg))
                    # these are individual-specific pregnancy probabilities
                    # for those who are not fertile this is forced to be zero
                        
                    # potential assets position of couple
                    
                    vpreg = self._shocks_single_preg[ind,t]
                    i_preg = (vpreg <= p_preg)
                    
                    
                    #i_pmat = i_pmat_p*(i_preg) + i_pmat_np*(~i_preg)
                    
                    
                    ic_out = matches['Not pregnant']['corresponding_iexo'][i_pmat_np]*(~i_preg) \
                                + matches['Pregnant']['corresponding_iexo'][i_pmat_p]*(i_preg)
                    
                    imatch = matches['Not pregnant']['corresponding_imatch'][i_pmat_np]*(~i_preg) \
                                + matches['Pregnant']['corresponding_imatch'][i_pmat_p]*(i_preg)
                    
                     
                    ia_out = matches['Not pregnant']['ia_c_table'][ia,i_pmat_np]*(~i_preg) \
                                + matches['Pregnant']['ia_c_table'][ia,i_pmat_p]*(i_preg)
                                
                    it_out = matches['Not pregnant']['itheta'][ia,i_pmat_np]*(~i_preg) + \
                                + matches['Pregnant']['itheta'][ia,i_pmat_p]*(i_preg)
                    
                    
                    iall, izf, izm, ipsi = self.Mlist[ipol].setup.all_indices(t,ic_out)
                    
                    
                    # compute for everyone
                    
                    
                    vmeet = self._shocks_single_meet[ind,t]
                    i_nomeet =  np.array( vmeet > pmeet )
                    
                    
                    v_abortion = self._shocks_outsm[ind,t]
                    i_abortion_allowed = np.array( v_abortion < p_abortion_access )
                    
                    
                    
                    i_pot_agree = matches['Not pregnant']['Decision'][ia,i_pmat_np]*(~i_preg) \
                                    + matches['Pregnant']['Decision'][ia,i_pmat_p]*(i_preg)
                                    
                    i_m_preferred = matches['Not pregnant']['Child immediately'][ia,i_pmat_np]*(~i_preg) \
                                    + matches['Pregnant']['Child immediately'][ia,i_pmat_p]*(i_preg)
                                    
                    i_abortion_preferred = matches['Not pregnant']['Abortion'][ia,i_pmat_np]*(~i_preg) \
                                    + matches['Pregnant']['Abortion'][ia,i_pmat_p]*(i_preg)
                                    
                    
                    # TODO: child immediately is a binary decision and does not
                    # respect infertility.
                    
                    
                    i_disagree = (~i_pot_agree)
                    
                    
                    
                    i_disagree_or_nomeet = (i_disagree) | (i_nomeet)
                    i_disagree_and_meet = (i_disagree) & ~(i_nomeet)
                    
                    i_abortion = (i_abortion_allowed) & (i_abortion_preferred) & (i_disagree_and_meet) & (i_preg) &  (sname=='Female, single')
                    i_kept_sm = (i_abortion_allowed) & (~i_abortion_preferred) & (i_disagree_and_meet) &  (i_preg) & (sname=='Female, single')
                    i_no_access_sm = ~(i_abortion_allowed) & (i_disagree_and_meet) &  (i_preg) &  (sname=='Female, single')
                    
                    
                    
                    
                    
                    
                    i_agree = ~i_disagree_or_nomeet
    
                    
                    i_agree_mar = (i_agree) & (i_m_preferred)
                    i_agree_coh = (i_agree) & (~i_m_preferred)
                    
                    
                    i_kept_kf = (i_agree_mar & ~(sname=='Female and child'))
                    
                    
                    self.disagreed[ind,t] = i_disagree_and_meet
                    self.met_a_partner[ind,t] = ~i_nomeet
                    self.unplanned_preg[ind,t] = (i_preg) & (~i_nomeet) & ~(sname=='Female and child')
                    self.aborted[ind,t] = i_abortion
                    
                    n_abortions = i_abortion.sum()
                    n_kept = i_kept_sm.sum() + i_no_access_sm.sum() + i_kept_kf.sum()
                    
                    
                    
                    
                    
                    if sname == 'Female, single':
                        self.share_aborted[t] = 100*n_abortions / (n_abortions + n_kept)
                        self.n_kept[t] = n_kept
                        self.n_aborted[t] = n_abortions
                    
                    if n_abortions>0 and self.verbose: print('{} abortions done at t = {} for {}'.format(n_abortions,t,sname))
                    if n_kept>0 and self.verbose: print('{} abortions refused at t = {} for {}'.format(n_kept,t,sname))
                    
                    
                    if not sname=='Female and child':
                        become_sm = (i_disagree_and_meet) & (i_preg) & self.female & ~(i_abortion)
                        become_single = (i_disagree_or_nomeet) & ~(become_sm)   
                        #assert not np.any(become_sm)
                    else:
                        become_sm = (i_disagree_or_nomeet) & (i_preg)
                        assert np.all(i_preg)
                        become_single = np.zeros_like(become_sm,dtype=np.bool)
                    assert np.all(i_disagree_or_nomeet == ((become_sm) | (become_single)))
                    
                    
                    self.agreed[ind,t] = (i_agree_mar) | (i_agree_coh)
                    self.planned_preg[ind,t] = (i_agree_mar) & ~(i_preg)
                    
                    
                    self.new_child[ind,t+1] = (i_kept_sm | i_no_access_sm) | (i_agree_mar & ~(sname=='Female and child'))
                    
                    
                    
                    assert np.all(~i_nomeet[i_agree])
                    
                    i_firstmar = (self.nmar[ind,t]==0)
                    
                    nmar, ncoh, ndis, nnom = np.sum(i_agree_mar),np.sum(i_agree_coh),np.sum(i_disagree_and_meet),np.sum(i_nomeet)
                    ntot = sum((nmar, ncoh, ndis, nnom))
                    
                    if self.verbose: print('for sname = {}: {} mar, {} coh,  {} disagreed, {} did not meet ({} total)'.format(sname,nmar,ncoh,ndis,nnom,ntot))
                    
                    
                    if np.any(i_agree_mar):
                        
                        self.itheta[ind[i_agree_mar],t+1] = it_out[i_agree_mar]
                        self.iexo[ind[i_agree_mar],t+1] = iall[i_agree_mar]
                        self.state[ind[i_agree_mar],t+1] = self.state_codes['Couple and child']
                        self.iassets[ind[i_agree_mar],t+1] = ia_out[i_agree_mar]
                        
                        self.agreed_k[ind[i_agree_mar],t] = True
                        self.agreed_unplanned[ind[i_agree_mar],t] = i_preg[i_agree_mar]*(sname!='Female and child')
                        
                        self.k_m[ind[i_agree_mar],(t+1):] = True
                        self.k_m_true[ind[i_agree_mar],(t+1):] = (i_preg[i_agree_mar][:,None] & (sname!='Female and child'))
                        
                        if sname=='Female and child':
                            self.has_step[ind[i_agree_mar],(t+1):] = True
                        
                        fls_policy = self.Mlist[ipol].V[t+1]['Couple and child']['fls']
                        
                        
                        def thti(*agrs): return np.round(self.tht_interpolate(*agrs)).astype(np.int8)
                        
                        self.ils_i[ind[i_agree_mar],t+1] = \
                            thti(fls_policy,(self.iassets[ind[i_agree_mar],t+1],self.iexo[ind[i_agree_mar],t+1]),self.itheta[ind[i_agree_mar],t+1])
                        
                        self.yaftmar[ind[i_agree_mar],t+1] = 0                        
                        self.nmar[ind[i_agree_mar],t+1:] += 1
                        
                    if np.any(i_agree_coh):
                        
                        assert not sname=='Female and child'
                        
                        self.itheta[ind[i_agree_coh],t+1] = it_out[i_agree_coh]
                        self.iexo[ind[i_agree_coh],t+1] = iall[i_agree_coh]
                        self.state[ind[i_agree_coh],t+1] = self.state_codes['Couple, no children']
                        self.iassets[ind[i_agree_coh],t+1] = ia_out[i_agree_coh]
                        
                        # FLS decision
                        #tg = self.Mlist[ipol].setup.v_thetagrid_fine
                        #fls_policy = self.V[t+1]['Couple, no children']['fls']
                        
                        def thti(*agrs): return np.round(self.tht_interpolate(*agrs)).astype(np.int8)
                        
                        
                        
                        fls_policy = self.Mlist[ipol].V[t+1]['Couple, no children']['fls']
                        
                        self.ils_i[ind[i_agree_coh],t+1] = \
                            thti(fls_policy,(self.iassets[ind[i_agree_coh],t+1],self.iexo[ind[i_agree_coh],t+1]),self.itheta[ind[i_agree_coh],t+1])
                        
                        self.yaftmar[ind[i_agree_coh],t+1] = 0                                                
                        self.nmar[ind[i_agree_coh],t+1:] = self.nmar[ind[i_agree_coh],t][:,None] + 1
                        
                        
                        
                    # determine izf and izm in case of disagreement
                    
                    p_csa = self.setup.pars['child_support_awarded_nm']
                    izf_cs, izm_cs = self.disagreement_with_child_support(t,ind,izf,izm,p_csa)
                    
                       
                        
                    if np.any(i_disagree_or_nomeet):
                        # do not touch assets
                        fls_policy = self.Mlist[ipol].V[t+1]['Female and child']['fls']
                        
                        
                        if sname=='Female and child': assert not np.any(become_single)
                        if sname=='Male, single': assert not np.any(become_sm)
                        
                        
                        cs_eligible = (become_sm & (sname != 'Female and child'))
                        
                        izf_no = izf*(~cs_eligible) + izf_cs*(cs_eligible)
                        izm_no = izm*(~cs_eligible) + izm_cs*(cs_eligible)
                        
                        iz = izf_no if self.female else izm_no
                        self.iexo[ind[i_disagree_or_nomeet],t+1] = iz[i_disagree_or_nomeet]
                        #self.state[ind[i_disagree_or_nomeet],t+1] = self.state_codes['Female, single']
                        self.state[ind[become_single],t+1] = self.state_codes[ss]
                        self.state[ind[become_sm],t+1] = self.state_codes['Female and child']
                        self.ils_i[ind[become_single],t+1] = self.ils_def
                        self.ils_i[ind[become_sm],t+1] = fls_policy[self.iassets[ind[become_sm],t+1],
                                                                       iz[become_sm]]
                        
                        self.yaftmar[ind[i_disagree_or_nomeet],t+1] = \
                            (self.yaftmar[ind[i_disagree_or_nomeet],t] + 1)*\
                                 (self.yaftmar[ind[i_disagree_or_nomeet],t]>=0)+\
                            (-1)*(self.yaftmar[ind[i_disagree_or_nomeet],t]<0)
                        
                
                
                    #print('next period kf: {}'.format(self.k_m[:,t+1].mean()))
                    #assert p_preg > 0.0
                    
                elif sname == "Couple and child" or sname == "Couple, no children":
                    
                    
                    ss = 'Female, single' if self.female else 'Male, single'
                    
                    decision = self.Mlist[ipol].decisions[t][sname]
    
                    haschild = (sname == "Couple and child")
                    
                    # by default keep the same theta and weights
                    
                    self.itheta[ind,t+1] = self.itheta[ind,t]
                    
                    nt = self.Mlist[ipol].setup.ntheta_fine
                    
                    
                    # initiate renegotiation
                    isc = self.iassets[ind,t+1]
                    iall, izf, izm, ipsi = self.Mlist[ipol].setup.all_indices(t+1,self.iexo[ind,t+1])
                    
                    if sname == "Couple and child":
                        
                        p_csa = self.setup.pars['child_support_awarded_div']
                        izf_div, izm_div = self.disagreement_with_child_support(t,ind,izf,izm,p_csa)
                        
                        
                    elif sname == "Couple, no children":
                        izf_div = izf
                        izm_div = izm
                    else:
                        assert False
                    
                    
                    itht = self.itheta[ind,t+1] 
                    agrid =  self.Mlist[ipol].setup.agrid_c 
                    agrid_s =  self.Mlist[ipol].setup.agrid_s
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
                    
                    assert np.all(self.yaftmar[ind,t]>=0)
                    self.yaftmar[ind,t+1] = self.yaftmar[ind,t] + 1 # extra year past marriage
                    
                    
                    pfert = self.setup.pars['fert_prob_t'][t]
                    
                    if np.any(i_div):
                        
                        
                        income_fem = np.exp(zf_grid[izf[i_div]])
                        income_mal = np.exp(zm_grid[izm[i_div]])
                        
                        income_share_fem = income_fem / (income_fem + income_mal)
                        
                        # !!!
                        costs = self.Mlist[ipol].setup.divorce_costs_k if sname == 'Couple and child' else self.Mlist[ipol].setup.divorce_costs_nk
                                   
                        share_f, share_m = costs.shares_if_split(income_share_fem)
                        
                        sf = share_f*sc[i_div]
                        sm = share_m*sc[i_div]
                        
                        self.just_divorced[ind,t] = i_div
                        
                        if haschild:
                            self.just_divorced_k[ind,t] = i_div
                        else:
                            self.just_divorced_nk[ind,t] = i_div
                        
                        shks = self._shocks_couple_a[ind[i_div],t]
                        
                        # FIXME: it should be agrid_s here
                        if self.female:
                            self.iassets[ind[i_div],t+1] = VecOnGrid(agrid_s,sf).roll(shocks=shks)
                        else:
                            self.iassets[ind[i_div],t+1] = VecOnGrid(agrid,sm).roll(shocks=shks)
                            
                        self.itheta[ind[i_div],t+1] = -1
                        iz = izf_div if self.female else izm_div
                        self.iexo[ind[i_div],t+1] = iz[i_div]
                        fls_policy = self.Mlist[ipol].V[t+1]['Female and child']['fls']
                        
                        if haschild and self.female:
                            self.state[ind[i_div],t+1] = self.state_codes['Female and child']
                            self.ils_i[ind[i_div],t+1] = fls_policy[self.iassets[ind[i_div],t+1],
                                                                       izf[i_div]]
                        else:
                            self.state[ind[i_div],t+1] = self.state_codes[ss]
                            self.ils_i[ind[i_div],t+1] = self.ils_def
                        
                        
                        
                        
                    if np.any(i_ren):
                        
                        self.itheta[ind[i_ren],t+1] = thts[i_ren]
                        
                        
                        #tg = self.setup.v_thetagrid_fine
                        
                        #Distinguish between marriage and cohabitation
                        if sname == "Couple and child":
                            self.state[ind[i_ren],t+1] = self.state_codes[sname]
                            
                            
                            ipick = (self.iassets[ind[i_ren],t+1],self.iexo[ind[i_ren],t+1],self.itheta[ind[i_ren],t+1])
                            def thti(*agrs): return np.round(self.tht_interpolate(*agrs)).astype(np.int8)
                            self.ils_i[ind[i_ren],t+1] = thti(self.Mlist[ipol].V[t+1][sname]['fls'],ipick[:-1],ipick[-1])
                        else:
                            
                            
                            i_birth = (decision['Give a birth'][isc,iall,thts]*pfert > self._shocks_planned_preg[ind,t])
                            i_birth1=i_birth[i_ren]
                            
                            self.planned_preg[ind[i_ren],t] = i_birth1
                            
                            self.m_k[ind[i_ren][i_birth1],(t+1):] = True
                            
                            self.new_child[ind[i_ren],t+1] = i_birth1
                                           
                            ipick = (self.iassets[ind[i_ren],t+1],self.iexo[ind[i_ren],t+1],self.itheta[ind[i_ren],t+1])
                            
                            def thti(*agrs): return np.round(self.tht_interpolate(*agrs)).astype(np.int8)

                            
                            ils_if_k = thti(self.Mlist[ipol].V[t+1]["Couple and child"]['fls'],ipick[:-1],ipick[-1])
                            ils_if_nk = thti(self.Mlist[ipol].V[t+1]["Couple, no children"]['fls'],ipick[:-1],ipick[-1])
                            
                            self.ils_i[ind[i_ren],t+1] = i_birth1*ils_if_k+(1-i_birth1)*ils_if_nk
                            self.state[ind[i_ren],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                            
                                
                           
                        
                    if np.any(i_sq):
                        self.state[ind[i_sq],t+1] = self.state_codes[sname]
                        # do not touch theta as already updated
                        
                        def thti(*agrs): return np.round(self.tht_interpolate(*agrs)).astype(np.int8)
                        
                        if sname == "Couple and child":
                            self.state[ind[i_sq],t+1] = self.state_codes[sname]
                            
                            ipick = (self.iassets[ind[i_sq],t+1],self.iexo[ind[i_sq],t+1],self.itheta[ind[i_sq],t+1])
                            self.ils_i[ind[i_sq],t+1] = thti(self.Mlist[ipol].V[t+1][sname]['fls'],ipick[:-1],ipick[-1])
                        else:
                            i_birth = (decision['Give a birth'][isc,iall,thts]*pfert > self._shocks_planned_preg[ind,t])
                            
                            
                            i_birth1=i_birth[i_sq]
                            self.m_k[ind[i_sq][i_birth1],(t+1):] = True                            
                            self.planned_preg[ind[i_sq],t] = i_birth1
                            self.state[ind[i_sq],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                            
                            self.new_child[ind[i_sq],t+1] = i_birth1
                            
                            ipick = (self.iassets[ind[i_sq],t+1],self.iexo[ind[i_sq],t+1],self.itheta[ind[i_sq],t+1])
                            
                            ils_if_k = thti(self.Mlist[ipol].V[t+1]["Couple and child"]['fls'],ipick[:-1],ipick[-1])
                            ils_if_nk = thti(self.Mlist[ipol].V[t+1]["Couple, no children"]['fls'],ipick[:-1],ipick[-1])
                           
                            self.ils_i[ind[i_sq],t+1] = i_birth1*ils_if_k+(1-i_birth1)*ils_if_nk
                            self.state[ind[i_sq],t+1] = i_birth1*self.state_codes["Couple and child"]+(1-i_birth1)*self.state_codes["Couple, no children"]
                
                else:
                    raise Exception('unsupported state?')
        
        assert not np.any(np.isnan(self.state[:,t+1]))   
        
        
    def disagreement_with_child_support(self,t,ind,izf,izm,p_award):
        
        transitions = self.setup.child_support_transitions[t+1]
        izf_div_0 = transitions['i_this_fem'][izf,izm]
        izf_div_1 = izf_div_0 + 1
        wzf_div_0 = transitions['w_this_fem'][izf,izm]
        pick_zf_0 = (self._shocks_child_support_fem[ind,t] < wzf_div_0)                        
        izf_div_cs = izf_div_0*(pick_zf_0) + izf_div_1*(~pick_zf_0)
        
        izm_div_0 = transitions['i_this_mal'][izm]
        izm_div_1 = izm_div_0 + 1
        wzm_div_0 = transitions['w_this_mal'][izm]
        pick_zm_0 = (self._shocks_child_support_mal[ind,t] < wzm_div_0)                        
        izm_div_cs = izm_div_0*(pick_zm_0) + izm_div_1*(~pick_zm_0)
        
        award_cs = (self._shocks_child_support_award[ind,t] < p_award) 
        
        izf_div = izf_div_cs*(award_cs) + izf*(~award_cs)
        izm_div = izm_div_cs*(award_cs) + izm*(~award_cs)
        
        if self.setup.pars['child_support_share'] < 1e-5 or p_award<1e-4:
            assert np.all(izf_div==izf)
            assert np.all(izm_div==izm)
        return izf_div, izm_div
        
        
        
    def compute_aux(self):
        # this should be ran after the last simulations
        
        
        
        couple_k = (self.state == self.state_codes['Couple and child'])
        single_k = (self.state == self.state_codes['Female and child'])
        self.labor_supply = np.ones((self.N,self.T),dtype=np.float32)
        self.psi_couple = -1000*np.ones((self.N,self.T),dtype=np.float32)
        self.theta_couple = -1*np.ones((self.N,self.T),dtype=np.float32)
        
        
        self.female_wage = -np.ones((self.N,self.T),dtype=np.float32)
        self.female_z = -1000*np.ones((self.N,self.T),dtype=np.float32)
        self.male_z = -1000*np.ones((self.N,self.T),dtype=np.float32)
        self.male_wage = -np.ones((self.N,self.T),dtype=np.float32)
        self.female_earnings = -np.ones((self.N,self.T),dtype=np.float32)
        self.male_earnings = -np.ones((self.N,self.T),dtype=np.float32)
        self.couple_earnings = -np.ones((self.N,self.T),dtype=np.float32)
        self.total_earnings = -np.ones((self.N,self.T),dtype=np.float32)
        self.total_resources = -np.ones((self.N,self.T),dtype=np.float32)
        self.taxes_paid = -1000*np.ones((self.N,self.T),dtype=np.float32)
        self.tax_rate = -1*np.ones((self.N,self.T),dtype=np.float32)
        
        
        self.total_expenditures = self.c + self.x + self.s
        
        
        self.savings_to_earnings = -np.ones((self.N,self.T),dtype=np.float32)
        
        self.labor_supply[couple_k] = self.setup.ls_levels['Couple and child'][self.ils_i[couple_k]]
        self.labor_supply[single_k] = self.setup.ls_levels['Female and child'][self.ils_i[single_k]]
        
        
        
        
        assert np.all(self.x>=0)
        
        
        # fill wage & earnings
        
        max_savings = -np.ones(len(self.state_codes))
        
        for t in range(self.T):
            tm = self.setup.pars['m_wage_trend'][t]
            tf = self.setup.pars['f_wage_trend'][t]
            R = self.setup.pars['R_t'][t]
            agrid_s = self.setup.agrid_s
            agrid_c = self.setup.agrid_c
            
            
            
            for state, i_state in self.state_codes.items():
                pick = (self.state[:,t] == i_state)
                
            
                
                if not np.any(pick): continue
                
                
                iexo = self.iexo[pick,t] # reshaped
                
                max_savings[i_state] = max(max_savings[i_state],np.max(self.s[pick,t]))
                
                if state == 'Female, single' or state == 'Female and child':
                    wage = np.exp(self.setup.exogrid.zf_t[t][iexo] + tf)
                    self.female_wage[pick,t] = wage
                    self.female_z[pick,t] = self.setup.exogrid.zf_t[t][iexo]
                    self.female_earnings[pick,t] = wage*self.labor_supply[pick,t]
                    self.savings_to_earnings[pick,t] = self.s[pick,t] / self.female_earnings[pick,t]
                    self.total_earnings[pick,t] = self.female_earnings[pick,t]
                    self.total_resources[pick,t] = self.total_earnings[pick,t] + R*agrid_s[self.iassets[pick,t]]
                    assert np.all(wage>0)
                elif state == 'Male, single':                    
                    wage = np.exp(self.setup.exogrid.zm_t[t][iexo] + tm)
                    self.male_wage[pick,t] = wage
                    self.male_z[pick,t] = self.setup.exogrid.zm_t[t][iexo]
                    self.male_earnings[pick,t] = wage # !!!
                    self.savings_to_earnings[pick,t] = self.s[pick,t] / self.male_earnings[pick,t]
                    self.total_earnings[pick,t] = self.male_earnings[pick,t]
                    self.total_resources[pick,t] = self.total_earnings[pick,t] + R*agrid_s[self.iassets[pick,t]]
                    assert np.all(wage>0)
                elif state == 'Couple and child' or state == 'Couple, no children':
                    iall, izf, izm, ipsi = self.setup.all_indices(t,iexo)
                    wage_f = np.exp(self.setup.exogrid.zf_t[t][izf] + tf)
                    wage_m = np.exp(self.setup.exogrid.zm_t[t][izm] + tm)
                    self.female_wage[pick,t] = wage_f
                    self.male_wage[pick,t] = wage_m
                    self.male_z[pick,t] = self.setup.exogrid.zm_t[t][izm]
                    self.female_z[pick,t] = self.setup.exogrid.zf_t[t][izf]                    
                    self.female_earnings[pick,t] = wage_f*self.labor_supply[pick,t]
                    self.male_earnings[pick,t] = wage_m
                    self.couple_earnings[pick,t] = self.female_earnings[pick,t] + self.male_earnings[pick,t]
                    self.savings_to_earnings[pick,t] = self.s[pick,t] / self.couple_earnings[pick,t]
                    self.theta_couple[pick,t] = self.setup.thetagrid_fine[self.itheta[pick,t]]
                    self.psi_couple[pick,t] = self.setup.exogrid.psi_t[t][ipsi]
                    
                    self.total_earnings[pick,t] = self.female_earnings[pick,t] + self.male_earnings[pick,t]
                    self.total_resources[pick,t] = self.total_earnings[pick,t] + R*agrid_c[self.iassets[pick,t]]
                    
                    
                    assert np.all(wage_f>0)
                    assert np.all(wage_m>0)
                else:
                    raise Exception('unsupported state code?')
                
                self.taxes_paid[pick,t] = self.total_resources[pick,t] - self.total_expenditures[pick,t]
                self.tax_rate[pick,t] = self.taxes_paid[pick,t] / self.total_earnings[pick,t]
                
        
        for state, i_state in self.state_codes.items():
            if max_savings[i_state] < 0: continue
            if self.verbose: print('Maximum savings for {} is {}'.\
                                   format(state,max_savings[i_state]))
            
        
    def marriage_stats(self):
        
        
        offer_single = np.zeros((self.T,),dtype=np.int32)
        offer_preg = np.zeros((self.T,),dtype=np.int32)
        offer_notpreg = np.zeros((self.T,),dtype=np.int32)
        offer_sm = np.zeros((self.T,),dtype=np.int32)
        
        total_single = np.zeros((self.T,),dtype=np.int32)
        total_sm = np.zeros((self.T,),dtype=np.int32)
        
        nmar_single = np.zeros((self.T,),dtype=np.int32)
        nmar_preg = np.zeros((self.T,),dtype=np.int32)
        nmar_notpreg = np.zeros((self.T,),dtype=np.int32)
        nmar_sm = np.zeros((self.T,),dtype=np.int32)
        
        
        
        for t in range(self.T):
            
            pick = (self.state[:,t] == self.state_codes[self.single_state])
            total_single[t] = np.sum(pick)
            offer_single[t] = self.met_a_partner[:,t][pick].sum() if np.any(pick) else 0
            offer_preg[t] = self.unplanned_preg[:,t][pick].sum() if np.any(pick) else 0
            offer_notpreg[t] = offer_single[t] - offer_preg[t]
            
            nmar_single[t] = self.agreed[:,t][pick].sum() if np.any(pick) else 0
            nmar_preg[t] = self.agreed_k[:,t][pick].sum() if np.any(pick) else 0
            nmar_notpreg[t] = nmar_single[t] - nmar_preg[t]
            
            
            
            
            
            pick = (self.state[:,t] == self.state_codes['Female and child'])
            total_sm[t] = np.sum(pick)
            offer_sm[t] = self.met_a_partner[:,t][pick].sum() if np.any(pick) else 0
            nmar_sm[t] = self.agreed[:,t][pick].sum() if np.any(pick) else 0
            
            
            
        total_all = total_single + total_sm
        offer_all = offer_single + offer_sm
        nmar_all = nmar_single + nmar_sm
        
        counts = {'Everyone':total_all,'Single':total_single,'SM':total_sm}
        
        offers = {'Everyone':offer_all,'Single':offer_single,'SM':offer_sm,
                 'Single, pregnant':offer_preg,
                 'Single, not pregnant':offer_notpreg}
        
        marriages = {'Everyone':nmar_all,'Single':nmar_single,'SM':nmar_sm,
                 'Single, pregnant':nmar_preg,
                 'Single, not pregnant':nmar_notpreg,'With kids':nmar_sm+nmar_preg}
            
        return counts, offers, marriages
    
    
    def sim_graphs(self):
        from sim_graphs import plot_sim_graphs
        plot_sim_graphs(self)
        
        
    def compute_cse(self):
        # this child consumption equivalent measure: 
        # average value given the same utility as Q of the child for
        # the first 18 years of life.
        
        bet = self.setup.pars['beta_t'][0] # assumes constant beta
        weights = np.array([bet**t for t in range(18)])
        coefs = np.cumsum(weights)
        coef_18 = coefs[-1]
        
        
        cse = np.zeros(self.N,dtype=np.float64)
        self.cse = np.zeros((self.N,self.T),dtype=np.float64)
        
        xi  = self.setup.pars['util_xi']
        lam = self.setup.pars['util_lam']
        kap = self.setup.pars['util_kap']  
        
        
        self.q = -1*np.ones((self.N,self.T),dtype=np.float64)
        ip = (self.x>0)
        self.q[ip] = (self.x[ip]**lam + kap*(1-self.labor_supply[ip])**lam)**(1/lam)
        uq = np.zeros_like(self.q)
        uq[~ip] = None
        uq[ip] = self.q[ip]**(1-xi)/(1-xi)
        
        
        for t in range(self.T-1):
            
            periods_ahead = (self.T-2) - t
            coef = coef_18 if periods_ahead >= 17 else coefs[periods_ahead]
            
            ind_ahead = min(18,periods_ahead+1)
            
            if t == 0:
                i_born = (self.x[:,t] > 0)
            else:
                i_born = (self.x[:,t] > 0) & (self.x[:,t-1] == 0)
            
            if not np.any(i_born):
                continue
            #else:
                #print((t,periods_ahead))
        
            assert np.all(cse[i_born] == 0)
            
            uqpick = uq[i_born,t:(t+ind_ahead)]
            
            wpick = weights[:ind_ahead][None,:]
            
            assert np.all(~np.isnan(uqpick))
            
            
            cse[i_born] = ((np.sum(uqpick*wpick,axis=1)*(1-xi)/coef)**(1/(1-xi)))
            
            self.cse[i_born,t:] = cse[i_born][:,None]
        
            
            
            
            
        
        
        
        
            
    def compute_moments(self):
        import moments
        return moments.compute_moments(self)
    
    def aux_moments(self):
        import moments
        return moments.aux_moments(self)
            
            
        
        
