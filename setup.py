#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for setting up the model
"""

import numpy as np

from rw_approximations import rouw_nonst
from mc_tools import combine_matrices_two_lists, int_prob,cut_matrix
from scipy.stats import norm
from collections import namedtuple
from gridvec import VecOnGrid

from scipy import sparse



class ModelSetup(object):
    def __init__(self,nogrid=False,divorce_costs_k='Default',divorce_costs_nk='Default',**kwargs): 
        p = dict()       
        T = 55
        Tret = 40 # first period when the agent is retired
        Tfert = 15 # first peroid when infertile
        Tdiv = 35 # first period when cannot divorce / renegotiate
        Tmeet = 25 # first period when stop meeting partners
        p['T'] = T
        p['Tret'] = Tret
        p['Tfert'] = Tfert
        p['Tsim'] = T
        p['sig_zf_0']  = 0.513 #0.25
        p['sig_zf']    = 0.084#np.sqrt(0.038)
        p['n_zf_t']      = [7]*Tret + [1]*(T-Tret)
        p['sig_zm_0']  = 0.525 #0.25
        p['sig_zm']    = 0.097 #np.sqrt(0.027)
        p['n_zm_t']      = [5]*Tret + [1]*(T-Tret)
        p['sigma_psi_init'] = 0.28
        p['sigma_psi']   = 0.11
        p['R_t'] = [1.04]*T
        p['n_psi_t']     = [12]*T
        p['beta_t'] = [0.98]*T
        p['A'] = 1.0 # consumption in couple: c = (1/A)*[c_f^(1+rho) + c_m^(1+rho)]^(1/(1+rho))
        p['crra_power'] = 1.5
        p['couple_rts'] = 0.23    
        p['sig_partner_a'] = 0.1
        p['sig_partner_z'] = 0.4
        p['m_bargaining_weight'] = 0.5
        p['pmeet'] = 0.5
        p['pmeet_t'] = 0.0
        
        p['poutsm'] = 1/3
        p['z_drift'] = -0.1
        p['no kids at meeting'] = True
        
        
        p['wret'] = 0.8
        p['uls'] = 0.2
        p['pls'] = 1.0
        
        
        p['u_shift_mar'] = 0.0
        p['u_shift_coh'] = 0.0
        
        
        
        '''
        p['f_wage_trend'] = np.array(
                            [0.0 + 
                             0.068386*(min(t,30) - 5)
                             - 0.0040512*((min(t,30)-5)**2)
                             + 0.0000784*((min(t,30)-5)**3)
                                         for t in range(T)]
                                    )
        
        p['m_wage_trend'] = np.array(
                                    [0.1457256 + 0.0676344*(min(t,30)-5) - 
                                         0.002752*((min(t,30)-5)**2) 
                                         + 0.0000784*((min(t,30)-5)**3)
                                         for t in range(T)]
                                    )
        '''
        '''
        se_05_m = 0.5596
        se_20_m = 0.6735
        
        sigma_m = np.sqrt( (se_20_m**2 - se_05_m**2)/15 )
        sigma_m0 = np.sqrt(se_05_m**2 - 4*sigma_m**2)
        
        se_05_f = 0.5397
        se_20_f = 0.6305
        
        sigma_f = np.sqrt( (se_20_f**2 - se_05_f**2)/15 )
        sigma_f0 = np.sqrt(se_05_f**2 - 4*sigma_f**2)
        
        '''
        
        
        
        
        p['f_wage_trend'] = np.array(
                            [0.0 + 
                             0.0844278*(min(t,30) - 5)
                             -0.0044622*((min(t,30)-5)**2)
                             + 0.0000788*((min(t,30)-5)**3)
                                         for t in range(T)]
                                    )
        
        p['m_wage_trend'] = np.array(
                                    [0.0776162 + 0.0810113*(min(t,30)-5)  
                                         - 0.0027093*((min(t,30)-5)**2) 
                                         + 0.0000298*((min(t,30)-5)**3)
                                         for t in range(T)]
                                    )
        
        
        
        p['util_lam'] = 0.4
        p['util_alp'] = 0.5
        p['util_xi'] = 1.5
        p['util_kap'] = 0.5
        p['util_qbar'] = 0.0
        
        
        p['preg_20'] = 0.2
        p['preg_28'] = 0.1
        
        
        for key, value in kwargs.items():
            assert (key in p), 'wrong name?'
            p[key] = value
        
        
        
        dp = (p['preg_28'] - p['preg_20'])/8
         
        p['preg_a0'] = p['preg_20']
        p['preg_at'] = dp
        p['preg_az'] = 0.00
        p['preg_azt'] = 0.00
        
        #Get the probability of meeting, adjusting for year-period
           
        
        p['is fertile'] = [True]*Tfert + [False]*(T-Tfert)
        p['can divorce'] = [True]*Tdiv + [False]*(T-Tdiv)        
        p['pmeet_t'] = [np.clip(p['pmeet'] + t*p['pmeet_t'],0.0,1.0) for t in range(Tmeet)] + [0.0]*(T-Tmeet)
        p['poutsm_t'] = [p['poutsm']]*T
        
        
        self.pars = p
        
        self.dtype = np.float32 # type for all floats
        
       
        
        # relevant for integration
        self.state_names = ['Female, single','Male, single','Female and child','Couple, no children','Couple and child']
        
        # female labor supply
        
        lmin = 0.2
        lmax = 1.0
        nl = 4
        
        ls = np.linspace(lmin,lmax,nl,dtype=self.dtype)
        ps = p['pls']*(1-np.linspace(0.0,2.0,nl,dtype=self.dtype))
        ps[-2:] = 0.0 # !!! think about this carefully
        
        self.ls_levels_nk = np.array([1.0],dtype=self.dtype)
        self.ls_levels_k = ls
        self.ls_levels_sk = ls
        
        #self.ls_utilities = np.array([p['uls'],0.0],dtype=self.dtype)
        self.ls_pdown_nk = np.array([0.0],dtype=self.dtype)
        self.ls_pdown_k = ps
        self.ls_pdown_sk = ps
        self.nls_k = len(self.ls_levels_k)
        self.nls_nk = len(self.ls_levels_nk)
        self.nls_sk = len(self.ls_levels_sk)
        
        
        
        
        
        #Cost of Divorce
        if divorce_costs_k == 'Default':
            # by default the costs are set in the bottom
            self.divorce_costs_k = DivorceCosts()
        else:
            if isinstance(divorce_costs_k,dict):
                # you can feed in arguments to DivorceCosts
                self.divorce_costs_k = DivorceCosts(**divorce_costs_k)
            else:
                # or just the output of DivorceCosts
                assert isinstance(divorce_costs_k,DivorceCosts)
                self.divorce_costs_k = divorce_costs_k
                
        #Cost of Separation
        if divorce_costs_nk == 'Default':
            # by default the costs are set in the bottom
            self.dov_costs_nok = DivorceCosts()
        else:
            if isinstance(divorce_costs_nk,dict):
                # you can feed in arguments to DivorceCosts
                self.divorce_costs_nk = DivorceCosts(**divorce_costs_nk)
            else:
                # or just the output of DivorceCosts
                assert isinstance(divorce_costs_nk,DivorceCosts)
                self.divorce_costs_nk = divorce_costs_nk
            
        # exogrid should be deprecated
        if not nogrid:
        
            exogrid = dict()
            
            
            # let's approximate three Markov chains
            # this sets up exogenous grid
            
            # FIXME: this uses number of points from 0th entry. 
            # in principle we can generalize this
            
            exogrid['zf_t'],  exogrid['zf_t_mat'] = rouw_nonst(p['T'],p['sig_zf'],p['sig_zf_0'],p['n_zf_t'][0])
            exogrid['zm_t'],  exogrid['zm_t_mat'] = rouw_nonst(p['T'],p['sig_zm'],p['sig_zm_0'],p['n_zm_t'][0])
            
            for t in range(Tret,T):
                exogrid['zf_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zm_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zf_t_mat'][t] = np.atleast_2d(1.0)
                exogrid['zm_t_mat'][t] = np.atleast_2d(1.0)
                
                
            # fix transition from non-retired to retired    
            exogrid['zf_t_mat'][Tret-1] = np.ones((p['n_zf_t'][Tret-1],1))
            exogrid['zm_t_mat'][Tret-1] = np.ones((p['n_zm_t'][Tret-1],1))
            
            exogrid['psi_t'], exogrid['psi_t_mat'] = rouw_nonst(p['T'],p['sigma_psi'],p['sigma_psi_init'],p['n_psi_t'][0])
            
            zfzm, zfzmmat = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], exogrid['zf_t_mat'], exogrid['zm_t_mat'])
            all_t, all_t_mat = combine_matrices_two_lists(zfzm,exogrid['psi_t'],zfzmmat,exogrid['psi_t_mat'])
            all_t_mat_sparse_T = [sparse.csc_matrix(D.T) if D is not None else None for D in all_t_mat]
            
            
            #Create a new bad version of transition matrix p(zf_t)
            
            
            
            zf_bad = [tauchen_drift(exogrid['zf_t'][t], exogrid['zf_t'][t+1], 
                                    1.0, p['sig_zf'], p['z_drift'])
                        for t in range(self.pars['T']-1) ] + [None]
            
            
            #zf_bad = [cut_matrix(exogrid['zf_t_mat'][t]) if t < Tret -1 
            #              else (exogrid['zf_t_mat'][t] if t < T - 1 else None) 
            #                  for t in range(self.pars['T'])]
            
            zf_t_mat_down = zf_bad
            
            zfzm, zfzmmat = combine_matrices_two_lists(exogrid['zf_t'], exogrid['zm_t'], zf_t_mat_down, exogrid['zm_t_mat'])
            all_t_down, all_t_mat_down = combine_matrices_two_lists(zfzm,exogrid['psi_t'],zfzmmat,exogrid['psi_t_mat'])
            all_t_mat_down_sparse_T = [sparse.csc_matrix(D.T) if D is not None else None for D in all_t_mat_down]
            
            
            
            all_t_mat_by_l_nk = [ [(1-p)*m + p*md if m is not None else None 
                                for m , md in zip(all_t_mat,all_t_mat_down)]
                               for p in self.ls_pdown_nk ]
            
            all_t_mat_by_l_spt_nk = [ [(1-p)*m + p*md if m is not None else None
                                    for m, md in zip(all_t_mat_sparse_T,all_t_mat_down_sparse_T)]
                               for p in self.ls_pdown_nk ]
            
            all_t_mat_by_l_k = [ [(1-p)*m + p*md if m is not None else None 
                                for m , md in zip(all_t_mat,all_t_mat_down)]
                               for p in self.ls_pdown_k ]
            
            all_t_mat_by_l_spt_k = [ [(1-p)*m + p*md if m is not None else None
                                    for m, md in zip(all_t_mat_sparse_T,all_t_mat_down_sparse_T)]
                               for p in self.ls_pdown_k ]
            
            zf_t_mat_by_l_sk = [ [(1-p)*m + p*md if md is not None else None 
                                for m , md in zip(exogrid['zf_t_mat'],zf_bad)]
                                   for p in self.ls_pdown_sk ]
            
            
            exogrid['all_t_mat_by_l_nk'] = all_t_mat_by_l_nk
            exogrid['all_t_mat_by_l_spt_nk'] = all_t_mat_by_l_spt_nk
            
            exogrid['all_t_mat_by_l_k'] = all_t_mat_by_l_k
            exogrid['all_t_mat_by_l_spt_k'] = all_t_mat_by_l_spt_k
            
            exogrid['zf_t_mat_by_l_sk'] = zf_t_mat_by_l_sk
            
            exogrid['all_t'] = all_t
            
            Exogrid_nt = namedtuple('Exogrid_nt',exogrid.keys())
            
            self.exogrid = Exogrid_nt(**exogrid)
            self.pars['nexo_t'] = [v.shape[0] for v in all_t]
            
            #assert False
            
            
            
        #Grid Couple
        self.na = 40
        self.amin = 0
        self.amax = 60
        self.agrid_c = np.linspace(self.amin,self.amax,self.na,dtype=self.dtype)
        tune=1.5
        #self.agrid_c = np.geomspace(self.amin+tune,self.amax+tune,num=self.na)-tune
        
        # this builds finer grid for potential savings
        s_between = 7 # default numer of points between poitns on agrid
        s_da_min = 0.001 # minimal step (does not create more points)
        s_da_max = 0.1 # maximal step (creates more if not enough)
        
        self.sgrid_c = build_s_grid(self.agrid_c,s_between,s_da_min,s_da_max)
        self.vsgrid_c = VecOnGrid(self.agrid_c,self.sgrid_c)
        
        
         
        #Grid Single
        self.amin_s = 0
        self.amax_s = self.amax/1.1
        self.agrid_s = np.linspace(self.amin_s,self.amax_s,self.na,dtype=self.dtype)
        tune_s=1.5
        #self.agrid_s = np.geomspace(self.amin_s+tune_s,self.amax_s+tune_s,num=self.na)-tune_s
        
        self.sgrid_s = build_s_grid(self.agrid_s,s_between,s_da_min,s_da_max)
        self.vsgrid_s = VecOnGrid(self.agrid_s,self.sgrid_s)
        
        # grid for theta
        self.ntheta = 11
        self.thetamin = 0.01
        self.thetamax = 0.99
        self.thetagrid = np.linspace(self.thetamin,self.thetamax,self.ntheta,dtype=self.dtype)
        
        
        
        
        
        
        # construct finer grid for bargaining
        ntheta_fine = 10*self.ntheta # actual number may be a bit bigger
        self.thetagrid_fine = np.unique(np.concatenate( (self.thetagrid,np.linspace(self.thetamin,self.thetamax,ntheta_fine,dtype=self.dtype)) ))
        self.ntheta_fine = self.thetagrid_fine.size
        
        i_orig = list()
        
        for theta in self.thetagrid:
            assert theta in self.thetagrid_fine
            i_orig.append(np.where(self.thetagrid_fine==theta)[0])
            
        assert len(i_orig) == self.thetagrid.size
        # allows to recover original gird points on the fine grid        
        self.theta_orig_on_fine = np.array(i_orig).flatten()
        self.v_thetagrid_fine = VecOnGrid(self.thetagrid,self.thetagrid_fine)
        # precomputed object for interpolation

            
        
        
        


        self.exo_grids = {'Female, single':exogrid['zf_t'],
                          'Male, single':exogrid['zm_t'],
                          'Female and child':exogrid['zf_t'],
                          'Couple and child':exogrid['all_t'],
                          'Couple, no children':exogrid['all_t']}
        self.exo_mats = {'Female, single':exogrid['zf_t_mat'],
                          'Male, single':exogrid['zm_t_mat'],
                          'Female and child':exogrid['zf_t_mat_by_l_sk'],
                          'Couple and child':exogrid['all_t_mat_by_l_k'],
                          'Couple, no children':exogrid['all_t_mat_by_l_nk']} # sparse version?
        
        
        self.utility_shifters = {'Female, single':0.0,
                                 'Male, single':0.0,
                                 'Female and child':p['u_shift_mar'],
                                 'Couple and child':p['u_shift_mar'],
                                 'Couple, no children':p['u_shift_coh']}
        
        
        # this pre-computes transition matrices for meeting a partner
        zf_t_partmat = [self.mar_mats_iexo(t,female=True) if t < p['T'] - 1 else None 
                            for t in range(p['T'])]
        zm_t_partmat = [self.mar_mats_iexo(t,female=False) if t < p['T'] - 1 else None 
                            for t in range(p['T'])]
        
        self.part_mats = {'Female, single':zf_t_partmat,
                          'Male, single':  zm_t_partmat,
                          'Female and cild':  None,
                          'Couple and child': None,
                          'Couple, no children': None} # last is added for consistency
        
        self.mar_mats_assets()
        
        self.mar_mats_combine()
        
        
        # building m grid
        ezfmin = min([np.min(np.exp(g+t)) for g,t in zip(exogrid['zf_t'],p['f_wage_trend'])])
        ezmmin = min([np.min(np.exp(g+t)) for g,t in zip(exogrid['zm_t'],p['m_wage_trend'])])
        ezfmax = max([np.max(np.exp(g+t)) for g,t in zip(exogrid['zf_t'],p['f_wage_trend'])])
        ezmmax = max([np.max(np.exp(g+t)) for g,t in zip(exogrid['zm_t'],p['m_wage_trend'])])
        
        
        mmin = 0.95*ezmmin
        mmax = ezfmax + ezmmax + np.max(self.pars['R_t'])*self.amax
        self.mgrid_c = np.linspace(mmin,mmax,600)
        
        mmin_s = 0.95*min(self.ls_levels_sk)*ezfmin
        mmax_s = ezfmax + ezmmax + np.max(self.pars['R_t'])*self.amax
        self.mgrid_s = np.linspace(mmin_s,mmax_s,600)
        
        self.u_precompute()
        self.unplanned_pregnancy_probability()
        
        
        
    def mar_mats_assets(self,npoints=4,abar=0.1):
        # for each grid point on single's grid it returns npoints positions
        # on (potential) couple's grid's and assets of potential partner 
        # (that can be off grid) and correpsonding probabilities. 
        
        
        na = self.agrid_s.size
        
        agrid_s = self.agrid_s
        agrid_c = self.agrid_c
        
        s_a_partner = self.pars['sig_partner_a']
        
        
        prob_a_mat = np.zeros((na,npoints),dtype=self.dtype)
        i_a_mat = np.zeros((na,npoints),dtype=np.int16)
        
        
        
        for ia, a in enumerate(agrid_s):
            lagrid_t = np.zeros_like(agrid_c)
            
            i_neg = (agrid_c <= max(abar,a) - 1e-6)
            
            lagrid_t[~i_neg] = np.log(2e-6 + (agrid_c[~i_neg] - a)/max(abar,a))
            lmin = lagrid_t[~i_neg].min()
            # just fill with very negative values so this is never chosen
            lagrid_t[i_neg] = lmin - s_a_partner*10 - \
                s_a_partner*np.flip(np.arange(i_neg.sum())) 
            
            # TODO: this needs to be checked
            p_a = int_prob(lagrid_t,mu=0,sig=s_a_partner,n_points=npoints)
            i_pa = (-p_a).argsort()[:npoints] # this is more robust then nonzero
            p_pa = p_a[i_pa]
            prob_a_mat[ia,:] = p_pa
            i_a_mat[ia,:] = i_pa
        
        
        self.prob_a_mat = prob_a_mat
        self.i_a_mat = i_a_mat
            
            
        
    
    def mar_mats_iexo(self,t,female=True,trim_lvl=0.001):
        # TODO: check timing
        # this returns transition matrix for single agents into possible couples
        # rows are single's states
        # columnts are couple's states
        # you have to transpose it if you want to use it for integration
        setup = self
        
        nexo = setup.pars['nexo_t'][t]
        sigma_psi_init = setup.pars['sigma_psi_init']
        sig_z_partner = setup.pars['sig_partner_z']
        psi_couple = setup.exogrid.psi_t[t+1]
        
        
        if female:
            nz_single = setup.exogrid.zf_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zf_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zm_t[t+1]
            zmat_own = setup.exogrid.zf_t_mat[t]
        else:
            nz_single = setup.exogrid.zm_t[t].shape[0]
            p_mat = np.empty((nexo,nz_single))
            z_own = setup.exogrid.zm_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zf_t[t+1]
            zmat_own = setup.exogrid.zm_t_mat[t]    
            
        def ind_conv(a,b,c): return setup.all_indices(t,(a,b,c))[0]
        
        
        for iz in range(n_zown):
            p_psi = int_prob(psi_couple,mu=0,sig=sigma_psi_init)
            if female:
                p_zm  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
                p_zf  = zmat_own[iz,:]
            else:
                p_zf  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
                p_zm  = zmat_own[iz,:]
            #sm = sf
        
            p_vec = np.zeros(nexo)
            
            for izf, p_zf_i in enumerate(p_zf):
                if p_zf_i < trim_lvl: continue
            
                for izm, p_zm_i in enumerate(p_zm):
                    if p_zf_i*p_zm_i < trim_lvl: continue
                
                    for ipsi, p_psi_i in enumerate(p_psi):                    
                        p = p_zf_i*p_zm_i*p_psi_i
                        
                        if p > trim_lvl:
                            p_vec[ind_conv(izf,izm,ipsi)] = p    
                            
            assert np.any(p_vec>trim_lvl), 'Everything is zero?'              
            p_vec = p_vec / np.sum(p_vec)
            p_mat[:,iz] = p_vec
            
        return p_mat.T
    
    
    def mar_mats_combine(self):
        # for time moment and each position in grid for singles (ia,iz)
        # it computes probability distribution over potential matches
        # this is relevant for testing and simulations
        
        
        pmat_a = self.prob_a_mat
        imat_a = self.i_a_mat
        
        self.matches = dict()
        
        for female in [True,False]:
            desc = 'Female, single' if female else 'Male, single'
            
            pmats = self.part_mats[desc] 
            
            
            match_matrix = list()
            
            
            for t in range(self.pars['T']-1):
                pmat_iexo = pmats[t] # nz X nexo
                # note that here we do not use transpose
                
                nz = pmat_iexo.shape[0]
                
                inds = np.where( np.any(pmat_iexo>0,axis=0) )[0]
                inds_fem = self.all_indices(t,inds)[1]
                
                npos_iexo = inds.size
                npos_a = pmat_a.shape[1]
                npos = npos_iexo*npos_a
                pmatch = np.zeros((self.na,nz,npos),dtype=self.dtype)
                iamatch = np.zeros((self.na,nz,npos),dtype=np.int32)
                iexomatch = np.zeros((self.na,nz,npos),dtype=np.int32)
                ifemmatch = np.zeros((self.na,nz,npos),dtype=np.int32)
                
                i_conv = np.zeros((npos_iexo,npos_a),dtype=np.int32)
                
                for ia in range(npos_a):
                    i_conv[:,ia] = np.arange(npos_iexo*ia,npos_iexo*(ia+1))
                 
                
                for iz in range(nz):
                    probs = pmat_iexo[iz,inds]
                    
                    for ia in range(npos_a):
                        
                        pmatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            (pmat_a[:,ia][:,None])*(probs[None,:])
                        iamatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            imat_a[:,ia][:,None]
                        iexomatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            inds[None,:]
                        ifemmatch[:,iz,(npos_iexo*ia):(npos_iexo*(ia+1))] = \
                            inds_fem[None,:]
                            
                            
                        
                assert np.allclose(np.sum(pmatch,axis=2),1.0)
                match_matrix.append({'p':pmatch,'ia':iamatch,'iexo':iexomatch,'iconv':i_conv,
                                         'ifem':ifemmatch})
                    
            self.matches[desc] = match_matrix
         
        
    
    
    def all_indices(self,t,ind_or_inds=None):
        
        # just return ALL indices if no argument is called
        if ind_or_inds is None: 
            ind_or_inds = np.array(range(self.pars['nexo_t'][t]))
        
        if isinstance(ind_or_inds,tuple):
            izf,izm,ipsi = ind_or_inds
            ind = izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t] + izm*self.pars['n_psi_t'][t] + ipsi
        else:
            ind = ind_or_inds
            izf = ind // (self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t])
            izm = (ind - izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t]) // self.pars['n_psi_t'][t]
            ipsi = ind - izf*self.pars['n_zm_t'][t]*self.pars['n_psi_t'][t] - izm*self.pars['n_psi_t'][t]
            
        return ind, izf, izm, ipsi

    
    # functions u_mult and c_mult are meant to be shape-perservings
    
    def u_mult(self,theta):
        assert np.all(theta > 0) and np.all(theta < 1)
        powr = (1+self.pars['couple_rts'])/(self.pars['couple_rts']+self.pars['crra_power'])
        tf = theta
        tm = 1-theta
        ces = (tf**powr + tm**powr)**(1/powr)
        umult = (self.pars['A']**(1-self.pars['crra_power']))*ces
        
        
        
        assert umult.shape == theta.shape
        
        return umult
    
    
    def c_mult(self,theta):
        assert np.all(theta > 0) and np.all(theta < 1)
        powr = (1+self.pars['couple_rts'])/(self.pars['couple_rts']+self.pars['crra_power'])
        irho = 1/(1+self.pars['couple_rts'])
        irs  = 1/(self.pars['couple_rts']+self.pars['crra_power'])
        tf = theta
        tm = 1-theta
        bottom = (tf**(powr) + tm**(powr))**irho 
        
        kf = self.pars['A']*(tf**(irs))/bottom
        km = self.pars['A']*(tm**(irs))/bottom
        
        assert kf.shape == theta.shape
        assert km.shape == theta.shape
        
        return kf, km
    
    def u(self,c):
        return u_aux(c,self.pars['crra_power'])#(c**(1-self.pars['crra_power']))/(1-self.pars['crra_power'])
    
    
    def u_pub(self,x,l):
        alp = self.pars['util_alp']
        xi = self.pars['util_xi']
        lam = self.pars['util_lam']
        kap = self.pars['util_kap']        
        return alp*(x**lam + kap*(1-l)**lam)**((1-xi)/lam)/(1-xi)
    
    
    def u_part_k(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        kf, km = self.c_mult(theta)   
        l = self.ls_levels_k[il]
        upub = self.u_pub(x,l) + ushift + psi
        return self.u(kf*c) + upub, self.u(km*c) + upub
    
    def u_couple_k(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta) 
        l = self.ls_levels_k[il]
        return umult*self.u(c) + self.u_pub(x,l) + ushift + psi
    
    def u_part_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        kf, km = self.c_mult(theta)   
        upub = ushift + psi
        return self.u(kf*c) + upub, self.u(km*c) + upub
    
    def u_couple_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta)         
        return umult*self.u(c) + ushift + psi
    
    def u_single_k(self,c,x,il,ushift):
        umult = 1.0
        l = self.ls_levels_sk[il]
        return umult*self.u(c) + self.u_pub(x,l) + ushift
    
    def vm_last_grid(self,ushift,haschild):
        # this returns value of vm on the grid corresponding to vm
        s = self.agrid_c[:,None]
        zm = self.exogrid.all_t[-1][:,1][None,:]
        zf = self.exogrid.all_t[-1][:,0][None,:]
        psi = self.exogrid.all_t[-1][:,2][None,:,None]
        theta = self.thetagrid[None,None,:]
        
        nl = self.nls_k if haschild else self.nls_nk
        
        na, nexo, ntheta = self.na, self.pars['nexo_t'][-1], self.ntheta
        
        shp = (na,nexo,ntheta,nl)
        
        u_couple_g, u_f_g, u_m_g, income_g, c_g, x_g =np.zeros((6,) + shp,dtype=self.dtype)
        
        
        ftrend = self.pars['f_wage_trend'][-1]
        mtrend = self.pars['m_wage_trend'][-1]
        
        
        if haschild:
            lsl, upc, ucoup, upart = \
                self.ls_levels_k, self.ucouple_precomputed_x_k, self.u_couple_k, self.u_part_k
        else:
            lsl, upc, ucoup, upart = \
                self.ls_levels_nk, self.ucouple_precomputed_x_nk, self.u_couple_nk, self.u_part_nk
        
        
        for il in range(len(lsl)):
           
            inc = self.pars['R_t'][-1]*s + np.exp(zm+mtrend) +  np.exp(zf+ftrend)*lsl[il]
            income_g[...,il]  = inc[...,None]
            
            
            for itheta in range(ntheta):
                
                vals = upc[:,itheta,il]
                x_g[...,itheta,il] = np.interp(inc,self.mgrid_c,vals)
                c_g[...,itheta,il] = inc - x_g[...,itheta,il]
            
            u_couple_g[...,il] = ucoup(c_g[...,il],x_g[...,il],il,theta,ushift,psi)
            u_f_g[...,il], u_m_g[...,il] = upart(c_g[...,il],x_g[...,il],il,theta,ushift,psi)
             
        #Get optimal FLS
        ls=np.argmax(u_couple_g,axis=3)
        lsi=ls[...,None]
        u_c, u_f, u_m, x, c = (np.take_along_axis(x,lsi,axis=3).squeeze(axis=3)
                                for x in (u_couple_g,u_f_g,u_m_g,x_g,c_g))
        
        V  = u_c 
        VM = u_m 
        VF = u_f 
        
        return V.astype(self.dtype), VF.astype(self.dtype), VM.astype(self.dtype), c.astype(self.dtype), x.astype(self.dtype), np.zeros_like(c).astype(self.dtype), ls.astype(np.int16), u_couple_g.astype(self.dtype)
    
    

    def vs_last(self,s,z_plus_trend,ushift,return_cs=False):  
        # generic last period utility for single agent
        income = self.pars['R_t'][-1]*s+np.exp(z_plus_trend) 
        if return_cs:
            return self.u(income).astype(self.dtype) + ushift, income.astype(self.dtype), np.zeros_like(income.astype(self.dtype))
        else:
            return self.u(income)
    
    def vs_last_grid(self,female,ushift,return_cs=False):
        # this returns value of vs on the grid corresponding to vs
        s_in = self.agrid_s[:,None]
        z_in = self.exogrid.zf_t[-1][None,:] if female else self.exogrid.zm_t[-1][None,:]
        trend = self.pars['f_wage_trend'][-1] if female else self.pars['m_wage_trend'][-1]        
        return self.vs_last(s_in,z_in+trend,ushift,return_cs)
        
    
    
    def vsk_last_grid(self,ushift):
        # this returns value of vm on the grid corresponding to vm
        s = self.agrid_s[:,None]
        zf = self.exogrid.zf_t[-1][None,:]
        
        nl = self.nls_sk
        
        na, nexo = self.na, self.pars['n_zf_t'][-1]
        
        shp = (na,nexo,nl)
        
        u_g, income_g, c_g, x_g =np.zeros((4,) + shp,dtype=self.dtype)
        
        
        ftrend = self.pars['f_wage_trend'][-1]
        
        
   
        lsl, upc, ufun = \
            self.ls_levels_sk, self.ucouple_precomputed_x_sk, self.u_single_k
        
        
        for il in range(len(lsl)):
           
            inc = self.pars['R_t'][-1]*s + np.exp(zf+ftrend)*lsl[il]
            income_g[...,il]  = inc
            
            vals = upc[...,il]
            x_g[...,il] = np.interp(inc,self.mgrid_s,vals)
            c_g[...,il] = inc - x_g[...,il]
            
                
            #(c,x,il,ushift)
            u_g[...,il] = ufun(c_g[...,il],x_g[...,il],il,ushift)
        assert np.all(c_g>0)
        #Get optimal FLS
        ls=np.argmax(u_g,axis=-1)
        lsi=ls[...,None]
        u_c, x, c = (np.take_along_axis(x,lsi,axis=-1).squeeze(axis=-1)
                                for x in (u_g,x_g,c_g))
        
        V  = u_c 
        
        return V.astype(self.dtype), c.astype(self.dtype), x.astype(self.dtype), np.zeros_like(c).astype(self.dtype), ls.astype(np.int16), u_g.astype(self.dtype)
    
    
    def u_precompute(self):
        from intratemporal import int_with_x
        sig = self.pars['crra_power']
        alp = self.pars['util_alp']
        xi = self.pars['util_xi']
        lam = self.pars['util_lam']
        kap = self.pars['util_kap']
        qbar = self.pars['util_qbar']
        
        nm = self.mgrid_c.size
        ntheta = self.ntheta
        nl = self.nls_k
        
        uout = np.empty((nm,ntheta,nl),dtype=np.float32)
        xout = np.empty((nm,ntheta,nl),dtype=np.float32)
        
        for il in range(nl):
            for itheta in range(ntheta):
                A = self.u_mult(self.thetagrid[itheta])
                ls = self.ls_levels_k[il]
                x, c, u, q = int_with_x(self.mgrid_c,A=A,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,
                                        qlb=qbar,lbr=ls)
                uout[:,itheta,il] = u
                xout[:,itheta,il] = x
                
                
        self.ucouple_precomputed_u_k = uout
        self.ucouple_precomputed_x_k = xout
        
        
        
        
        from intratemporal import int_no_x
        
        nm = self.mgrid_c.size
        ntheta = self.ntheta
        nl = self.nls_nk
        
        uout = np.empty((nm,ntheta,nl),dtype=np.float32)
        xout = np.empty((nm,ntheta,nl),dtype=np.float32)
        
        for il in range(nl):
            for itheta in range(ntheta):
                A = self.u_mult(self.thetagrid[itheta])
                ls = self.ls_levels_nk[il]
                x, c, u, q = int_no_x(self.mgrid_c,A=A,sig=sig)
                uout[:,itheta,il] = u
                xout[:,itheta,il] = x
                
                
        self.ucouple_precomputed_u_nk = uout
        self.ucouple_precomputed_x_nk = xout
        
        
        nm = self.mgrid_s.size
        nl = self.nls_sk
        
        uout_s = np.empty((nm,nl),dtype=np.float32)
        xout_s = np.empty((nm,nl),dtype=np.float32)
        
        
        for il in range(nl):        
            A = 1.0
            ls = self.ls_levels_sk[il]
            x, c, u, q = int_with_x(self.mgrid_s,A=A,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,
                                    qlb=qbar,lbr=ls)
            uout_s[:,il] = u
            xout_s[:,il] = x

                
        self.ucouple_precomputed_u_sk = uout_s
        self.ucouple_precomputed_x_sk = xout_s
    
    def _unplanned_pregnancy_probability_fun(self,t,z):
        p = self.pars['preg_a0'] + self.pars['preg_at']*t + \
            self.pars['preg_az']*z + self.pars['preg_azt']*t*z
        return np.clip(p,0.0,1.0)
    
    def unplanned_pregnancy_probability(self):
        self.upp_precomputed = list()
        for t in range(self.pars['T']):
            zf = self.exogrid.zf_t[t]
            pf = self._unplanned_pregnancy_probability_fun(t,zf)
            self.upp_precomputed.append(pf)
            
        
                
    

#from numba import jit
#@jit(nopython=True)
def u_aux(c,sigma):
    # this is pretty important not to have (c^sigma - 1) here as it is hard to 
    # keep everywhere and occasionally this generates nasty mistakes
    if sigma!=1:
        return (c**(1-sigma))/(1-sigma)
    else:
        return np.log(c)

    


class DivorceCosts(object):
    # this is something that regulates divorce costs
    # it aims to be fully flexible
    def __init__(self, 
                 unilateral_divorce=True, # whether to allow for unilateral divorce
                 assets_kept = 1.0, # how many assets of couple are splited (the rest disappears)
                 u_lost_m=0.0,u_lost_f=0.0, # pure utility losses b/c of divorce
                 money_lost_m=0.0,money_lost_f=0.0, # pure money (asset) losses b/c of divorce
                 money_lost_m_ez=0.0,money_lost_f_ez=0.0, # money losses proportional to exp(z) b/c of divorce
                 eq_split=0.0 #The more of less equal way assets are split within divorce
                 ): # 
        
        self.unilateral_divorce = unilateral_divorce # w
        self.assets_kept = assets_kept
        self.u_lost_m = u_lost_m
        self.u_lost_f = u_lost_f
        self.money_lost_m = money_lost_m
        self.money_lost_f = money_lost_f
        self.money_lost_m_ez = money_lost_m_ez
        self.money_lost_f_ez = money_lost_f_ez
        self.eq_split = eq_split
        
    def shares_if_split(self,income_share_f):
        
        
        shf=(0.5*self.eq_split + income_share_f*(1-self.eq_split))
        share_f = self.assets_kept*shf - self.money_lost_f
        share_m = self.assets_kept*(1-shf) - self.money_lost_m
        
        return share_f, share_m
        
from rw_approximations import normcdf_tr

def tauchen_drift(z_now,z_next,rho,sigma,mu):
    z_now = np.atleast_1d(z_now)
    z_next = np.atleast_1d(z_next)
    if z_next.size == 1:
        return np.ones((z_now.size,1),dtype=np.float32)
    
    d = np.diff(z_next)
    assert np.ptp(d) < 1e-5, 'Step size should be fixed'
    
    h_half = d[0]/2
    
    Pi = np.zeros((z_now.size,z_next.size),dtype=np.float32)
    
    ez = rho*z_now + mu
    
    Pi[:,0] = normcdf_tr( ( z_next[0] + h_half - ez )/sigma)
    Pi[:,-1] = 1 - normcdf_tr( (z_next[-1] - h_half - ez ) / sigma )
    for j in range(1,z_next.size - 1):
        Pi[:,j] = normcdf_tr( ( z_next[j] + h_half - ez )/sigma) - \
                    normcdf_tr( ( z_next[j] - h_half - ez )/sigma)
    return Pi



def build_s_grid(agrid,n_between,da_min,da_max):
    sgrid = np.array([0.0],np.float64)
    for j in range(agrid.size-1):
        step = (agrid[j+1] - agrid[j])/n_between
        if step >= da_min and step <= da_max:
            s_add = np.linspace(agrid[j],agrid[j+1],n_between)[:-1]
        elif step < da_min:
            s_add = np.arange(agrid[j],agrid[j+1],da_min)
        elif step > da_max:
            s_add = np.arange(agrid[j],agrid[j+1],da_max)
        sgrid = np.concatenate((sgrid,s_add))
    
    sgrid = np.concatenate((sgrid,np.array([agrid[-1]])))
            
    if sgrid[0] == sgrid[1]: 
        sgrid = sgrid[1:]
        
    return sgrid
