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
        Tret = 45 # first period when the agent is retired
        Tfert = 18 # first peroid when infertile
        Tdiv = 44 # first period when cannot divorce / renegotiate
        Tmeet = 25 # first period when stop meeting partners
        p['T'] = T
        p['Tret'] = Tret
        p['Tfert'] = Tfert
        p['Tsim'] = T
        p['n_zf_t']      = [7]*Tret + [1]*(T-Tret)
        p['n_zm_t']      = [5]*Tret + [1]*(T-Tret)
        p['sigma_psi_mult'] = 0.28
        p['sigma_psi']   = 0.11
        p['R_t'] = [1.025]*T
        p['n_psi_t']     = [11]*T
        p['beta_t'] = [0.98]*T
        p['A'] = 1.0 # consumption in couple: c = (1/A)*[c_f^(1+rho) + c_m^(1+rho)]^(1/(1+rho))
        p['crra_power'] = 1.5
        p['couple_rts'] = 0.23    
        p['sig_partner_a'] = 0.1
        p['mu_partner_a_female'] = 0.00
        p['mu_partner_a_male'] = -0.00
        p['dump_factor_z'] = 0.85
        p['sig_partner_z'] = 1.2
        p['mu_partner_z_male'] = -0.02
        p['mu_partner_z_female'] = 0.02
        p['m_bargaining_weight'] = 0.5
        p['pmeet_0'] = 0.5
        p['pmeet_t'] = 0.0
        p['pmeet_t2'] = 0.0
        
        p['m_zf'] = 1.0
        p['m_zf0'] = 1.0
        
        p['z_drift'] = -0.09
        p['no kids at meeting'] = True
        p['high education'] = True # what trend to pick
        
        p['wret'] = 0.8
        p['uls'] = 0.2
        p['pls'] = 1.0
        
        
        p['preg_mult'] = 1.0
        
        p['u_shift_mar'] = 0.0
        p['u_shift_coh'] = 0.0
        p['sm_shift'] = 0.0
        
        p['disutil_marry_sm_fem_coef'] = 0.0
        p['disutil_marry_sm_mal_coef'] = 10.0
        p['disutil_shotgun_coef'] = 2.0
        p['pmeet_multiplier_fem'] = 1.0
        p['p_to_meet_sm_if_mal'] = 0.1
        
        p['taste_shock_mult'] = 1.0
        
       
        
        p['child_a_cost'] = 0.0
        
        
        '''
            Condition is h_male==1
            Age 24, h_male==1: .44558678
            Age 30, h_male==1: .45272711
            Implied standard deviation if h_male==1 is .03269622
            Implied initial if h_male==1 is .44197336
            Condition is h_female==1
            Age 24, h_female==1: .41251325
            Age 30, h_female==1: .42495642
            Implied standard deviation if h_female==1 is .04167489
            Implied initial if h_female==1 is .40614873
        '''
        
        if p['high education']:            
            p['sig_zm']    = 0.03269622
            p['sig_zm_0']  = 0.44197336
            
            
            p['sig_zf']    = p['m_zf']*0.04167489
            p['sig_zf_0']  = p['m_zf0']*0.40614873
            
        else:
            p['sig_zm']    = 0.09294824
            p['sig_zm_0']  = 0.40376647
            
            
            p['sig_zf']    = 0.08942736
            p['sig_zf_0']  = 0.39882978
            
        
        
        # college
            
            
        
        
        p['util_lam'] = 0.7
        p['util_alp'] = 0.5
        p['util_xi'] = 1.5
        p['util_kap'] = 0.5
        p['util_qbar'] = 0.0
        
        p['util_out_lf'] = 0.0
        p['preg_a0'] = 0.05
        p['preg_at'] = 0.01
        p['preg_at2'] = -0.001
        
        
        for key, value in kwargs.items():
            assert (key in p), 'wrong name?'
            p[key] = value
        
        '''
              h_female_T |   .0685814   .0004133   165.92   0.000     .0677713    .0693915
             h_female_T2 |  -.0038631   .0000449   -86.03   0.000    -.0039511   -.0037751
             h_female_T3 |   .0000715   1.30e-06    55.20   0.000      .000069    .0000741
                  h_male |   .0818515   .0018308    44.71   0.000     .0782632    .0854398
                h_male_T |   .0664369   .0005108   130.07   0.000     .0654358    .0674379
               h_male_T2 |  -.0032096   .0000528   -60.77   0.000    -.0033131   -.0031061
               h_male_T3 |   .0000529   1.48e-06    35.68   0.000       .00005    .0000558
                l_female |  -.3715448   .0015389  -241.43   0.000    -.3745611   -.3685285
              l_female_T |   .0258538   .0003048    84.83   0.000     .0252564    .0264512
             l_female_T2 |  -.0012931   .0000378   -34.20   0.000    -.0013672    -.001219
             l_female_T3 |   .0000273   1.15e-06    23.79   0.000     .0000251    .0000296
                  l_male |  -.2587716   .0014324  -180.66   0.000     -.261579   -.2559641
                l_male_T |   .0348409   .0002597   134.14   0.000     .0343319      .03535
               l_male_T2 |  -.0013866   .0000328   -42.27   0.000    -.0014509   -.0013223
               l_male_T3 |    .000024   1.01e-06    23.79   0.000     .0000221     .000026
        '''
        
        if p['high education']:
            p['f_wage_trend'] = np.array(
                                [0.0 + 
                                 0.0685814*(min(t,30) - 5)
                                 -.0038631*((min(t,30)-5)**2)
                                 + 0.0000715*((min(t,30)-5)**3)
                                             for t in range(T)]
                                        )
            
            p['m_wage_trend'] = np.array(
                                        [0.0818515 + 0.0664369*(min(t,30)-5)  
                                             -.0032096*((min(t,30)-5)**2) 
                                             + 0.0000529*((min(t,30)-5)**3)
                                             for t in range(T)]
                                        )
        
        else:
        # no college
        
            p['f_wage_trend'] = np.array(
                                [-0.3668214 + 
                                 0.0264887*(min(t,30) - 5)
                                 -0.0012464*((min(t,30)-5)**2)
                                 +0.0000251*((min(t,30)-5)**3)
                                             for t in range(T)]
                                        )
            
            p['m_wage_trend'] = np.array(
                                        [-0.2424105 + 0.037659*(min(t,30)-5)  
                                             -0.0015337*((min(t,30)-5)**2) 
                                             + 0.000026*((min(t,30)-5)**3)
                                             for t in range(T)]
                                        )
        
        
        
        # derivative parameters
        p['sigma_psi_init'] = p['sigma_psi_mult']*p['sigma_psi']
        
        
        p['disutil_marry_sm_mal'] = p['disutil_marry_sm_mal_coef']*p['u_shift_mar']
        p['disutil_marry_sm_fem'] = p['disutil_marry_sm_fem_coef']*p['u_shift_mar']
        p['disutil_shotgun'] =  p['disutil_shotgun_coef']*p['sigma_psi_init']
        
        p['preg_az'] =  0.00
        p['preg_azt'] = 0.00
        
        #Get the probability of meeting, adjusting for year-period
           
        
        p['taste_shock'] = p['taste_shock_mult']*p['sigma_psi']
        
        p['is fertile'] = [True]*Tfert + [False]*(T-Tfert)
        p['can divorce'] = [True]*Tdiv + [False]*(T-Tdiv)        
        p['pmeet_t'] = [np.clip(p['pmeet_0'] + (t-9)*p['pmeet_t'] + ((t-9)**2)*p['pmeet_t2'],0.0,1.0) for t in range(Tmeet)] + [0.0]*(T-Tmeet)
        #p['poutsm_t'] = [p['poutsm']]*T
        
        
        self.pars = p
        
        self.dtype = np.float64 # type for all floats
        
       
        
        # relevant for integration
        self.state_names = ['Female, single','Male, single','Female and child','Couple, no children','Couple and child']
        
        # female labor supply
        
        lmin = 0.2
        lmax = 1.0
        nl = 2
        
        ls = np.array([0.2,1.0]) #np.linspace(lmin,lmax,nl,dtype=self.dtype)
        ps = np.array([p['pls'],0.0])
        ls_ushift = np.array([p['util_out_lf'],0.0])
        
        
        self.ls_levels = dict()
        self.ls_levels['Couple, no children'] = np.array([1.0],dtype=self.dtype)
        self.ls_levels['Female, single'] = np.array([1.0],dtype=self.dtype)
        self.ls_levels['Male, single'] = np.array([1.0],dtype=self.dtype)
        self.ls_levels['Couple and child'] = ls
        self.ls_levels['Female and child'] = ls
        
        self.ls_ushift = dict()
        self.ls_ushift['Couple, no children'] = np.array([0.0],dtype=self.dtype)
        self.ls_ushift['Female, single'] = np.array([0.0],dtype=self.dtype)
        self.ls_ushift['Male, single'] = np.array([0.0],dtype=self.dtype)
        self.ls_ushift['Couple and child'] = ls_ushift
        self.ls_ushift['Female and child'] = ls_ushift
        
        
        
        
        
        
        #self.ls_utilities = np.array([p['uls'],0.0],dtype=self.dtype)
        self.ls_pdown = dict()
        self.ls_pdown['Couple, no children'] = np.array([0.0],dtype=self.dtype)
        self.ls_pdown['Female, single'] = np.array([0.0],dtype=self.dtype)
        self.ls_pdown['Male, single']   = np.array([0.0],dtype=self.dtype)
        self.ls_pdown['Female and child'] = ps
        self.ls_pdown['Couple and child'] = ps
        self.nls = dict()
        self.nls['Couple and child'] = len(self.ls_levels['Couple and child'])
        self.nls['Couple, no children'] = len(self.ls_levels['Couple, no children'])
        self.nls['Female and child'] = len(self.ls_levels['Female and child'])
        self.nls['Female, single'] = len(self.ls_levels['Female, single'])
        self.nls['Male, single'] = len(self.ls_levels['Male, single'])
        
        
        
        
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
                               for p in self.ls_pdown['Couple, no children'] ]
            
            all_t_mat_by_l_spt_nk = [ [(1-p)*m + p*md if m is not None else None
                                    for m, md in zip(all_t_mat_sparse_T,all_t_mat_down_sparse_T)]
                               for p in self.ls_pdown['Couple, no children'] ]
            
            all_t_mat_by_l_k = [ [(1-p)*m + p*md if m is not None else None 
                                for m , md in zip(all_t_mat,all_t_mat_down)]
                               for p in self.ls_pdown['Couple and child'] ]
            
            all_t_mat_by_l_spt_k = [ [(1-p)*m + p*md if m is not None else None
                                    for m, md in zip(all_t_mat_sparse_T,all_t_mat_down_sparse_T)]
                               for p in self.ls_pdown['Couple and child'] ]
            
            zf_t_mat_by_l_sk = [ [(1-p)*m + p*md if md is not None else None 
                                for m , md in zip(exogrid['zf_t_mat'],zf_bad)]
                                   for p in self.ls_pdown['Female and child'] ]
            
            
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
        self.amax = 50
        self.agrid_c = np.linspace(self.amin**0.5,self.amax**0.5,self.na,dtype=self.dtype)**2
        #tune=1.5
        #self.agrid_c = np.geomspace(self.amin+tune,self.amax+tune,num=self.na)-tune
        
        # this builds finer grid for potential savings
        s_between = 7 # default numer of points between poitns on agrid
        s_da_min = 0.01 # minimal step (does not create more points)
        s_da_max = 0.1 # maximal step (creates more if not enough)
        
        self.sgrid_c = build_s_grid(self.agrid_c,s_between,s_da_min,s_da_max)
        self.vsgrid_c = VecOnGrid(self.agrid_c,self.sgrid_c)
        
        
         
        #Grid Single
        self.amin_s = 0
        self.amax_s = self.amax/2.0
        self.agrid_s = self.agrid_c/2.0
        #tune_s=1.5
        #self.agrid_s = np.geomspace(self.amin_s+tune_s,self.amax_s+tune_s,num=self.na)-tune_s
        
        self.sgrid_s = build_s_grid(self.agrid_s,s_between,s_da_min,s_da_max)
        self.vsgrid_s = VecOnGrid(self.agrid_s,self.sgrid_s)
        
        # grid for theta
        self.ntheta = 11
        self.thetamin = 0.01
        self.thetamax = 0.99
        self.thetagrid = np.linspace(self.thetamin,self.thetamax,self.ntheta,dtype=self.dtype)
        
        
        self.child_a_cost_single = np.minimum(self.agrid_s,self.pars['child_a_cost'])
        self.child_a_cost_couple = np.minimum(self.agrid_c,self.pars['child_a_cost'])
        
        self.vagrid_child_single = VecOnGrid(self.agrid_s, self.agrid_s - self.child_a_cost_single)
        self.vagrid_child_couple = VecOnGrid(self.agrid_c, self.agrid_c - self.child_a_cost_couple)
        
        
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
                                 'Female and child':p['u_shift_mar'] + p['sm_shift'],
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
        
        self.money_min = 0.95*min(self.ls_levels['Female and child'])*min(ezmmin,ezfmin) # cause FLS can be up to 0
        mmin = self.money_min
        mmax = ezfmax + ezmmax + np.max(self.pars['R_t'])*self.amax
        mint = (ezfmax + ezmmax) # poin where more dense grid begins
        
        ndense = 600
        nm = 1500
        
        gsparse = np.linspace(mint,mmax,nm-ndense)
        gdense = np.linspace(mmin,mint,ndense+1) # +1 as there is a common pt
        
        self.mgrid = np.zeros(nm,dtype=self.dtype)
        self.mgrid[ndense:] = gsparse
        self.mgrid[:(ndense+1)] = gdense
        self.mgrid_c = self.mgrid
        self.mgrid_s = self.mgrid
        assert np.all(np.diff(self.mgrid)>0)
        
        self.u_precompute()
        self.unplanned_pregnancy_probability()
        
        
        
    def _mar_mats_assets(self,npoints=4,female=True,upp=False,abar=0.1):
        # for each grid point on single's grid it returns npoints positions
        # on (potential) couple's grid's and assets of potential partner 
        # (that can be off grid) and correpsonding probabilities. 
        
        
        na = self.agrid_s.size
        
        agrid_s = self.agrid_s
        agrid_c = self.agrid_c
        
        s_a_partner = self.pars['sig_partner_a']
        mu_a_partner = self.pars['mu_partner_a_female'] if female else self.pars['mu_partner_a_male']
        
        
        prob_a_mat = np.zeros((na,npoints),dtype=self.dtype)
        i_a_mat = np.zeros((na,npoints),dtype=np.int16)
        
        if upp: assert female
        
        aloss = 0.0 if not upp else self.pars['child_a_cost']
        
        for ia, a in enumerate(agrid_s):
            lagrid_t = np.zeros_like(agrid_c)
            
            i_neg = (agrid_c <= max(abar,a) - 1e-6)
            
            lagrid_t[~i_neg] = np.log(2e-6 + (agrid_c[~i_neg] - a)/max(abar,a))
            lmin = lagrid_t[~i_neg].min()
            # just fill with very negative values so this is never chosen
            lagrid_t[i_neg] = lmin - s_a_partner*10 - \
                s_a_partner*np.flip(np.arange(i_neg.sum())) 
            
            # TODO: this needs to be checked
            p_a = int_prob(lagrid_t,mu=mu_a_partner-aloss,sig=s_a_partner,n_points=npoints)
            i_pa = (-p_a).argsort()[:npoints] # this is more robust then nonzero
            p_pa = p_a[i_pa]
            prob_a_mat[ia,:] = p_pa
            i_a_mat[ia,:] = i_pa
        
        return prob_a_mat, i_a_mat
            
    def mar_mats_assets(self,npoints=4,abar=0.1):
        self.prob_a_mat_female_noupp, self.i_a_mat_female_noupp = self._mar_mats_assets(npoints=npoints,female=True,upp=False,abar=abar)
        self.prob_a_mat_male, self.i_a_mat_male = self._mar_mats_assets(npoints=npoints,female=False,abar=abar)
        self.prob_a_mat_female_upp, self.i_a_mat_female_upp = self._mar_mats_assets(npoints=npoints,female=True,upp=True,abar=abar)
        
    
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
        mu_z_partner = setup.pars['mu_partner_z_female']  if female else setup.pars['mu_partner_z_male']
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
        
        
        df = setup.pars['dump_factor_z']
        
        for iz in range(n_zown):
            p_psi = int_prob(psi_couple,mu=0,sig=sigma_psi_init)
            if female:
                p_zm  = int_prob(z_partner, mu=df*z_own[iz] + mu_z_partner,sig=sig_z_partner)
                p_zf  = zmat_own[iz,:]
            else:
                p_zf  = int_prob(z_partner, mu=df*z_own[iz] + mu_z_partner,sig=sig_z_partner)
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
        
        
        self.matches = dict()
        
        for female, upp in [(True,False),(True,True),(False,False)]:
            desc = 'Female, single, upp' if (female and upp) else 'Female, single, no upp' if (female and not upp) else 'Male, single'
            
            
            if female and upp:
                pmat_a, imat_a = self.prob_a_mat_female_upp, self.i_a_mat_female_upp
            elif female and not upp:
                pmat_a, imat_a = self.prob_a_mat_female_noupp, self.i_a_mat_female_noupp
            else:
                pmat_a, imat_a = self.prob_a_mat_male, self.i_a_mat_male
            
            
            pmats = self.part_mats['Female, single'] if female else self.part_mats['Male, single']
            
            
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
        l = self.ls_levels['Couple and child'][il]
        us = self.ls_ushift['Couple and child'][il]
        upub = self.u_pub(x,l) + ushift + psi
        return self.u(kf*c) + upub + us, self.u(km*c) + upub
    
    def u_couple_k(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta) 
        l = self.ls_levels['Couple and child'][il]
        us = theta*self.ls_ushift['Couple and child'][il] 
        return umult*self.u(c) + self.u_pub(x,l) + us + ushift + psi
    
    def u_part_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        kf, km = self.c_mult(theta)   
        upub = ushift + psi
        return self.u(kf*c) + upub, self.u(km*c) + upub
    
    def u_couple_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta)         
        return umult*self.u(c) + ushift + psi
    
    def u_single_k(self,c,x,il,ushift):
        umult = 1.0
        l = self.ls_levels['Female and child'][il]
        us = self.ls_ushift['Female and child'][il]
        return umult*self.u(c) + self.u_pub(x,l) + us + ushift
    
    def u_precompute(self):
        
        
        self.u_precomputed = dict()
        self.u_precomputed['Couple and child'] = dict()
        self.u_precomputed['Couple, no children'] = dict()
        self.u_precomputed['Female, single'] = dict()
        self.u_precomputed['Male, single'] = dict()
        self.u_precomputed['Female and child'] = dict()
        
        from intratemporal import int_with_x
        sig = self.pars['crra_power']
        alp = self.pars['util_alp']
        xi = self.pars['util_xi']
        lam = self.pars['util_lam']
        kap = self.pars['util_kap']
        qbar = self.pars['util_qbar']
        
        
        
        # couple and child
        nm = self.mgrid_c.size
        ntheta = self.ntheta
        nl = self.nls['Couple and child']
        
        uout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        xout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        
        for il in range(nl):
            for itheta in range(ntheta):
                A = self.u_mult(self.thetagrid[itheta])
                ls = self.ls_levels['Couple and child'][il]
                x, c, u, q = int_with_x(self.mgrid_c,A=A,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,
                                        qlb=qbar,lbr=ls)
                tht = self.thetagrid[itheta]
                
                uout[:,itheta,il] = u + tht*self.ls_ushift['Couple and child'][il]
                xout[:,itheta,il] = x
                
                
        self.u_precomputed['Couple and child']['u'] = uout
        self.u_precomputed['Couple and child']['x'] = xout
        
        
        
        
        
        nm = self.mgrid_s.size
        nl = self.nls['Female and child']
        
        uout_s = np.empty((nm,1,nl),dtype=self.dtype)
        xout_s = np.empty((nm,1,nl),dtype=self.dtype)
        
        for il in range(nl):        
            A = 1.0
            ls = self.ls_levels['Female and child'][il]
            x, c, u, q = int_with_x(self.mgrid_s,A=A,alp=alp,sig=sig,xi=xi,lam=lam,kap=kap,
                                    qlb=qbar,lbr=ls)
            uout_s[:,0,il] = u + self.ls_ushift['Female and child'][il]
            xout_s[:,0,il] = x
               
        self.u_precomputed['Female and child']['u'] = uout_s
        self.u_precomputed['Female and child']['x'] = xout_s
        
        
        
        from intratemporal import int_no_x
        
        nm = self.mgrid_c.size
        ntheta = self.ntheta
        nl = self.nls['Couple, no children']
        
        uout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        xout = np.empty((nm,ntheta,nl),dtype=self.dtype)
        
        for il in range(nl):
            for itheta in range(ntheta):
                A = self.u_mult(self.thetagrid[itheta])
                ls = self.ls_levels['Couple, no children'][il]
                x, c, u, q = int_no_x(self.mgrid_c,A=A,sig=sig)
                uout[:,itheta,il] = u
                xout[:,itheta,il] = x
                
                
        self.u_precomputed['Couple, no children']['u'] = uout
        self.u_precomputed['Couple, no children']['x'] = xout
        
        
        
        from intratemporal import int_no_x
        
        
        nm = self.mgrid_s.size
        ntheta = self.ntheta
        nl = self.nls['Female, single'] # !!!
        
        uout = np.empty((nm,1,nl),dtype=self.dtype)
        xout = np.empty((nm,1,nl),dtype=self.dtype)
        
        for il in range(nl):
            A = 1
            ls = self.ls_levels['Female, single'][il]
            x, c, u, q = int_no_x(self.mgrid_s,A=A,sig=sig)
            uout[:,0,il] = u
            xout[:,0,il] = x
                
            
                
        self.u_precomputed['Female, single']['u'] = uout
        self.u_precomputed['Female, single']['x'] = xout
        
        self.u_precomputed['Male, single']['u'] = uout
        self.u_precomputed['Male, single']['x'] = xout
        
        
        
        
        
    
    def _unplanned_pregnancy_probability_fun(self,t,z):
        p = self.pars['preg_a0'] + self.pars['preg_at']*(t-9) + \
            self.pars['preg_at2']*((t-9)**2) + \
            self.pars['preg_az']*z + self.pars['preg_azt']*t*z
        #p = self.pars['preg_a0'] + self.pars['preg_at']*t + \
        #    self.pars['preg_az']*z + self.pars['preg_azt']*t*z
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
                 assets_kept = 0.9, # how many assets of couple are splited (the rest disappears)
                 u_lost_m=0.0,u_lost_f=0.0, # pure utility losses b/c of divorce
                 money_lost_m=0.0,money_lost_f=0.0, # pure money (asset) losses b/c of divorce
                 money_lost_m_ez=0.0,money_lost_f_ez=0.0, # money losses proportional to exp(z) b/c of divorce
                 eq_split=1.0 #The more of less equal way assets are split within divorce
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
        return np.ones((z_now.size,1),dtype=z_now.dtype)
    
    d = np.diff(z_next)
    assert np.ptp(d) < 1e-5, 'Step size should be fixed'
    
    h_half = d[0]/2
    
    Pi = np.zeros((z_now.size,z_next.size),dtype=z_now.dtype)
    
    ez = rho*z_now + mu
    
    Pi[:,0] = normcdf_tr( ( z_next[0] + h_half - ez )/sigma)
    Pi[:,-1] = 1 - normcdf_tr( (z_next[-1] - h_half - ez ) / sigma )
    for j in range(1,z_next.size - 1):
        Pi[:,j] = normcdf_tr( ( z_next[j] + h_half - ez )/sigma) - \
                    normcdf_tr( ( z_next[j] - h_half - ez )/sigma)
    return Pi



def build_s_grid(agrid,n_between,da_min,da_max):
    sgrid = np.array([0.0],agrid.dtype)
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
