#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains things relevant for setting up the model
"""

import numpy as np

from rw_approximations import tauchen_nonst, rouw_nonst
from mc_tools import combine_matrices_two_lists, int_prob,cut_matrix
from scipy.stats import norm
from collections import namedtuple
from gridvec import VecOnGrid
from aux_routines import prob_polyfit
from itertools import product
from tiktak import filer
from nopar_assets_dist import get_estimates

from scipy import sparse

try:
    import cupy as cp
    q = cp.array([1,2,3])
    assert cp.allclose(q.mean(),2.0)
    print('using cupy!')
    cupy = True
except:
    print('using cpu!')
    cupy = False



class ModelSetup(object):
    def __init__(self,nogrid=False,divorce_costs_k='Default',divorce_costs_nk='Default',**kwargs): 
        
        p = dict()       
        T = 55
        Tret = 45 # first period when the agent is retired
        Tfert = 18 # first peroid when infertile
        Tdiv = 44 # first period when cannot divorce / renegotiate
        Tmeet = 43 # first period when stop meeting partners
        Tinc = 25 # first period where stop tracking income process and assume it to be fixed
        p['T'] = T
        p['Tret'] = Tret
        p['Tfert'] = Tfert
        p['Tsim'] = T
        p['n_zf_t']      = [7]*Tret + [1]*(T-Tret)
        p['n_zm_t']      = [5]*Tret + [1]*(T-Tret)
        p['sigma_psi_init'] = 0.28
        p['sigma_psi']   = 0.11
        p['R_t'] = [1.025]*T
        p['n_psi'] = 15
        p['beta_t'] = [0.98]*T
        p['A'] = 1.0 # consumption in couple: c = (1/A)*[c_f^(1+rho) + c_m^(1+rho)]^(1/(1+rho))
        p['crra_power'] = 1.5
        p['couple_rts'] = 0.23    
        p['sig_partner_a'] = 0.1
        p['mu_partner_a_female'] = 0.00
        p['mu_partner_a_male'] = -0.00
        p['dump_factor_z'] = 0.75
        p['sig_partner_z'] = 0.25
        p['mu_partner_z_male'] = -0.02
        p['mu_partner_z_female'] = 0.02
        p['m_bargaining_weight'] = 0.5
        
        
        p['pmeet_21'] = 0.1
        p['pmeet_30'] = 0.2
        p['pmeet_40'] = 0.1
        
        p['m_zf'] = 0.9
        p['m_zf0'] = 1.0
        
        p['z_drift'] = -0.09
        p['no kids at meeting'] = True
        p['high education'] = True # what trend to pick
        p['any kids'] = True
        
        p['wret'] = 0.8
        p['uls'] = 0.2
        p['pls'] = 1.0
        
        p['income_sd_mult'] = 1.0
        p['pay_gap'] = True
        
        
        p['preg_mult'] = 1.0
        
        p['u_shift_mar'] = 0.0
        p['u_shift_coh'] = 0.0
        p['sm_shift'] = 0.0
        
        p['disutil_marry_sm_fem'] = 0.0
        p['disutil_marry_sm_mal'] = 10.0
        p['disutil_shotgun'] = 2.0
        p['pmeet_multiplier_fem'] = 1.0
        p['p_to_meet_sm_if_mal'] = 0.1
        
        p['taste_shock_mult'] = 1.0
        
        p['p_abortion_access'] = 0.5
        p['abortion_costs'] = 10.0
        
        p['u_lost_divorce'] = 0.0 
       
        
        p['child_a_cost'] = 0.0
        
        p['child_support_share'] = 0.2  
        p['child_support_awarded_nm'] = 0.284
        p['child_support_awarded_div'] = 0.461
            
        
        
        p['util_lam'] = 0.7
        p['util_alp'] = 0.5
        p['util_xi'] = 1.5
        p['util_kap'] = 0.5
        p['util_qbar'] = 0.0
        
        p['util_out_lf'] = 0.0
        
        
        
        p['preg_21'] = 0.01
        p['preg_28'] = 0.5
        p['preg_35'] = 0.3
        
        
        for key, value in kwargs.items():
            assert (key in p), 'wrong name?'
            p[key] = value
        
        
        
        
        if p['high education']:            
            p['sig_zm']    = p['income_sd_mult']*0.16138593
            p['sig_zm_0']  = p['income_sd_mult']*0.41966813 
            
            p['sig_zf']    = p['income_sd_mult']*p['m_zf']*0.19571624
            p['sig_zf_0']  = p['income_sd_mult']*p['m_zf0']*0.43351219
        else:
            p['sig_zm']    = p['income_sd_mult']*0.17195085
            p['sig_zm_0']  = p['income_sd_mult']*0.2268650
            
            
            p['sig_zf']    = p['income_sd_mult']*p['m_zf']*0.1762148
            p['sig_zf_0']  = p['income_sd_mult']*p['m_zf0']*0.1762148
            
        
        p['sm_init'] = (0.02 if p['high education'] else 0.25) if p['any kids'] else 0.0 # initial share of single moms
        
        
        # college
        
        if p['high education']: 
            
            
            m_trend_data = [0.0,0.11316599,.2496034,.31260625,.37472204,.4268548,.48067884,.52687573,
                            .57293878,.60941412,.65015743,.6685226,.72482815,.74446455,.76712521,.78038137,.79952806,
                            .80092523,.81972567,.82913486,.83849471,.84308452,.84646086,.85437072,.85499576]
            
            
            f_trend_data = [0.0,0.06715984,.21149606,.32283002,.46885336,.52712037,.58302632,.63348555,.68024646,
                             .71450132,.74246337,.77044807,.79946406,.80640353,.83799304,.85356081,.86832235
                             ,.87407447,.87820755,.86840901,.87630054,.8765972,.87894493,.87800553,.87932908]
            
            
            nm = len(m_trend_data)-1
            nf = len(f_trend_data)-1
            
            t0 = 4
            gap = 3.0340077 - 2.8180354 # male - female
            
            c_female = -f_trend_data[t0]
            c_male = gap - m_trend_data[t0]
            
            p['m_wage_trend'] = np.array(
                                        [c_male + m_trend_data[min(t,nm)]
                                             for t in range(T)]
                                        )
            p['f_wage_trend'] = np.array(
                                [c_female + f_trend_data[min(t,nf)]
                                             for t in range(T)]
                                        )
            
        else:
        # no college
            m_trend_data = [0.0,0.03824801,.13410532,.1663402,.18535393,.20804419,.22880115,.23963687,.26544877,.27445022,.28828461,.29908889,.3242355,.34399191,.35703786,.36160155,.37354454,.37365049,.38967079,.39410233,.40492857,.40538787,.42001778,.43326506,.43527713]
            
            
            f_trend_data = [0.0,0.03709545,.07178513,.09427489,.18766845,.20733048,.21432513,.22962527,.24421213,.25502674,.26330492,.26669114,.271962,.2775313,.29847667,.29413686,.30664712,.30294726,.31538057,.31768117,.32177537,.32804634,.32827188,.33866797,.34713842]
            
            
            nm = len(m_trend_data)-1
            nf = len(f_trend_data)-1
            
            t0 = 4
            gap = 2.5449782 - 2.3597982 # male - female
            
            c_female = -f_trend_data[t0]
            c_male = gap - m_trend_data[t0]
            
            p['m_wage_trend'] = np.array(
                                        [c_male + m_trend_data[min(t,nm)]
                                             for t in range(T)]
                                        )
            p['f_wage_trend'] = np.array(
                                [c_female + f_trend_data[min(t,nf)]
                                             for t in range(T)]
                                        )
            
        
        
        if not p['pay_gap']:
            p['sig_zf'], p['sig_zf_0'] = p['sig_zm'], p['sig_zm_0']
            p['f_wage_trend'] = p['m_wage_trend']
        
        
        p['preg_az'] =  0.00
        p['preg_azt'] = 0.00
        
        #Get the probability of meeting, adjusting for year-period
           
        
        p['taste_shock'] = 0.0 #*p['taste_shock_mult']*0.0#p['sigma_psi']
        
        p['is fertile'] = [p['any kids']]*Tfert + [False]*(T-Tfert)
        p['can divorce'] = [True]*Tdiv + [False]*(T-Tdiv)        
        #p['poutsm_t'] = [p['poutsm']]*T
        
        
        
        p['pmeet_0'],  p['pmeet_t'], p['pmeet_t2'] = prob_polyfit(
                    (p['pmeet_21'],0),(p['pmeet_30'],7),(p['pmeet_40'],14),
                                                                   max_power=2)
        
        p['preg_a0'],  p['preg_at'], p['preg_at2'] = prob_polyfit(
                    (p['preg_21'],0),(p['preg_28'],7),(p['preg_35'],14),
                                                                   max_power=2)
        
        
        p['pmeet_t'] = [np.clip(p['pmeet_0'] + t*p['pmeet_t'] + (t**2)*p['pmeet_t2'],0.0,1.0) for t in range(Tmeet)] + [0.0]*(T-Tmeet)
        
        
        
        p['n_psi_t'] = [p['n_psi']]*T
        
        p['psi_clip'] = 2.5*p['sigma_psi_init']
        
        
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
            self.divorce_costs_k = DivorceCosts(u_lost_m=self.pars['u_lost_divorce'],
                                                u_lost_f=self.pars['u_lost_divorce'])
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
            self.divorce_costs_nk = DivorceCosts(u_lost_m=self.pars['u_lost_divorce'],
                                              u_lost_f=self.pars['u_lost_divorce'])
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
            
            
            for t in range(Tinc,Tret):
                for key in ['zf_t','zf_t_mat','zm_t','zm_t_mat']:
                    exogrid[key][t] = exogrid[key][Tinc]
            
            
            for t in range(Tret,T):
                exogrid['zf_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zm_t'][t] = np.array([np.log(p['wret'])])
                exogrid['zf_t_mat'][t] = np.atleast_2d(1.0)
                exogrid['zm_t_mat'][t] = np.atleast_2d(1.0)
                
                
            # fix transition from non-retired to retired    
            exogrid['zf_t_mat'][Tret-1] = np.ones((p['n_zf_t'][Tret-1],1))
            exogrid['zm_t_mat'][Tret-1] = np.ones((p['n_zm_t'][Tret-1],1))
            
            #exogrid['psi_t'], exogrid['psi_t_mat'] = rouw_nonst(p['T'],p['sigma_psi'],p['sigma_psi_init'],p['n_psi_t'][0])
            exogrid['psi_t'], exogrid['psi_t_mat'] = tauchen_nonst(p['T'],p['sigma_psi'],p['sigma_psi_init'],p['n_psi_t'][0],nsd=2.5,fix_0=False)
            
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
            
            self.compute_child_support_transitions(child_support_share=p['child_support_share'])
            #assert False
            
            
            
        #Grid Couple
        self.na = 40
        self.amin = 0
        self.amax = 100
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
        
        assert self.pars['child_a_cost']<1e-3, 'not implemented'
        #self.vagrid_child_single = VecOnGrid(self.agrid_s, self.agrid_s - self.child_a_cost_single)
        #self.vagrid_child_couple = VecOnGrid(self.agrid_c, self.agrid_c - self.child_a_cost_couple)
        
        
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
       
        
        name_fem_pkl = 'az_dist_fem.pkl' if p['high education'] else 'az_dist_fem_noc.pkl'
        name_mal_pkl = 'az_dist_mal.pkl' if p['high education'] else 'az_dist_mal_noc.pkl'
        name_fem_csv = 'income_assets_distribution_male.csv' if p['high education'] else 'ia_male_noc.csv'
        name_mal_csv = 'income_assets_distribution_female.pkl' if p['high education'] else 'ia_female_noc.csv'
        # this is not an error, things are switched
        
        
        try:
            self.partners_distribution_fem = filer(name_fem_pkl,0,0,repeat=False)
            self.partners_distribution_mal = filer(name_mal_pkl,0,0,repeat=False)
        except:
            print('recreating estimates...')
            
            est_fem = get_estimates(fname=name_fem_csv,
                                    age_start=23,age_stop=42,
                                    zlist=self.exogrid.zm_t[2:])
            filer(name_fem_pkl,est_fem,True,repeat=False)
            self.partners_distribution_fem = est_fem
            est_mal = get_estimates(fname=name_mal_csv,
                                    age_start=21,age_stop=40,
                                    zlist=self.exogrid.zf_t[0:])
            filer(name_mal_pkl,est_mal,True,repeat=False)
            self.partners_distribution_mal = est_mal
            
        self.build_matches()
        
        
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
        self.compute_taxes()
        
        self.cupyfy()
        
        
    
    def _mar_mats(self,t,female=True):
        # this returns transition matrix for single agents into possible couples
        # rows are single's states
        # columnts are couple's states
        # you have to transpose it if you want to use it for integration
        setup = self
        
        nexo = setup.pars['nexo_t'][t]
        sigma_psi_init = setup.pars['sigma_psi_init']
        psi_couple = setup.exogrid.psi_t[t+1]
        
        
        if female:
            nz_single = setup.exogrid.zf_t[t].shape[0]
            p_mat = np.empty((nz_single,nexo))
            z_own = setup.exogrid.zf_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zm_t[t+1]
            zmat_own = setup.exogrid.zf_t_mat[t]
            pz_precomputed = self.partners_distribution_fem
        else:
            nz_single = setup.exogrid.zm_t[t].shape[0]
            p_mat = np.empty((nz_single,nexo))
            z_own = setup.exogrid.zm_t[t]
            n_zown = z_own.shape[0]
            z_partner = setup.exogrid.zf_t[t+1]
            zmat_own = setup.exogrid.zm_t_mat[t] 
            pz_precomputed = self.partners_distribution_mal
            
        def ind_conv(a,b,c): return setup.all_indices(t,(a,b,c))[0]
        
        pz_all = pz_precomputed['prob_z']
        
        
        pick = t if t < len(pz_all) else -1
        
        pz = pz_precomputed['prob_z'][pick]         
        pa = pz_precomputed['prob_a_by_z'][pick] 
        va = pz_precomputed['val_a_by_z'][pick] 
        
        
        na_matches = pa.shape[-1]
        
        p_mat_pa = np.empty((nexo,) + (na_matches,),dtype=setup.dtype)
        p_mat_ia = np.empty((nexo,) + (setup.na,) + (na_matches,),dtype=np.int16)
        
        agrid_s = self.agrid_s
        
        for iz in range(n_zown):
            p_psi = int_prob(psi_couple,mu=0.0,sig=sigma_psi_init)
            if female:
                p_zm  = np.array(pz)
                p_zf  = zmat_own[iz,:]
            else:
                p_zf  = np.array(pz)
                p_zm  = zmat_own[iz,:]
            
            
            p_vec = np.zeros(nexo)
            ie, izf, izm, ipsi = setup.all_indices(t)
            
            izpnext = izm if female else izf
            
            p_vec = p_zm[izm]*p_zf[izf]*p_psi[ipsi]
            assert np.allclose(p_vec.sum(),1.0)                            
            p_mat[iz,:] = p_vec
        
        a_resuling = np.clip(agrid_s[None,:,None] + va[:,None,:],0.0,self.agrid_c.max()-1e-5)
        ia_resulting = np.clip(np.searchsorted(setup.agrid_c,a_resuling) - 1,0,setup.na-1)
        p_resulting = pa
        assert np.allclose(pa.sum(axis=-1),1.0)
        # this slightly brings things down
        
        p_mat_ia[:,:,:] = ia_resulting[izpnext,:,:]
        p_mat_pa[:,:] = p_resulting[izpnext,:]
            
        nexo_ext = nexo*na_matches
        p_exo_ext = np.empty((nz_single,nexo_ext),dtype=self.dtype)
        ia_table = np.empty((setup.na,nexo_ext),dtype=np.int16) # resulting ia for each na and nexo_ext
        
        
        
        corresponding_iexo = np.empty(nexo_ext,dtype=np.int16)
        corresponding_imatch = np.empty(nexo_ext,dtype=np.int8)
        
        for im in range(na_matches):
            p_exo_ext[:,(im*nexo):((im+1)*nexo)] = p_mat*(p_mat_pa[:,im][None,:])
            ia_table[:,(im*nexo):((im+1)*nexo)] = p_mat_ia[:,:,im].T
            corresponding_iexo[(im*nexo):((im+1)*nexo)] = ie
            corresponding_imatch[(im*nexo):((im+1)*nexo)] = im
        
        
            
        return {'p_mat_iexo':p_mat,
                'p_mat_extended':p_exo_ext,
                'ia_c_table':ia_table,
                'corresponding_iexo':corresponding_iexo,
                'corresponding_imatch':corresponding_imatch}
    
        
        
    
    
    def build_matches(self):
        # for time moment and each position in grid for singles (ia,iz)
        # it computes probability distribution over potential matches
        # this is relevant for testing and simulations
        
        out_m = []
        out_f = []
        for t in range(self.pars['Tret'] - 1):
            out_m.append(self._mar_mats(t,female=False))
            out_f.append(self._mar_mats(t,female=True))
            
        self.matches_fem = out_f
        self.matches_mal = out_m
        
        
         
        
    
    
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
        #l = self.ls_levels['Couple and child'][il]
        #us = self.ls_ushift['Couple and child'][il]
        l, us = self._get_l_and_us('Couple and child',il)
        upub = self.u_pub(x,l) + ushift + psi
        return self.u(kf*c) + upub + us, self.u(km*c) + upub
    
    def u_couple_k(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta) 
        l, us = self._get_l_and_us('Couple and child',il)
        tus = theta*us
        return umult*self.u(c) + self.u_pub(x,l) + tus + ushift + psi
        
    def u_part_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        kf, km = self.c_mult(theta)   
        upub = ushift + psi
        return self.u(kf*c) + upub, self.u(km*c) + upub
    
    def u_couple_nk(self,c,x,il,theta,ushift,psi): # this returns utility of each partner out of some c
        umult = self.u_mult(theta)         
        return umult*self.u(c) + ushift + psi
    
    def u_single_k(self,c,x,il,ushift):
        umult = 1.0
        #l = self.ls_levels['Female and child'][il]
        #us = self.ls_ushift['Female and child'][il]
        l, us = self._get_l_and_us('Female and child',il)
        return umult*self.u(c) + self.u_pub(x,l) + us + ushift
    
    def _get_l_and_us(self,key,il):
        try:
            l = self.ls_levels[key][il]
            us = self.ls_ushift[key][il] 
        except:
            l = 0.0*il
            us = 0.0*il
            for i, (li, usi) in enumerate(zip(self.ls_levels[key],self.ls_ushift[key])):
                mask = (il==i)
                l += mask*li
                us += mask*usi
        return l, us

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
        p = self.pars['preg_a0'] + self.pars['preg_at']*t + \
            self.pars['preg_at2']*(t**2) + \
            self.pars['preg_az']*z + self.pars['preg_azt']*t*z
        #p = self.pars['preg_a0'] + self.pars['preg_at']*t + \
        #    self.pars['preg_az']*z + self.pars['preg_azt']*t*z
        return np.clip(p,0.0,1.0)
    
    def unplanned_pregnancy_probability(self):
        self.upp_precomputed_fem = list()
        for t in range(self.pars['T']):
            zf = self.exogrid.zf_t[t]
            pf = self._unplanned_pregnancy_probability_fun(t,zf)
            self.upp_precomputed_fem.append(pf)
            
        self.upp_precomputed_mal = list()
        for t in range(self.pars['T']):
            zm = self.exogrid.zm_t[t]
            pm = self._unplanned_pregnancy_probability_fun(t,zm)
            self.upp_precomputed_mal.append(pm)    
            
            
    def _tax_fun(self,lam,tau,avg_inc):
        def tax(income):
            return 1 - lam*((income/avg_inc)**(-tau))
        return tax
    
    def compute_taxes(self):
        self.taxes = dict()
        self.taxes['Female, single'] =      self._tax_fun(0.882,0.036,1.1086085)
        self.taxes['Male, single'] =        self._tax_fun(0.882,0.036,1.2123768)
        self.taxes['Female and child'] =    self._tax_fun(0.926,0.042,1.5605029) # use one child
        self.taxes['Couple, no children'] = self._tax_fun(0.903,0.058,2.3169997)
        self.taxes['Couple and child'] =    self._tax_fun(0.925,0.070,2.3169997)
        
    def compute_child_support_transitions(self,*,child_support_share):
        from interp_np import interp
        # this computes sets of transition matrices when in married couple
        # that breaks male gives away part of his income to wife. In this 
        # setup that means that productivity of male decreases and 
        # productivity of female increases.
        
        # structure: for each t there are indices of zf, probabilities of
        # zf, indices of zm, probabilities of zm
        
        p = self.pars
        transitions = list()
        
        for t in range(p['T']):
            trans_t = dict()
            
            
            
            
            trend_f = p['f_wage_trend'][t]
            trend_m = p['m_wage_trend'][t]
            zf_all = self.exogrid.zf_t[t]
            zm_all = self.exogrid.zm_t[t]
            
            if zf_all.size > 1 and zm_all.size>1:
                
                income_fem = np.exp(trend_f + zf_all)
                income_mal = np.exp(trend_m + zm_all)
                
                
                child_support = child_support_share*income_mal
                income_fem_adj = income_fem[:,None] + child_support[None,:]
                income_mal_adj = income_mal - child_support
                
                z_fem_adj = np.log(income_fem_adj) - trend_f
                z_mal_adj = np.log(income_mal_adj) - trend_m
                
                #von_fem = VecOnGrid(zf_all,z_fem_adj)
                trans_t['i_this_fem'], trans_t['w_this_fem'] = \
                    interp(zf_all,z_fem_adj,return_wnext=False,trim=True)
                    
                
                #von_mal = VecOnGrid(zm_all,z_mal_adj)
                trans_t['i_this_mal'], trans_t['w_this_mal'] = \
                    interp(zm_all,z_mal_adj,return_wnext=False,trim=True)
                    
            else:
                fill = np.array([-1],dtype=np.int16), np.array([0.0],dtype=self.dtype)
                trans_t['i_this_fem'], trans_t['w_this_fem'] = fill
                trans_t['i_this_mal'], trans_t['w_this_mal'] = fill
                
            
            transitions.append(trans_t)
        
        self.child_support_transitions = transitions
    
    def cupyfy(self):
        # this sends copies of important objects to cupy
        # only things that are going to be reused are worth storing
        
        
        if not cupy:
            self.cupy = None
            return
        
        stuff = dict()
        
        stuff['theta_orig_on_fine'] = cp.array(self.theta_orig_on_fine,dtype=self.dtype)
        stuff['thetagrid'] = cp.array(self.thetagrid,dtype=self.dtype)
        stuff['thetagrid_fine'] = cp.array(self.thetagrid_fine,dtype=self.dtype)
        stuff['v_thetagrid_fine'] = VecOnGrid(stuff['thetagrid'],stuff['thetagrid_fine'],force_cupy=True) 
        stuff['agrid_c'] = cp.array(self.agrid_c,dtype=self.dtype)
        stuff['agrid_s'] = cp.array(self.agrid_s,dtype=self.dtype)
        stuff['sgrid_c'] = cp.array(self.sgrid_c,dtype=self.dtype)
        stuff['sgrid_s'] = cp.array(self.sgrid_s,dtype=self.dtype)
        stuff['vsgrid_c'] = VecOnGrid(self.agrid_c,self.sgrid_c,force_cupy=True)
        stuff['vsgrid_s'] = VecOnGrid(self.agrid_s,self.sgrid_s,force_cupy=True)
        stuff['mgrid_c'] = cp.array(self.mgrid_c,dtype=self.dtype)
        stuff['mgrid_s'] = cp.array(self.mgrid_s,dtype=self.dtype)
        
        uu = self.u_precomputed
        upc = {key : {e: cp.array(uu[key][e]) for e in uu[key]} for key in uu}
        stuff['u_precomputed'] = upc
        
        self.cupy = namedtuple('cupy',stuff.keys())(**stuff)
        
        
        
        
    

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
                 z_change_f=0.0,z_change_m=0.0,
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
        self.z_change_f = z_change_f
        self.z_change_m = z_change_m
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

