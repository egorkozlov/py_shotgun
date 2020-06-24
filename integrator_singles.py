#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
#import dill as pickle

from marriage import v_mar_igrid, v_no_mar
    



def ev_single(setup,V,sown,female,t,trim_lvl=0.001):
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet_t'][t]
    
    skip_mar = (pmeet < 1e-5)
    
    
    if female:
        cangiveabirth = setup.pars['is fertile'][t]
    else:
        cangiveabirth = setup.pars['is fertile'][t-2] if t>=2 else False
    
    
    EV_meet_np, dec_np = ev_single_meet(setup,V,sown,female,t,
                                  skip_mar=skip_mar,trim_lvl=trim_lvl,
                                      unplanned_preg=False)
    
    
    
    
    
    if cangiveabirth:
        
        ppreg0 = setup.upp_precomputed_fem[t] if female else setup.upp_precomputed_mal[t]
        
        ppreg = ppreg0[None,:]
        
        EV_meet_p, dec_p = ev_single_meet(setup,V,sown,female,t,
                                      skip_mar=skip_mar,trim_lvl=trim_lvl,
                                          unplanned_preg=True)
        
        EV_meet = EV_meet_np*(1-ppreg) + EV_meet_p*ppreg
    else:
        dec_p = dec_np
        EV_meet = EV_meet_np
    
    
    
    
    if female:
        M = setup.exogrid.zf_t_mat[t].T
        EV_nomeet =  np.dot(V['Female, single']['V'],M)
    else:
        M = setup.exogrid.zm_t_mat[t].T
        EV_nomeet =  np.dot(V['Male, single']['V'],M)
    
    
    dec = {'Not pregnant':dec_np, 'Pregnant':dec_p}
    
    
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet, dec



def ev_single_k(setup,V,sown,t,trim_lvl=0.001):
    # behave as if meet & pregnant
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet_t'][t]*setup.pars['pmeet_multiplier_fem']
    
    female = True
    
    skip_mar = (pmeet < 1e-5)

    EV_meet, dec = ev_single_meet(setup,V,sown,female,t,
                                  skip_mar=skip_mar,trim_lvl=trim_lvl,
                                      unplanned_preg=True,single_mom=True)
        
    
    nl = len(setup.ls_levels['Female and child'])
    
    # FIXME: when people meet their skill depreciation stops for one period.
    # This may be minor but is a bit inconsistent
    
    EV_stay = np.zeros(EV_meet.shape + (nl,),dtype=setup.dtype)
    
    for il in range(nl):
         M = setup.exogrid.zf_t_mat_by_l_sk[il][t]         
         EV_stay[...,il] = np.dot(V['Female and child']['V'],M.T)
    
    
    dec = {'Not pregnant':dec, 'Pregnant':dec}
    
    return (1-pmeet)*EV_stay + pmeet*EV_meet[...,None], dec
    


def ev_single_meet(setup,V,sown,female,t,skip_mar=False,
                   unplanned_preg=False,single_mom=False,trim_lvl=0.001):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
    # single_mom emulates unplanned pregnancy but add utility punishments
    
    nexo = setup.pars['nexo_t'][t]
    if female:
        cangiveabirth = setup.pars['is fertile'][t]
    else:
        cangiveabirth = setup.pars['is fertile'][t-2] if t>=2 else False
    ns = sown.size
    no_kids_at_meet = setup.pars['no kids at meeting']
    
    
    uloss_fem = setup.pars['disutil_marry_sm_fem'] if single_mom else 0.0
    uloss_mal = setup.pars['disutil_marry_sm_mal'] if single_mom else 0.0
    
    
    
    uloss_fem_single = 0.0
    uloss_mal_single = 0.0
    
        
    V_next = np.ones((ns,nexo))*(-1e20)
    
    
    EV = 0.0
    
    # this is a shallow copy: it does not actually copy big matrices
    if (female and not unplanned_preg):# or (female and single_mom):
        matches = setup.matches['Female meets male, no upp'][t].copy()
    elif female and unplanned_preg:
        matches = setup.matches['Female meets male, upp'][t].copy()
    elif female and single_mom:
        matches = setup.matches['Single mother meets male'][t].copy()
    else:
        matches = setup.matches['Male meets female, no upp'][t].copy()
    
        
        
    p_mat = matches['iexo_matrix'].T # !!! transpose here 
    i_assets_c, p_assets_c = matches['i_a_mat'], matches['p_a_mat']
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    npart = i_assets_c.shape[1]
    
    
    dec = np.zeros(matches['iexo'].shape,dtype=np.bool)
    abn = np.zeros(matches['iexo'].shape,dtype=np.bool)
    morc = np.zeros(matches['iexo'].shape,dtype=np.bool)
    tht = -1*np.ones(matches['iexo'].shape,dtype=np.int32)
    iconv = matches['iconv']
    
    for i in range(npart):
        
        
        
        if not skip_mar:
            if not unplanned_preg and not single_mom:
                # compare whether to give or not to give a birth
                res_c = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                         female=female,giveabirth=False,
                                         uloss_fem=uloss_fem,uloss_mal=uloss_mal,
                                         uloss_fem_single=uloss_fem_single,
                                         uloss_mal_single=uloss_mal_single)
                
                if cangiveabirth and (not no_kids_at_meet):
                    # maybe cannot give a birth at all
                    res_m = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                         female=female,giveabirth=True,
                                         uloss_fem=uloss_fem,uloss_mal=uloss_mal,
                                         uloss_fem_single=uloss_fem_single,
                                         uloss_mal_single=uloss_mal_single)
                else:            
                    res_m = res_c
            else:
                # no choices
                if unplanned_preg and (not single_mom): assert cangiveabirth
                res_m = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                    unplanned_pregnancy=True,
                                         female=female,giveabirth=True,
                                         uloss_fem=uloss_fem,uloss_mal=uloss_mal,
                                         uloss_fem_single=uloss_fem_single,
                                         uloss_mal_single=uloss_mal_single)
                res_c = res_m
                
        else:
            # skip
            res_c = v_no_mar(setup,t,V,i_assets_c[:,i],inds,
                                     female=female,giveabirth=False)
            res_m = res_c
            
        
        (vfoutc, vmoutc), nprc, decc, thtc, abnc =  res_c['Values'], res_c['NBS'], res_c['Decision'], res_c['theta'], res_c['Abortion']
        (vfoutm,vmoutm), nprm, decm, thtm, abnm = res_m['Values'], res_m['NBS'], res_m['Decision'], res_m['theta'], res_m['Abortion']
        
        # choice is made based on Nash Surplus value
        i_birth = (nprm>nprc) 
        
        if not cangiveabirth: assert not np.any(i_birth)
        
        if female:
            vout = i_birth*vfoutm + (1-i_birth)*vfoutc
        else:
            vout = i_birth*vmoutm + (1-i_birth)*vmoutc
            
        dec[:,:,iconv[:,i]] = (i_birth*decm + (1-i_birth)*decc)[:,None,:]
        tht[:,:,iconv[:,i]] = (i_birth*thtm + (1-i_birth)*thtc)[:,None,:]
        abn[:,:,iconv[:,i]] = (i_birth*abnm + (1-i_birth)*abnc)[:,None,:]
        
        if not unplanned_preg and not single_mom:
            morc[:,:,iconv[:,i]] = i_birth[:,None,:]
        else:
            morc[:,:,iconv[:,i]] = True
        
        V_next[:,inds] = vout
        
        EV += (p_assets_c[:,i][:,None])*np.dot(V_next,p_mat)

    mout = matches#.copy() # this is here for reason. otherwise things get inserted into
    mout['Decision'] = dec
    mout['Child immediately'] = morc
    mout['theta'] = tht
    
    mout['Abortion'] = abn if female and unplanned_preg else False*abn
    
    return EV, mout