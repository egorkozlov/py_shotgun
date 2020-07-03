#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
#import dill as pickle

from marriage import v_mar
    



def ev_single(setup,V,female,t):
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet_t'][t]
    
    skip_mar = (pmeet < 1e-5)
    
    if female:
        M = setup.exogrid.zf_t_mat[t].T
        EV_nomeet =  np.dot(V['Female, single']['V'],M)
    else:
        M = setup.exogrid.zm_t_mat[t].T
        EV_nomeet =  np.dot(V['Male, single']['V'],M)
        
    
    if skip_mar:
        return EV_nomeet, {}
        
        
    # else do the meeting thing
    
    ti = t if female else (t-2 if t>=2 else -1)
    upp_possible = setup.pars['is fertile'][ti]
    
    EV_meet_np, dec_np = ev_single_meet(setup,V,female,t,match_type='Regular')
    # possible choice to have children immediately is done inside
    
    if upp_possible:
        
        ppreg0 = setup.upp_precomputed_fem[t] if female else setup.upp_precomputed_mal[t]
        ppreg = ppreg0[None,:]
        EV_meet_p, dec_p = ev_single_meet(setup,V,female,t,match_type='Unplanned pregnancy')
        EV_meet = EV_meet_np*(1-ppreg) + EV_meet_p*ppreg
    
    else:
        
        dec_p = dec_np
        EV_meet = EV_meet_np
    
    
    
    dec = {'Not pregnant':dec_np, 'Pregnant':dec_p}


        
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet, dec



def ev_single_k(setup,V,t):
    # behave as if meet & pregnant
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet_t'][t]*setup.pars['pmeet_multiplier_fem']
    
    female = True
    
    skip_mar = (pmeet < 1e-5)
        
    nz = setup.pars['n_zf_t'][t] if female else setup.pars['n_zm_t'][t]
    nl = len(setup.ls_levels['Female and child'])
    
    # FIXME: when people meet their skill depreciation stops for one period.
    # This may be minor but is a bit inconsistent
    
    EV_stay = np.zeros((setup.na,nz) + (nl,),dtype=setup.dtype)
    
    for il in range(nl):
         M = setup.exogrid.zf_t_mat_by_l_sk[il][t]         
         EV_stay[...,il] = np.dot(V['Female and child']['V'],M.T)
    
    
    if not skip_mar:
        EV_meet, dec = ev_single_meet(setup,V,female,t,match_type='Single mother')
    else:
        EV_meet = EV_stay
        dec = {}
    
    dec = {'Not pregnant':dec, 'Pregnant':dec}
    
    if not skip_mar:
        return (1-pmeet)*EV_stay + pmeet*EV_meet[...,None], dec
    else:
        return EV_stay, {}
    


def ev_single_meet(setup,V,female,t,*,match_type):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
    
    
        
    try:
        matches = setup.matches_fem[t] if female else setup.matches_mal[t]
    except:
        matches = setup.matches_fem[-1] if female else setup.matches_mal[-1]
    
    #p_mat_iexo = matches['p_mat_iexo']
    p_mat_ext  = matches['p_mat_extended']
    icouple = matches['ia_c_table']
    inds_all = matches['corresponding_iexo']
    

    
    
    
    out = v_mar(setup,V,t,icouple,inds_all,match_type=match_type,female=female)
    
    
    
    
    if (match_type=='Single mother' or match_type=='Unplanned pregnancy'):
        pick_k = np.ones(out['Agree'].shape,dtype=np.bool_)
    else:
        pick_k = np.zeros(out['Agree'].shape,dtype=np.bool_)
    
    
    if match_type=='Regular' and (not setup.pars['no kids at meeting']):
        # additional choice to have kids immediately
        out_k = v_mar(setup,V,t,icouple,inds_all,match_type='Child immediately',female=female)
        pick_k = out_k['NBS'] > out['NBS'] # > is important as 0 are possible
        out_nk = out.copy()
        out = {key: np.array(out_k[key]*pick_k + out_nk[key]*(~pick_k),dtype=out_k[key].dtype) for key in out_nk}
    
    # parse output
        
    V_out = out['V_fem'] if female else out['V_mal']
    EV = np.dot( V_out, p_mat_ext.T )
    
    
    

      # this is here for reason. otherwise things get inserted into
    mout = matches.copy()
    
    mout['Decision'] = out['Agree']
    mout['Child immediately'] = pick_k
    mout['itheta'] = out['itheta']

    abn = out['Abortion']
    mout['Abortion'] = abn if female and match_type=='Unplanned pregnancy' else (False & abn)
    
    return EV, mout
