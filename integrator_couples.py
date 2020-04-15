#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is integrator for couples

"""

import numpy as np
#from renegotiation import v_last_period_renegotiated, v_renegotiated_loop
    
from renegotiation_unilateral import v_ren_uni, v_no_ren

def ev_couple_m_c(setup,Vpostren,t,haschild,use_sparse=True):
    # computes expected value of couple entering the next period with an option
    # to renegotiate or to break up
    
    canswitch = setup.pars['is fertile'][t]
    can_divorce = setup.pars['can divorce'][t] # !! no divorce = no fertility here
    
    
    if can_divorce:
        out = v_ren_uni(setup,Vpostren,haschild,canswitch,t)
        
    else:
        out = v_no_ren(setup,Vpostren,haschild,canswitch,t)
        
        
        
        
    _Vren2 = out.pop('Values') 
    #_Vren2=out['Values']
    dec = out
    
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine]
    
    Vren = {'M':{'VR':tk(_Vren2[0]),'VC':tk(_Vren2[1]), 'VF':tk(_Vren2[2]),'VM':tk(_Vren2[3])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}
    
    '''
    # this is false: instead of V0.min should be minimum of many things
    V0 = Vpostren['Couple and child']['V']
    V1 = _Vren2[1]
    
    if haschild:
        assert V1.max() <= V0.max() + 1e-4
        assert V1.min() >= V0.min() - 1e-4
        
    '''
    
        
    # accounts for exogenous transitions
    
    EVr, EVc, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,haschild,use_sparse,down=False)
    
    
    return (EVr, EVc, EVf, EVm), dec


def ev_couple_exo(setup,Vren,t,haschild,use_sparse=True,down=False):
    
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    
    def mmult(a,b):
        if use_sparse:
            return a*b
        else:
            return np.dot(a,b.T)
        
        
    if haschild:
        mat_sp = setup.exogrid.all_t_mat_by_l_spt_k
        mat_nsp = setup.exogrid.all_t_mat_by_l_k
    else:
        mat_sp = setup.exogrid.all_t_mat_by_l_spt_nk
        mat_nsp = setup.exogrid.all_t_mat_by_l_nk
        
    
    nl = len(mat_sp)
    
    na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta 
    
    
    Vr, Vc, Vf, Vm = Vren['VR'], Vren['VC'], Vren['VF'], Vren['VM']
    EVr, EVc, EVf, EVm = np.zeros((4,na,nexo,ntheta,nl),dtype=setup.dtype)
    
    
    for il in range(nl):
        
        M = mat_sp[il][t] if use_sparse else mat_nsp[il][t]
        
        
        
        for itheta in range(ntheta):
            EVr[...,itheta,il]  = mmult(Vr[...,itheta],M)
            EVc[...,itheta,il]  = mmult(Vc[...,itheta],M)
            EVf[...,itheta,il] = mmult(Vf[...,itheta],M)             
            EVm[...,itheta,il] = mmult(Vm[...,itheta],M)             
            

    #assert not np.allclose( EV[...,0], EV[...,1])
    
    
    return EVr, EVc, EVf, EVm