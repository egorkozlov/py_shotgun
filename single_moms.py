#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:51:35 2020

@author: egorkozlov
"""


import numpy as onp
try:
    import cupy as cp
    #import cupyx
    gpu = True
except:
    gpu = False


from optimizers import v_optimize_couple

nogpu = (not gpu)

np = onp if not gpu else cp
    

from singles import co



def solve_single_mom(model,t,Vnext,ushift,verbose=False):
    setup = model.setup
    
    EV, dec = ev_single_k(setup,Vnext,t) # Vnext can be None
    model.time('Integration (s)')
    
    V_tuple = v_iter_single_mom(setup,t,EV,ushift)
    model.time('Optimization (s)')
    
    return V_tuple, dec





def v_iter_single_mom(setup,t,EV,ushift):
    
    agrid_s = setup.agrid_s if nogpu else setup.cupy.agrid_s
    sgrid_s = setup.sgrid_s if nogpu else setup.cupy.sgrid_s
    ls = setup.ls_levels['Female and child']
    
    
    taxfun = setup.taxes['Female and child']
    dtype = setup.dtype
        
    
    zvals = co(setup.exogrid.zf_t[t])
    ztrend = setup.pars['f_wage_trend'][t]
    R = setup.pars['R_t'][t]
    beta = setup.pars['beta_t'][t]
    
    wm0 = np.zeros_like(zvals + ztrend)
    money_t = (R*agrid_s,np.exp(zvals + ztrend),wm0)
    
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size,ls.size),dtype=setup.dtype)
    
    
    vs = setup.vsgrid_s if nogpu else setup.cupy.vsgrid_s
    # I add virtual theta axis
    EV_t = (vs, EV[...,None,:])
    
    uu = setup.u_precomputed['Female and child']['u'] if nogpu else setup.cupy.u_precomputed['Female and child']['u']
    ux = setup.u_precomputed['Female and child']['x'] if nogpu else setup.cupy.u_precomputed['Female and child']['x'] 
    
    V, c, x, s, i_opt, ils, V_all_l = \
        v_optimize_couple(money_t,sgrid_s,EV_t,setup.mgrid_s,uu,ux,
                              ls,beta,ushift,taxfun=taxfun,dtype=dtype)
    
    
    # remove the virtual axis
    V_old, c, x, s, i_opt, ils, V_all_l = (q.squeeze(axis=2) for q in (V, c, x, s, i_opt, ils, V_all_l))
    
    
    u = setup.u_single_k(c,x,ils,ushift)
    
    
    EVexp = vs.apply_preserve_shape(EV)    
    V_refine = u + beta*np.take_along_axis(np.take_along_axis(EVexp,i_opt[...,None],0),ils[...,None],2).squeeze(axis=2)
    #print('After retirement (sm) max diff is {}'.format(np.max(np.abs(V_old-V_refine))))
    
    def r(x): return x.astype(dtype,copy=False)
    
    
    return r(V_refine), r(c), r(x), r(s), r(ils)


from singles import ev_single_meet 

def ev_single_k(setup,V,t):
    # behave as if meet & pregnant
    # expected value of single person meeting a partner with a chance pmeet
    if V is None: return None, {}
    
    pmeet = setup.pars['pmeet_t'][t]*setup.pars['pmeet_multiplier_fem']
    
    female = True
    
    skip_mar = (pmeet < 1e-5)
        
    nz = setup.pars['n_zf_t'][t] if female else setup.pars['n_zm_t'][t]
    nl = len(setup.ls_levels['Female and child'])
    
    # FIXME: when people meet their skill depreciation stops for one period.
    # This may be minor but is a bit inconsistent
    
    EV_stay = np.zeros((setup.na,nz) + (nl,),dtype=setup.dtype)
    
    for il in range(nl):
         M = co(setup.exogrid.zf_t_mat_by_l_sk[il][t])      
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
    
