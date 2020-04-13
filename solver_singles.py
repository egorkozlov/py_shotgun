#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
#from scipy.optimize import fminbound

#from opt_test import build_s_grid, sgrid_on_agrid, get_EVM

from optimizers import v_optimize_couple

def v_iter_single(setup,t,EV,female,ushift):
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    
    
    dtype = setup.dtype
    
    
    zvals = setup.exogrid.zf_t[t] if female else setup.exogrid.zm_t[t]
    ztrend = setup.pars['f_wage_trend'][t] if female else setup.pars['m_wage_trend'][t]
    beta = setup.pars['beta_t'][t]
    R = setup.pars['R_t'][t]
    
    
    money_t = (R*agrid_s,np.exp(zvals + ztrend),np.zeros_like(zvals))
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size),dtype=setup.dtype)    
    
    sname = 'Female, single' if female else 'Male, single'
    
    uu = setup.u_precomputed[sname]['u']
    ux = setup.u_precomputed[sname]['x']
    ls = setup.ls_levels[sname]
    
    
    V, c, x, s, i_opt, ils, V_all_l = \
        v_optimize_couple(money_t,sgrid_s,(setup.vsgrid_s,EV[...,None,None]),setup.mgrid_s,uu,ux,
                              ls,beta,ushift,dtype=dtype)
    
    
    assert np.all(x<1e-5) # no x and no labor supply
    
    V_ret, c_opt, s_opt, i_opt = (q.squeeze(axis=2) for q in (V, c, s, i_opt))
    
    
    # reconstruct V for precision
    
    c_refine = R*agrid_s[:,None] + np.exp(zvals + ztrend)[None,:] - s_opt
     
    EVexp = setup.vsgrid_s.apply_preserve_shape(EV)
    V_refine = setup.u(c_refine) + ushift + beta*np.take_along_axis(EVexp,i_opt,0)
    #print('after refinement V difference is {}'.format(np.abs(np.max(V_ret-V_refine))))
    
    
    def r(x): return x.astype(dtype,copy=False)
    
    return r(V_refine), r(c_opt), r(s_opt)#, s_opt/money




def v_iter_single_mom(setup,t,EV,ushift):
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    ls = setup.ls_levels['Female and child']
    
    
    
    dtype = setup.dtype
        
    
    zvals = setup.exogrid.zf_t[t]
    ztrend = setup.pars['f_wage_trend'][t]
    R = setup.pars['R_t'][t]
    beta = setup.pars['beta_t'][t]
    
    wm0 = np.zeros_like(zvals + ztrend)
    money_t = (R*agrid_s,np.exp(zvals + ztrend),wm0)
    
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size,ls.size),dtype=setup.dtype)
    
    # I add virtual theta axis
    EV_t = (setup.vsgrid_s,EV[...,None,:])
    
    uu = setup.u_precomputed['Female and child']['u']
    ux = setup.u_precomputed['Female and child']['x']
    
    V, c, x, s, i_opt, ils, V_all_l = \
        v_optimize_couple(money_t,sgrid_s,EV_t,setup.mgrid_s,uu,ux,
                              ls,beta,ushift,dtype=dtype)
    
    
    # remove the virtual axis
    V_old, c, x, s, i_opt, ils, V_all_l = (q.squeeze(axis=2) for q in (V, c, x, s, i_opt, ils, V_all_l))
    
    
    u = setup.u_single_k(c,x,ils,ushift)
    
    
    EVexp = setup.vsgrid_s.apply_preserve_shape(EV)    
    V_refine = u + beta*np.take_along_axis(np.take_along_axis(EVexp,i_opt[...,None],0),ils[...,None],2).squeeze(axis=2)
    #print('After retirement (sm) max diff is {}'.format(np.max(np.abs(V_old-V_refine))))
    
    def r(x): return x.astype(dtype,copy=False)
    
    
    return r(V_refine), r(c), r(x), r(s), r(ils), r(V_all_l)
