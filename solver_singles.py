#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
#from scipy.optimize import fminbound

#from opt_test import build_s_grid, sgrid_on_agrid, get_EVM
from optimizers import v_optimize_single



def v_iter_single(setup,t,EV,female,ushift):
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    
    
    dtype = setup.dtype
    
    ind, p = setup.vsgrid_s.i, setup.vsgrid_s.wthis
    
    
    zvals = setup.exogrid.zf_t[t] if female else setup.exogrid.zm_t[t]
    ztrend = setup.pars['f_wage_trend'][t] if female else setup.pars['m_wage_trend'][t]
    sigma = setup.pars['crra_power']
    beta = setup.pars['beta_t'][t]
    R = setup.pars['R_t'][t]
    
    
    money_t = (R*agrid_s,np.exp(zvals + ztrend))
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size),dtype=setup.dtype)    
    
    V_ret, c_opt, s_opt = v_optimize_single(money_t,sgrid_s,(ind,p,EV),sigma,beta,ushift,dtype=dtype)
    
    def r(x): return x.astype(dtype)
    
    return r(V_ret), r(c_opt), r(s_opt)#, s_opt/money



from optimizers import v_optimize_couple

def v_iter_single_mom(setup,t,EV,ushift):
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    ls = setup.ls_levels_sk
    
    
    
    dtype = setup.dtype
    
    ind, p = setup.vsgrid_s.i, setup.vsgrid_s.wthis
    
    
    zvals = setup.exogrid.zf_t[t]
    ztrend = setup.pars['f_wage_trend'][t]
    R = setup.pars['R_t'][t]
    beta = setup.pars['beta_t'][t]
    
    wm0 = np.zeros_like(zvals + ztrend)
    money_t = (R*agrid_s,np.exp(zvals + ztrend),wm0)
    
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size,ls.size),dtype=setup.dtype)
    
    # I add virtual theta axis
    EV_t = (ind,p,EV[...,None,:])
    
    uu, ux = setup.ucouple_precomputed_u_sk[:,None,:],setup.ucouple_precomputed_x_sk[:,None,:]
    
    V, c, x, s, i_opt, ils, V_all_l = \
        v_optimize_couple(money_t,sgrid_s,EV_t,setup.mgrid_s,uu,ux,
                              ls,beta,ushift,dtype=dtype)
       
    # remove the virtual axis
    def r(x): return x.squeeze(axis=2).astype(dtype)
    
    return r(V), r(c), r(x), r(s), r(ils), r(V_all_l)
