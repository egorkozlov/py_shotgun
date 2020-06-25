#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np
from timeit import default_timer


from optimizers import v_optimize_couple

from platform import system

if system() != 'Darwin' and system() != 'Windows':    
    nbatch_def = 500
    use_cp = True
    
elif system() == 'Windows':
    
    nbatch_def = 17
    use_cp = True
    
else:
    
    nbatch_def = 17
    use_cp = False

def v_iter_couple(setup,t,EV_tuple,ushift,haschild,nbatch=nbatch_def,verbose=False):
    
    if verbose: start = default_timer()
    
    agrid = setup.agrid_c
    sgrid = setup.sgrid_c
    
    dtype = setup.dtype
    
    
    
    
    
    key = 'Couple and child' if haschild else 'Couple, no children'
    
    
    taxfun = setup.taxes[key] if t < setup.pars['Tret'] else lambda x : 0.0*x
    
    ls = setup.ls_levels[key]
    uu, ux = setup.u_precomputed[key]['u'],setup.u_precomputed[key]['x']
    upart, ucouple = (setup.u_part_k, setup.u_couple_k) if haschild else (setup.u_part_nk, setup.u_couple_nk)

    nls = len(ls)
    
    
    # type conversion is here
    
    zf  = setup.exogrid.all_t[t][:,0]
    zm  = setup.exogrid.all_t[t][:,1]
    zftrend = setup.pars['f_wage_trend'][t]
    zmtrend = setup.pars['m_wage_trend'][t]

    psi = setup.exogrid.all_t[t][:,2]
    beta = setup.pars['beta_t'][t]
    sigma = setup.pars['crra_power']
    R = setup.pars['R_t'][t]

    
    
    wf = np.exp(zf + zftrend)
    wm = np.exp(zm + zmtrend)
    
    #labor_income = np.exp(zf) + np.exp(zm)
    
    #money = R*agrid[:,None] + wf[None,:] 
    
    
    nexo = setup.pars['nexo_t'][t]
    shp = (setup.na,nexo,setup.ntheta)
    
    if EV_tuple is None:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = np.zeros((4,) + shp + (nls,),dtype=setup.dtype)
    else:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = EV_tuple    
    
    # type conversion
    sgrid,sigma,beta = (dtype(x) for x in (sgrid,sigma,beta))
    
    V_couple, c_opt, s_opt, x_opt = np.empty((4,)+shp,dtype)
    i_opt, il_opt = np.empty(shp,np.int16), np.empty(shp,np.int16)
    
    V_all_l = np.empty(shp+(nls,),dtype=dtype)
    
    theta_val = dtype(setup.thetagrid)
    
    # the original problem is max{umult*u(c) + beta*EV}
    # we need to rescale the problem to max{u(c) + beta*EV_resc}
    
    istart = 0
    ifinish = nbatch if nbatch < nexo else nexo
    
    # this natually splits everything onto slices
    
    for ibatch in range(int(np.ceil(nexo/nbatch))):
        #money_i = money[:,istart:ifinish]
        assert ifinish > istart
        
        money_t = (R*agrid, wf[istart:ifinish], wm[istart:ifinish])
        EV_t = (setup.vsgrid_c,EVr_by_l[:,istart:ifinish,:,:])
        
        
        
        
        
        V_pure_i, c_opt_i, x_opt_i, s_opt_i, i_opt_i, il_opt_i, V_all_l_i = \
           v_optimize_couple(money_t,sgrid,EV_t,setup.mgrid_c,uu,ux,
                                 ls,beta,ushift,taxfun=taxfun,dtype=dtype)
           
        V_ret_i = V_pure_i + psi[None,istart:ifinish,None]
        
        
        
        
        
        V_couple[:,istart:ifinish,:] = V_ret_i # this estimate of V can be improved
        c_opt[:,istart:ifinish,:] = c_opt_i 
        s_opt[:,istart:ifinish,:] = s_opt_i
        i_opt[:,istart:ifinish,:] = i_opt_i
        x_opt[:,istart:ifinish,:] = x_opt_i
        il_opt[:,istart:ifinish,:] = il_opt_i
        V_all_l[:,istart:ifinish,:,:] = V_all_l_i # we need this for l choice so it is ok
        
        
        istart = ifinish
        ifinish = ifinish+nbatch if ifinish+nbatch < nexo else nexo
        
        if verbose: print('Batch {} done at {} sec'.format(ibatch,default_timer()-start))
    
    
    assert np.all(c_opt > 0)
    
    psi_r = psi[None,:,None].astype(setup.dtype,copy=False)
    
    # finally obtain value functions of partners
    uf, um = upart(c_opt,x_opt,il_opt,theta_val[None,None,:],ushift,psi_r)
    uc = ucouple(c_opt,x_opt,il_opt,theta_val[None,None,:],ushift,psi_r)
    
    
    EVf_all, EVm_all, EVc_all  = (setup.vsgrid_c.apply_preserve_shape(x) for x in (EV_fem_by_l, EV_mal_by_l,EVc_by_l))
    
    V_fem = uf + beta*np.take_along_axis(np.take_along_axis(EVf_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_mal = um + beta*np.take_along_axis(np.take_along_axis(EVm_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_all = uc + beta*np.take_along_axis(np.take_along_axis(EVc_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    def r(x): return x
    
    '''
    if V_all.shape[-1] > 1:
        V_weighted = V_fem*setup.thetagrid[None,None,:] + V_mal*(1-setup.thetagrid[None,None,:])
        vd_max = np.max(np.abs(V_weighted - V_all))
        vd_mean = np.mean(np.abs(V_weighted - V_all))
        print('max diff is {}, mean diff is {}'.format(vd_max,vd_mean))
    '''
    assert V_all.dtype==EVc_all.dtype==V_couple.dtype
    
    
    '''
    if uf.shape[-1] > 1:
        assert np.all(np.diff(uf,axis=-1)>0), 'Monotonicity is violated for females in u'
        assert np.all(np.diff(um,axis=-1)<0), 'Monotonicity is violated for males in u'
    
    if V_fem.shape[-1] > 1:
        assert np.all(np.diff(V_fem,axis=-1)>0), 'Monotonicity is violated for females in V'
        assert np.all(np.diff(V_mal,axis=-1)<0), 'Monotonicity is violated for males in V'
        
    '''
    return r(V_all), r(V_fem), r(V_mal), r(c_opt), r(x_opt), r(s_opt), il_opt, r(V_all_l)




