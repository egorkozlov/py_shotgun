#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 17:03:01 2020

@author: egorkozlov
"""

# this is a joined routine obitaining things for couples


import numpy as onp
from timeit import default_timer

try:
    import cupy as cp
    #import cupyx
    gpu = True
except:
    gpu = False


from optimizers import v_optimize_couple
    
    
def solve_couples(model,t,Vnext,ushift,haschild,verbose=False):
    setup = model.setup
    EV_tuple, dec = ev_couple_m_c(model,Vnext,t,haschild,use_sparse=True)
    model.time('Integration (c)')
    V_tuple = v_iter_couple(model,t,EV_tuple,ushift,haschild)
    model.time('Optimization (c)')
    return V_tuple, dec
    


from renegotiation_unilateral import v_ren_uni, v_no_ren

def ev_couple_m_c(model,Vpostren,t,haschild,use_sparse=True):
    # computes expected value of couple entering the next period with an option
    # to renegotiate or to break up
    
    setup = model.setup
    
    if Vpostren is None: return None, {}
    
    canswitch = setup.pars['is fertile'][t]
    can_divorce = setup.pars['can divorce'][t] # !! no divorce = no fertility here
    
    
    if can_divorce:
        out = v_ren_uni(model,Vpostren,haschild,canswitch,t)
    else:
        out = v_no_ren(model,Vpostren,haschild,canswitch,t)
        
        
        
        
    _Vren2 = out.pop('Values') 
    #_Vren2=out['Values']
    dec = out
    
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine]
    
    Vren = {'M':{'VR':tk(_Vren2[0]),'VC':tk(_Vren2[1]), 'VF':tk(_Vren2[2]),'VM':tk(_Vren2[3])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}
    
        
    # accounts for exogenous transitions
    
    EVr, EVc, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,haschild,use_sparse,down=False)
    
    
    return (EVr, EVc, EVf, EVm), dec


def ev_couple_exo(setup,Vren,t,haschild,use_sparse=True,down=False):

    if gpu:
        np = cp
        use_sparse = False # cupy-sparse is bad somehow
    else:
        np = onp
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    
    def mmult(a,b):
        if use_sparse:
            return a*b
        else:
            return np.dot(a,b.T)
        
        
    if haschild:
        if use_sparse:
            mat = setup.exogrid.all_t_mat_by_l_spt_k
        elif gpu:
            mat = setup.exogrid.all_t_mat_by_l_k
    else:
        if use_sparse:
            mat = setup.exogrid.all_t_mat_by_l_spt_nk
        else:
            mat = setup.exogrid.all_t_mat_by_l_nk
        
        
    
    
    nl = len(mat)
    
    
    
    na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta 
    
    
    Vr, Vc, Vf, Vm = Vren['VR'], Vren['VC'], Vren['VF'], Vren['VM']
    EVr, EVc, EVf, EVm = np.empty((4,na,nexo,ntheta,nl),dtype=setup.dtype)
    
    
    for il in range(nl):
        
        M = np.array(mat[il][t],copy=False) if gpu else mat[il][t]
        # when np is cp this puts stuff onto device
        # there is no point in doing it in advance unless you compute all of
        # them on gpu
        
        for itheta in range(ntheta):            
            EVr[...,itheta,il]  = mmult(Vr[...,itheta],M)
            EVc[...,itheta,il]  = mmult(Vc[...,itheta],M)
            EVf[...,itheta,il] = mmult(Vf[...,itheta],M)             
            EVm[...,itheta,il] = mmult(Vm[...,itheta],M)             
            
        del M
        
    
    
    return EVr, EVc, EVf, EVm





def v_iter_couple(model,t,EV_tuple,ushift,haschild,nbatch=500,verbose=False):
    
    
    setup = model.setup
    
    if verbose: start = default_timer()
    
    if gpu:
        np = cp
    else:
        np = onp
        
    nogpu = not gpu
    
    agrid = setup.agrid_c if nogpu else setup.cupy.agrid_c
    sgrid = setup.sgrid_c if nogpu else setup.cupy.sgrid_c
    
    dtype = setup.dtype
    
    
    
    
    key = 'Couple and child' if haschild else 'Couple, no children'
    
    
    taxfun = setup.taxes[key] if t < setup.pars['Tret'] else lambda x : 0.0*x
    
    ls = setup.ls_levels[key]
    if nogpu:
        uu, ux, mgrid  = setup.u_precomputed[key]['u'], setup.u_precomputed[key]['x'], setup.mgrid_c
    else:
        uu, ux, mgrid  = setup.cupy.u_precomputed[key]['u'], setup.cupy.u_precomputed[key]['x'], setup.cupy.mgrid_c
    
    upart, ucouple = (setup.u_part_k, setup.u_couple_k) if haschild else (setup.u_part_nk, setup.u_couple_nk)

    nls = len(ls)
    
    
    # type conversion is here
    
    all_t = np.array(setup.exogrid.all_t[t],copy=False)
    
    
    zf  = all_t[:,0]
    zm  = all_t[:,1]
    zftrend = setup.pars['f_wage_trend'][t]
    zmtrend = setup.pars['m_wage_trend'][t]

    psi = all_t[:,2]
    beta = setup.pars['beta_t'][t]
    sigma = setup.pars['crra_power']
    R = setup.pars['R_t'][t]

    
    
    wf = np.exp(zf + zftrend)
    wm = np.exp(zm + zmtrend)
    
    
    nexo = setup.pars['nexo_t'][t]
    shp = (setup.na,nexo,setup.ntheta)
    
    if EV_tuple is None:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = np.zeros((4,) + shp + (nls,),dtype=setup.dtype)
    else:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = EV_tuple    
    
    # type conversion
    sigma,beta = (dtype(x) for x in (sigma,beta))
    
    V_couple, c_opt, s_opt, x_opt = np.empty((4,)+shp,dtype)
    i_opt, il_opt = np.empty(shp,np.int16), np.empty(shp,np.int16)
    
    sgrid = sgrid.astype(dtype,copy=False)
    theta_val = dtype(setup.thetagrid) if nogpu else setup.cupy.thetagrid
    
    # the original problem is max{umult*u(c) + beta*EV}
    # we need to rescale the problem to max{u(c) + beta*EV_resc}
    
    istart = 0
    ifinish = nbatch if nbatch < nexo else nexo
    
    # this natually splits everything onto slices
    
    
    
    vs = setup.vsgrid_c if nogpu else setup.cupy.vsgrid_c
    
    for ibatch in range(int(np.ceil(nexo/nbatch))):
        #money_i = money[:,istart:ifinish]
        assert ifinish > istart
        
        money_t = (R*agrid, wf[istart:ifinish], wm[istart:ifinish])
        EV_t = (vs,EVr_by_l[:,istart:ifinish,:,:])
        
        
        
        
        
        V_pure_i, c_opt_i, x_opt_i, s_opt_i, i_opt_i, il_opt_i, V_all_l_i = \
           v_optimize_couple(money_t,sgrid,EV_t,mgrid,uu,ux,
                                 ls,beta,ushift,taxfun=taxfun,dtype=dtype)
           
        V_ret_i = V_pure_i + psi[None,istart:ifinish,None]
        
        
        
        
        
        V_couple[:,istart:ifinish,:] = V_ret_i # this estimate of V can be improved
        c_opt[:,istart:ifinish,:] = c_opt_i 
        s_opt[:,istart:ifinish,:] = s_opt_i
        i_opt[:,istart:ifinish,:] = i_opt_i
        x_opt[:,istart:ifinish,:] = x_opt_i
        il_opt[:,istart:ifinish,:] = il_opt_i
        #V_all_l[:,istart:ifinish,:,:] = V_all_l_i # we need this for l choice so it is ok
        
        
        istart = ifinish
        ifinish = ifinish+nbatch if ifinish+nbatch < nexo else nexo
        
        if verbose: print('Batch {} done at {} sec'.format(ibatch,default_timer()-start))
    
    
    assert np.all(c_opt > 0)
    
    psi_r = psi[None,:,None].astype(setup.dtype,copy=False)
    
    # finally obtain value functions of partners
    uf, um = upart(c_opt,x_opt,il_opt,theta_val[None,None,:],ushift,psi_r)
    uc = ucouple(c_opt,x_opt,il_opt,theta_val[None,None,:],ushift,psi_r)
    
    
    EVf_all, EVm_all, EVc_all  = (vs.apply_preserve_shape(x) for x in (EV_fem_by_l, EV_mal_by_l,EVc_by_l))
    
    V_fem = uf + beta*np.take_along_axis(np.take_along_axis(EVf_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_mal = um + beta*np.take_along_axis(np.take_along_axis(EVm_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_all = uc + beta*np.take_along_axis(np.take_along_axis(EVc_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    def r(x): return x
    
    assert V_all.dtype==EVc_all.dtype==V_couple.dtype
    
    
    return r(V_all), r(V_fem), r(V_mal), r(c_opt), r(x_opt), r(s_opt), il_opt
