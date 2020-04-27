#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 18:53:06 2020

@author: egorkozlov
"""


import numpy as np
from numba import njit, prange, cuda
from gridvec import VecOnGrid


from platform import system
if system() != 'Darwin':
    try:
        g = cuda.device_array((2,5))
        del(g)
        ugpu = True
        upar = True
    except:
        print('no gpu mode')
        ugpu = False
        upar = False
else:
    ugpu = False
    upar = True
    

from renegotiation_unilateral_gpu import v_ren_gpu_oneopt
from renegotiation_unilateral_gpu import v_ren_gpu_twoopt

def v_ren_uni(setup,V,haschild,canswitch,t,return_extra=False,return_vdiv_only=False,rebuild=True,rescale=True):
    # this returns value functions for couple that entered the period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    # 
    # combine = True creates matrix (n_sc-by-n_inds)
    # combine = False assumed that n_sc is the same shape as n_inds and creates
    # a flat array.
     
    #Get Divorce or Separation Costs
    if haschild:
        dc = setup.divorce_costs_k
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
        V_fem = V['Female and child']['V']
    else:
        dc = setup.divorce_costs_nk
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
        V_fem = V['Female, single']['V']
        
    assert is_unil
    
    ind, izf, izm, ipsi = setup.all_indices(t+1)
    
    zfgrid = setup.exo_grids['Female, single'][t+1]
    zmgrid = setup.exo_grids['Male, single'][t+1]
    
    share=(np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) )
    relat=np.ones(share.shape)*0.5
    income_share_f=(1.0*share+0.0*relat).squeeze()
    #income_share_f =(np.exp(zfgrid[izf]+setup.pars['f_wage_trend'][t+1]) / ( np.exp(zmgrid[izm]+setup.pars['m_wage_trend'][t+1]) + np.exp(zfgrid[izf]+setup.pars['f_wage_trend'][t+1]) ) ).squeeze()
    
    share_f, share_m = dc.shares_if_split(income_share_f)
   
    
    sc = setup.agrid_c
    
    # values of divorce
    vf_n, vm_n = v_div_byshare(
        setup, dc, t, sc, share_f, share_m,
        V['Male, single']['V'], V_fem,
        izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m)
    
    assert vf_n.dtype == vm_n.dtype
    
    if return_vdiv_only:
        return {'Value of Divorce, male': vm_n,
                'Value of Divorce, female': vf_n}
    

    
    itht = setup.v_thetagrid_fine.i
    wntht = setup.v_thetagrid_fine.wnext
    thtgrid = setup.thetagrid_fine

    sig = setup.pars['taste_shock']
    
        
    if haschild:        
        
        if not ugpu:
            v_out, vf_out, vm_out, itheta_out, _, _ = \
             v_ren_core_two_opts_with_int(V['Couple and child']['V'][None,...],
                                          V['Couple and child']['VF'][None,...], 
                                          V['Couple and child']['VM'][None,...], 
                                          vf_n, vm_n,
                                          itht, wntht, thtgrid, sig)
            
        else:
            v_out, vf_out, vm_out, itheta_out = \
                            v_ren_gpu_oneopt(V['Couple and child']['V'],
                                          V['Couple and child']['VF'], 
                                          V['Couple and child']['VM'], 
                                          vf_n, vm_n,
                                          itht, wntht, thtgrid, sig)
        
            
        assert v_out.dtype == setup.dtype
         
    else:
        if canswitch:
            vcc_v, vcc_vf, vcc_vm = \
            [setup.vagrid_child_couple.apply_preserve_shape(x)
                for x in [V['Couple and child'][y] for y in ['V','VF','VM']]]
                
            if not ugpu:                
                
                    
                v_out, vf_out, vm_out, itheta_out, switch, pchoice = \
                    v_ren_core_two_opts_with_int(
                               np.stack([V['Couple, no children']['V'],vcc_v]),
                               np.stack([V['Couple, no children']['VF'],vcc_vf]), 
                               np.stack([V['Couple, no children']['VM'],vcc_vm]), 
                                        vf_n, vm_n,
                                        itht, wntht, thtgrid, sig)
                    
                
            else:    
                v_out, vf_out, vm_out, itheta_out, switch, pchoice = \
                       v_ren_gpu_twoopt(
                                       V['Couple, no children']['V'], vcc_v,
                                       V['Couple, no children']['VF'],vcc_vf, 
                                       V['Couple, no children']['VM'],vcc_vm, 
                                        vf_n, vm_n,
                                        itht, wntht, thtgrid, sig)     
        else:
            if not ugpu:
                v_out, vf_out, vm_out, itheta_out, switch, _ = \
                    v_ren_core_two_opts_with_int(
                              V['Couple, no children']['V'][None,...],
                              V['Couple, no children']['VF'][None,...], 
                              V['Couple, no children']['VM'][None,...], 
                                        vf_n, vm_n,
                                        itht, wntht, thtgrid, sig)
                    
                pchoice = np.zeros_like(itheta_out,dtype=setup.dtype)
                
            else:    
                v_out, vf_out, vm_out, itheta_out = \
                       v_ren_gpu_oneopt(
                                       V['Couple, no children']['V'],
                                       V['Couple, no children']['VF'], 
                                       V['Couple, no children']['VM'], 
                                        vf_n, vm_n,
                                        itht, wntht, thtgrid,sig)  
                pchoice = np.zeros_like(itheta_out,dtype=setup.dtype)              
                
        assert v_out.dtype == setup.dtype
        
        
    # optional rescale
    def v_rescale(v,it_out):
    
        vo = v.copy()
        itheta_in = np.broadcast_to(np.arange(thtgrid.size,dtype=np.int16)[None,None,:],it_out.shape)
        stay = (it_out!=-1)
        
        decrease = (it_out < itheta_in) & stay
        f_dec = ((thtgrid[itheta_in[decrease]])/(thtgrid[it_out[decrease]]))
        vo[decrease] = f_dec*vo[decrease]
        assert np.all(f_dec>1)
        increase = (it_out > itheta_in) & stay
        f_inc = ((1-thtgrid[itheta_in[increase]])/(1-thtgrid[it_out[increase]]))
        assert np.all(f_inc>1)
        vo[increase] = f_inc*vo[increase]
        
        return vo
    
    v_resc = v_rescale(v_out,itheta_out) if rescale else v_out.copy()
    
    if rebuild:
        v_out = thtgrid[None,None,:]*vf_out + \
            (setup.dtype(1)-thtgrid[None,None,:])*vm_out
    
    
    
    
        
    def r(x): return x
        
    
    #v_weighted = thtgrid[None,None,:]*vf_out + (1-thtgrid[None,None,:])*vm_out
    #maxd = np.max(np.abs(v_weighted - v_out))
    #meand = np.mean(np.abs(v_weighted - v_out))
    
    
    #print('t is {}, hachild is {}, canswitch is {}, max diff if {}, mean diff is {}'.\
    #      format(t,haschild,canswitch,maxd,meand))
    
    decision = np.any((itheta_out>=0),axis=2)
    result =  {'Decision': decision, 'thetas': itheta_out,
                'Values': (r(v_resc), r(v_out), r(vf_out), r(vm_out)),'Divorce':(vf_n,vm_n)}
    
    
    if not haschild:
        result['Probability of a birth'] = pchoice
        
        
       
    extra = {'Values':result['Values'],
             'Value of Divorce, male': vm_n, 'Value of Divorce, female': vf_n}
    
    
    if not return_extra:
        return result
    else:
        return result, extra




shrs_def = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.40,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]#[0.2,0.35,0.5,0.65,0.8]
def v_div_allsplits(setup,dc,t,sc,Vmale,Vfemale,izm,izf,
                        shrs=None,cost_fem=0.0,cost_mal=0.0):
    if shrs is None: shrs = shrs_def # grid on possible assets divisions    
    shp  =  (sc.size,izm.size,len(shrs))  
    Vm_divorce_M = np.zeros(shp) 
    Vf_divorce_M = np.zeros(shp)
    
    # find utilities of divorce for different divisions of assets
    for i, shr in enumerate(shrs):
        sv_m = VecOnGrid(setup.agrid_s,shr*sc - cost_mal)
        sv_f = VecOnGrid(setup.agrid_s,shr*sc - cost_fem)
        
        Vm_divorce_M[...,i] = sv_m.apply(Vmale,    axis=0,take=(1,izm),reshape_i=True) - dc.u_lost_m
        Vf_divorce_M[...,i] = sv_f.apply(Vfemale,  axis=0,take=(1,izf),reshape_i=True) - dc.u_lost_f
    
    return Vm_divorce_M, Vf_divorce_M
    


def v_div_byshare(setup,dc,t,sc,share_fem,share_mal,Vmale,Vfemale,izf,izm,
                  shrs=None,cost_fem=0.0,cost_mal=0.0):
    # this produces value of divorce for gridpoints given possibly different
    # shares of how assets are divided. 
    # Returns Vf_divorce, Vm_divorce -- values of singles in case of divorce
    # matched to the gridpionts for couples
    
    # optional cost_fem and cost_mal are monetary costs of divorce
    if shrs is None: shrs = shrs_def
    
    Vm_divorce_M, Vf_divorce_M = v_div_allsplits(setup,dc,t,sc,
                                                 Vmale,Vfemale,izm,izf,
                                shrs=shrs,cost_fem=cost_fem,cost_mal=cost_mal)
    
    # share of assets that goes to the female
    # this has many repetative values but it turns out it does not matter much
    
    
    
    
    fem_gets = VecOnGrid(np.array(shrs),share_fem)
    mal_gets = VecOnGrid(np.array(shrs),share_mal)
    
    i_fem = fem_gets.i
    wn_fem = fem_gets.wnext
    
    i_mal = mal_gets.i
    wn_mal = mal_gets.wnext
    
    inds_exo = np.arange(setup.pars['nexo_t'][t+1])
    
    
    
    Vf_divorce = (1-wn_fem[None,:])*Vf_divorce_M[:,inds_exo,i_fem] + \
                wn_fem[None,:]*Vf_divorce_M[:,inds_exo,i_fem+1]
    
    Vm_divorce = (1-wn_mal[None,:])*Vm_divorce_M[:,inds_exo,i_mal] + \
                wn_mal[None,:]*Vm_divorce_M[:,inds_exo,i_mal+1]
                
    
                
    return Vf_divorce, Vm_divorce



def v_no_ren(setup,V,haschild,canswitch,t):
    
    # this works live v_ren_new but does not actually run renegotiation
    
    expnd = lambda x : setup.v_thetagrid_fine.apply(x,axis=2)
    
    
    if haschild:
        v_y = expnd(V['Couple and child']['V'])
        vf_y = expnd(V['Couple and child']['VF'])
        vm_y = expnd(V['Couple and child']['VM'])
    else:
        v_y_nk = expnd(V['Couple, no children']['V'])
        vf_y_nk = expnd(V['Couple, no children']['VF'])
        vm_y_nk = expnd(V['Couple, no children']['VM'])
        # make a baby
        
        if canswitch:
            v_y_k = expnd(V['Couple and child']['V'])
            vf_y_k = expnd(V['Couple and child']['VF'])
            vm_y_k = expnd(V['Couple and child']['VM'])
        else:
            v_y_k = v_y_nk
            vf_y_k = vf_y_nk
            vm_y_k = vm_y_nk
            
        # switching criterion        
        switch = (v_y_k > v_y_nk)
        
        if not canswitch: assert not np.any(switch)
        
        v_y = switch*v_y_k + (~switch)*v_y_nk
        vf_y = switch*vf_y_k + (~switch)*vf_y_nk
        vm_y = switch*vm_y_k + (~switch)*vm_y_nk
     
        
    def r(x): return x
    
    shape_notheta = v_y.shape[:-1]
    yes = np.full(shape_notheta,True)
    ntheta = setup.ntheta_fine
    i_theta_out = np.broadcast_to(np.arange(ntheta,dtype=np.int16)[None,None,:],v_y.shape).copy()
        
    vf_n, vm_n = np.full((2,) + shape_notheta,-np.inf,dtype=setup.dtype)
    
        
    result =  {'Decision': yes, 'thetas': i_theta_out,
            'Values': (r(v_y), r(v_y), r(vf_y), r(vm_y)),'Divorce':(vf_n,vm_n)}
    
    if not haschild:
        result['Give a birth'] = switch
        
        
    return result


ofval = 14.0

@njit(parallel=upar)
def v_ren_core_two_opts_with_int(v_y_ni, vf_y_ni, vm_y_ni, vf_no, vm_no, itht, wntht, thtgrid, sig):
    # this takes values with no interpolation and interpolates inside
    # this also makes a choice of mar / coh
    # choice is based on comparing v_y_ni_0 vs v_y_ni_1 in the interpolated pt


    # this code is not really elegant but @njit requires some dumb things
    # note that v_y_ni has either 1 or 2 elements at 0th dimension
    # (so if two functions are passed, v_y_ni is np.stack(v_y_c,v_y_m)),
    # otherwise it is just v_y_m[None,...]. x[0] is equivalent to x[0,...].
    
    #nu = 0.5772156649
    correction = 0 #*sig*nu
    
    
    if v_y_ni.shape[0] == 2:
        nochoice = False
        v_y_ni_0, v_y_ni_1 = v_y_ni[0], v_y_ni[1]
        vf_y_ni_0, vf_y_ni_1 = vf_y_ni[0], vf_y_ni[1]
        vm_y_ni_0, vm_y_ni_1 = vm_y_ni[0], vm_y_ni[1]
    else:
        nochoice = True
        v_y_ni_0 = v_y_ni[0]
        vf_y_ni_0 = vf_y_ni[0]
        vm_y_ni_0 = vm_y_ni[0]
        
    
    dtype = v_y_ni.dtype
    
        
    na, ne, nt_coarse = v_y_ni_0.shape
    nt = thtgrid.size
    
    shp = (na,ne,nt)
    
    v_out = np.empty(shp,dtype=dtype)
    vm_out = np.empty(shp,dtype=dtype)
    vf_out = np.empty(shp,dtype=dtype)
    
    itheta_out = np.full(v_out.shape,-1,dtype=np.int16)
    ichoice_out = np.zeros(v_out.shape,dtype=np.bool_)
    pchoice_out = np.zeros(v_out.shape,dtype=dtype)
    
    
    f1 = np.float32(1)
    
    
    for ia in prange(na):
        for ie in prange(ne):
            # first we form value functions and choices
            # then we do renegotiation
            # this saves lots of operations
            
            v_opt = np.empty((nt,),dtype=dtype)
            vf_opt = np.empty((nt,),dtype=dtype)
            vm_opt = np.empty((nt,),dtype=dtype)
            
            vf_no_ae = vf_no[ia,ie]
            vm_no_ae = vm_no[ia,ie]
            
            # this part does all interpolations and maximization
            
            found_ren = False
            
            
            is_good = np.zeros((nt,),dtype=np.bool_)
            
            logit = False if sig <= 1e-3 else True
            
            
            for it in range(nt):
                it_c = itht[it]
                it_cp = it_c+1
                wn_c = wntht[it]
                wt_c = f1 - wn_c
                
                def wsum(x):
                    return x[ia,ie,it_c]*wt_c + x[ia,ie,it_cp]*wn_c
                
                def L(x):
                    if not logit:
                        return (1.0 if x>=0 else 0.0)
                    else:
                        emx = np.exp(-x)
                        return 1.0/(1.0+emx)
                
                    
                
                v_y_0 = wsum(v_y_ni_0)
                
                
                if not nochoice:
                    v_y_1 = wsum(v_y_ni_1)     
                    
                    if not logit:
                    
                        pick_1 = (v_y_1 > v_y_0)
                        
                        ichoice_out[ia,ie,it] = pick_1 
                        pchoice_out[ia,ie,it] = 1.0
                        
                        if pick_1:
                            vf_opt[it] = wsum(vf_y_ni_1)
                            vm_opt[it] = wsum(vm_y_ni_1)
                            v_opt[it] = v_y_1
                        else:
                            vf_opt[it] = wsum(vf_y_ni_0)
                            vm_opt[it] = wsum(vm_y_ni_0)
                            v_opt[it] = v_y_0
                    else:
                        
                        pick_1 = (v_y_1 > v_y_0)
                        
                        ichoice_out[ia,ie,it] = pick_1
                        
                        v_diff_scaled = (v_y_1 - v_y_0)/sig
                        
                        if (v_diff_scaled <= ofval) and (v_diff_scaled >=-ofval):
                            p1 = L(v_diff_scaled)
                            p0 = f1 - p1
                            
                            
                            v_smax = v_y_0 + sig*np.logaddexp(0.0,v_diff_scaled) - correction
                            
                            vf_0 = wsum(vf_y_ni_0)
                            vm_0 = wsum(vm_y_ni_0)
                            vf_1 = wsum(vf_y_ni_1)
                            vm_1 = wsum(vm_y_ni_1)
                            
                            #vf_smax = v_smax + p1*(vf_1-v_y_1) + p0*(vf_0-v_y_0)
                            #vm_smax = v_smax + p1*(vm_1-v_y_1) + p0*(vm_0-v_y_0)
                            v_pure = v_smax - p1*v_y_1 - p0*v_y_0
                            vf_smax = v_pure + p1*vf_1 + p0*vf_0
                            vm_smax = v_pure + p1*vm_1 + p0*vm_0
                            
                        elif (v_diff_scaled >= ofval):
                            p1 = f1
                            v_smax = v_y_1 - correction
                            vf_smax = wsum(vf_y_ni_1) - correction
                            vm_smax = wsum(vm_y_ni_1) - correction
                        elif (v_diff_scaled <= -ofval):
                            p1 = 0.0
                            v_smax = v_y_0 - correction
                            vf_smax = wsum(vf_y_ni_0) - correction
                            vm_smax = wsum(vm_y_ni_0) - correction
                        
                        
                        pchoice_out[ia,ie,it] = p1
                        vf_opt[it] = vf_smax
                        vm_opt[it] = vm_smax
                        v_opt[it] = v_smax
                        
                        
                else:                    
                    vf_opt[it] = wsum(vf_y_ni_0)
                    vm_opt[it] = wsum(vm_y_ni_0)
                    v_opt[it] = v_y_0
                
                if vf_opt[it] >= vf_no_ae and vm_opt[it] >= vm_no_ae:
                    is_good[it] = True
                    found_ren = True
                
            
            if not found_ren:
                # no search
                for it in range(nt):
                    tht = thtgrid[it]
                    v_out[ia,ie,it] = tht*vf_no_ae + (f1-tht)*vm_no_ae
                    vf_out[ia,ie,it] = vf_no_ae
                    vm_out[ia,ie,it] = vm_no_ae
                    itheta_out[ia,ie,it] = -1
            else:
                it_right = -np.ones((nt,),dtype=np.int16)
                it_left  = -np.ones((nt,),dtype=np.int16)
                if is_good[-1]: it_right[-1] = nt-1
                
                for it in range(nt-2,-1,-1):
                    it_right[it] = it if is_good[it] else it_right[it+1]
                    
                if is_good[0]: it_left[0] = 0
                for it in range(1,nt):
                    it_left[it] = it if is_good[it] else it_left[it-1]
                    
                it_best  = -np.ones((nt,),dtype=np.int16)
                
                for it in range(nt):
                    # pick the best and fill the values
                    
                    if is_good[it]:
                        it_best[it] = it
                    else:
                        if it_right[it] >= 0 and it_left[it] >= 0:
                            dist_right = it_right[it] - it
                            dist_left = it - it_left[it]
                            assert dist_right>0
                            assert dist_left>0
                            
                            if dist_right < dist_left:
                                it_best[it] = it_right[it]
                            elif dist_right > dist_left:
                                it_best[it] = it_left[it]
                            else:                                
                                # tie breaker
                                drc = 2*it_right[it] - nt
                                if drc<0: drc = -drc
                                dlc = 2*it_left[it] - nt
                                if dlc<0: dlc = -dlc                                
                                it_best[it] = it_left[it] if \
                                    dlc <= drc else it_right[it]
                        elif it_right[it] >= 0:
                            it_best[it] = it_right[it]
                        elif it_left[it] >= 0:
                            it_best[it] = it_left[it]
                        else:
                            assert False, 'this should not happen'
                    
                    itb = it_best[it]
                    v_out[ia,ie,it] = v_opt[itb]
                    vf_out[ia,ie,it] = vf_opt[itb]
                    vm_out[ia,ie,it] = vm_opt[itb]
                    itheta_out[ia,ie,it] = itb
                    
                    assert vf_out[ia,ie,it] >= vf_no_ae
                    assert vm_out[ia,ie,it] >= vm_no_ae
            
                
    return v_out, vf_out, vm_out, itheta_out, ichoice_out, pchoice_out







