#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for renegotiation and marriage

Possibly inefficient but very scalable is the goal

"""


#from trans_unif import transition_uniform
import numpy as np
from gridvec import VecOnGrid



def v_mar_igrid(setup,t,V,icouple,ind_or_inds,*,female,giveabirth,
                uloss_fem,uloss_mal,uloss_fem_single=0.0,uloss_mal_single=0.0,
                unplanned_pregnancy=False,interpolate=True,return_all=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    
    
    if giveabirth:
        coup = 'Couple and child'
    else:
        coup = 'Couple, no children'
    
    
    if unplanned_pregnancy:
        assert giveabirth
        Vfem = V['Female and child']['V']
    else:
        Vfem = V['Female, single']['V']
    
    
    
    # import objects
    agrid_c = setup.agrid_c
    agrid_s = setup.agrid_s
    gamma = setup.pars['m_bargaining_weight']   
    
    
    VMval_single, VFval_single = V['Male, single']['V'] - uloss_mal_single, Vfem - uloss_fem_single
    VMval_postren, VFval_postren = V[coup]['VM'][icouple,...] - uloss_mal, V[coup]['VF'][icouple,...] - uloss_fem
    
    
    
    
    # substantial part
    ind, izf, izm, ipsi = setup.all_indices(t,ind_or_inds)
    
    
    if unplanned_pregnancy:
        pass
    
    
    # using trim = True implicitly trims things on top
    # so if sf is 0.75*amax and sm is 0.75*amax then sc is 1*amax and not 1.5
    
    #sc = sf+sm # savings of couple
    s_partner = agrid_c[icouple] - agrid_s # we assume all points on grid
    
    
    # this implicitly trims negative or too large values
    s_partner_v = VecOnGrid(agrid_s,s_partner,trim=True) 
    
    
    # this applies them
    
    if female:
        Vfs = VFval_single[:,izf]
        Vms = s_partner_v.apply(VMval_single,axis=0,take=(1,izm))
    else:
        Vms = VMval_single[:,izm]
        Vfs = s_partner_v.apply(VFval_single,axis=0,take=(1,izf))
        
        
    
    expnd = lambda x : setup.v_thetagrid_fine.apply(x,axis=2)
    
    
    Vmm, Vfm = (expnd(x[:,ind,:]) for x in 
                     (VMval_postren,VFval_postren))
    
   
    ins = [Vfm,Vmm,Vfs,Vms,gamma]
    ins = [setup.dtype(x,copy=False) for x in ins] # optional type conversion
    vfout, vmout, nbsout, agree, ithetaout = mar_mat(*ins)
    
    if not return_all:
        return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree}
    else:
        return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree, 'ins':ins}
    


def v_no_mar(setup,t,V,icouple,ind_or_inds,*,female,giveabirth):
    # emulates v_mar_igrid but with no marriage
    
    
    ind, izf, izm, ipsi = setup.all_indices(t,ind_or_inds)
    
    vmout, vfout = V['Male, single']['V'][:,izm], V['Female, single']['V'][:,izf]
    
    
    nbsout = np.zeros_like(vmout,dtype=setup.dtype)
    ithetaout = -np.ones_like(vmout,dtype=np.int16)
    agree = np.zeros_like(vmout,dtype=np.bool_)
    
    return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree}



def mar_mat(vfy,vmy,vfn,vmn,gamma):

    sf = vfy - np.expand_dims(vfn,vfn.ndim)
    sm = vmy - np.expand_dims(vmn,vmn.ndim)
    
    
    
    vfout = vfn.copy()
    vmout = vmn.copy()
    
    
    agree = (sf>0) & (sm>0)
    any_agree = np.any(agree,axis=-1)
    
    # this reshapes things
    n_agree = np.sum(any_agree)
    
    nbsout = np.zeros(vfn.shape,dtype=vfy.dtype)
    ithetaout = -1*np.ones(vfn.shape,dtype=np.int32)
    
    
    if n_agree > 0:       
        
        sf_a = sf[any_agree,:]
        sm_a = sm[any_agree,:]
        nbs_a = np.zeros(sf_a.shape,dtype=vfy.dtype)
        
        a_pos = (sf_a>0) & (sm_a>0)
        
        nbs_a[a_pos] = (sf_a[a_pos]**gamma) * (sm_a[a_pos]**(1-gamma))
        inds_best = np.argmax(nbs_a,axis=1)
        
        take = lambda x : np.take_along_axis(x,inds_best[:,None],axis=1).reshape((n_agree,))
        
        nbsout[any_agree] = take(nbs_a) 
        assert np.all(nbsout[any_agree] > 0)
        ithetaout[any_agree] = inds_best
        vfout[any_agree] = take(vfy[any_agree,:])
        vmout[any_agree] = take(vmy[any_agree,:])
        
        
        
        
    return vfout, vmout, nbsout, any_agree, ithetaout





