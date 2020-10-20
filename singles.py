#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:27:27 2020

@author: egorkozlov
"""



import numpy as onp

try:
    import cupy as cp
    gpu = True
except:
    gpu = False

nogpu = (not gpu)

from optimizers import v_optimize_couple


from marriage import v_mar

np = onp if not gpu else cp
    
        
def co(x):
    if nogpu:
        return x
    else:
        return cp.asarray(x)


def solve_singles(model,t,Vnext,ushift,female,verbose=False):
    setup = model.setup
    
    if Vnext is None: assert t > model.setup.pars['T'] - 3
        
    EV, dec = ev_single(setup,Vnext,female,t) # Vnext can be None
    model.time('Integration (s)')
    
    V_tuple = v_iter_single(setup,t,EV,female,ushift)
    model.time('Optimization (s)')
    
    return V_tuple, dec



def v_iter_single(setup,t,EV,female,ushift):
    
    agrid_s = setup.agrid_s if nogpu else setup.cupy.agrid_s
    sgrid_s = setup.sgrid_s if nogpu else setup.cupy.sgrid_s
    
    
    dtype = setup.dtype

    
    zvals = co(setup.exogrid.zf_t[t]) if female else co(setup.exogrid.zm_t[t])
    ztrend = setup.pars['f_wage_trend'][t] if female else setup.pars['m_wage_trend'][t]
    beta = setup.pars['beta_t'][t]
    R = setup.pars['R_t'][t]
    
    
    money_t = (R*agrid_s,np.exp(zvals + ztrend),np.zeros_like(zvals))
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size),dtype=setup.dtype)    
    
    sname = 'Female, single' if female else 'Male, single'
    
    taxfun = setup.taxes[sname][t] if t < setup.pars['Tret'] else lambda x : 0.0*x
    
    uu = setup.u_precomputed[sname]['u'] if nogpu else setup.cupy.u_precomputed[sname]['u']
    ux = setup.u_precomputed[sname]['x'] if nogpu else setup.cupy.u_precomputed[sname]['x']
    mgrid = setup.mgrid_s if nogpu else setup.cupy.mgrid_s
    ls = setup.ls_levels[sname]
    
    vsgrid = setup.vsgrid_s if nogpu else setup.cupy.vsgrid_s
    
    
    V, c, x, s, i_opt, ils, V_all_l = \
        v_optimize_couple(money_t,sgrid_s,(vsgrid,EV[...,None,None]),mgrid,uu,ux,
                              ls,beta,ushift,dtype=dtype,taxfun=taxfun)
    
    
    assert np.all(x<1e-5) # no x and no labor supply
    
    V_ret, c_opt, s_opt, i_opt = (q.squeeze(axis=2) for q in (V, c, s, i_opt))
    
    
    c_refine = R*agrid_s[:,None] + np.exp(zvals + ztrend)[None,:] - s_opt
     
    vs = setup.vsgrid_s if nogpu else setup.cupy.vsgrid_s
    
    EVexp = vs.apply_preserve_shape(EV)
    V_refine = setup.u(c_refine) + ushift + beta*np.take_along_axis(EVexp,i_opt,0)    
    
    def r(x): return x.astype(dtype,copy=False)
    
    return r(V_refine), r(c_opt), r(s_opt)





def ev_single(setup,V,female,t):
    if V is None: return None, {}
    
    pmeet = setup.pars['pmeet_t'][t]
    
    skip_mar = (t >= setup.pars['Tmeet'])
    
    if female:
        M = co(setup.exogrid.zf_t_mat[t].T)
        EV_nomeet =  np.dot(V['Female, single']['V'],M)
    else:
        M = co(setup.exogrid.zm_t_mat[t].T)
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
        ppreg = co(ppreg0[None,:])
        EV_meet_p, dec_p = ev_single_meet(setup,V,female,t,match_type='Unplanned pregnancy')
        EV_meet = EV_meet_np*(1-ppreg) + EV_meet_p*ppreg
    
    else:
        
        dec_p = dec_np
        EV_meet = EV_meet_np
    
    
    
    dec = {'Not pregnant':dec_np, 'Pregnant':dec_p}


        
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet, dec



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
    p_mat_ext  = co(matches['p_mat_extended'])
    icouple = co(matches['ia_c_table'])
    inds_all = co(matches['corresponding_iexo'])
    

    
    
    
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
    
    # in cupy case this is a mixture of normal and gpu arrays but it all gets converted later
    
    return EV, mout




    