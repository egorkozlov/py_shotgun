#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 10:59:28 2020

@author: egorkozlov
"""

import numpy as np
from gridvec import VecOnGrid
from numba import jit, prange

def v_mar(setup,V,t,iassets_couple,iexo_couple,*,match_type,female):
    
    # this builds matrix for all matches specified by grid positions
    # iassets_couple (na X nmatches) and iexo_couple (nmatches)
    
    iexo, izf, izm, ipsi = setup.all_indices(t,iexo_couple)
    
    vals = pick_values(setup,V,match_type=match_type)
    V_fem, V_mal = vals['V_fem'], vals['V_mal']
    
    # obtain partner's value
    assets_partner = np.clip(setup.agrid_c[iassets_couple] - setup.agrid_s[:,None],0.0,setup.agrid_s.max())
    # this can be fastened by going over repeated values of ipsi
    v_assets_partner = VecOnGrid(setup.agrid_s,assets_partner)
    i, wnext, wthis = v_assets_partner.i, v_assets_partner.wnext, v_assets_partner.wthis
    
    if female:
        V_f_no = V_fem[:,izf]
        V_m_no = V_mal[i,izm[None,:]]*wthis + V_mal[i+1,izm[None,:]]*wnext
    else:
        V_f_no = V_fem[i,izf[None,:]]*wthis + V_fem[i+1,izf[None,:]]*wnext
        V_m_no = V_mal[:,izm]
        
    V_fem_mar, V_mal_mar = vals['V_fem_mar'], vals['V_mal_mar']
    
    ia = iassets_couple
    
    V_f_yes = V_fem_mar[ia,iexo[None,:],:]
    V_m_yes = V_mal_mar[ia,iexo[None,:],:]
    
    assert V_f_yes.shape[:-1] == V_f_no.shape
    assert V_m_yes.shape[:-1] == V_m_no.shape
    
    # fill abortion decisions 
    if match_type == 'Unplanned pregnancy' and female:
        do_abortion = vals['i_abortion'][:,izf]
    else:
        do_abortion = np.zeros(V_f_no.shape,dtype=np.bool_)
    
    
    it, wnt = setup.v_thetagrid_fine.i, setup.v_thetagrid_fine.wnext
    
    v_f, v_m, agree, nbs, itheta = get_marriage_values(V_f_yes,V_m_yes,V_f_no,V_m_no,it,wnt)
    
    return {'V_fem':v_f,'V_mal':v_m,'Agree':agree,'NBS':nbs,'itheta':itheta,'Abortion':do_abortion}

@jit(nopython=True)
def get_marriage_values(vfy,vmy,vfn,vmn,ithtgrid,wnthtgrid):
    # this is the core function that does bargaining
    
    def nbs_fun(x,y):
        if x > 0.0 and y > 0.0:
            return x*y
        else:
            return 0.0
    
    na = vfy.shape[0]
    ne = vfy.shape[1]
    #nt_crude = vfy.shape[2]
    nt = ithtgrid.size
    
    v_f = np.empty((na,ne),dtype=vfy.dtype)
    v_m = np.empty((na,ne),dtype=vfy.dtype)
    
    agree = np.zeros((na,ne),dtype=np.bool_)
    nbs = np.zeros((na,ne),dtype=vfy.dtype)
    itheta = np.empty((na,ne),dtype=np.int16)
    
    
    for ia in range(na):
        for ie in range(ne):
            
            vfy_store = np.empty((nt,),dtype=vfy.dtype)
            vmy_store = np.empty((nt,),dtype=vmy.dtype)
            
            vfn_store = vfn[ia,ie]
            vmn_store = vmn[ia,ie]
            
            # interpolate
            for it in range(nt):
                it_c = ithtgrid[it]
                it_cp = it_c+1
                wn_c = wnthtgrid[it]
                wt_c = 1.0 - wn_c
                
                vfy_store[it] = vfy[ia,ie,it_c]*wt_c + vfy[ia,ie,it_cp]*wn_c
                vmy_store[it] = vmy[ia,ie,it_c]*wt_c + vmy[ia,ie,it_cp]*wn_c
                
            # optimize  
            nbs_save = 0.0
            itheta_save = -1
            yes = False
            
            for it in range(nt):
                
                nbs_cand = nbs_fun( vfy_store[it] - vfn_store, vmy_store[it] - vmn_store )
                if nbs_cand > nbs_save:
                    nbs_save = nbs_cand
                    itheta_save = it
                    yes = True
                 
            itheta[ia,ie] = itheta_save
            nbs[ia,ie] = nbs_save            
            agree[ia,ie] = yes
            
            v_f[ia,ie] = vfy_store[it] if yes else vfn_store
            v_m[ia,ie] = vmy_store[it] if yes else vmn_store
    
    
    return v_f, v_m, agree, nbs, itheta
                    
            
    
    
    
    
    


def pick_values(setup,V,*,match_type):
    
    if match_type == 'Regular':
        
        V_fem = V['Female, single']['V']
        V_mal = V['Male, single']['V']
        
        V_fem_mar = V['Couple, no children']['VF']
        V_mal_mar = V['Couple, no children']['VM']
        
        return {'V_fem':V_fem,'V_mal':V_mal,
                'V_fem_mar':V_fem_mar,'V_mal_mar':V_mal_mar}
        
    elif match_type == 'Child immediately':
        
        V_fem = V['Female, single']['V']
        V_mal = V['Male, single']['V']
        
        V_fem_mar = V['Couple and child']['VF']
        V_mal_mar = V['Couple and child']['VM']
        
        return {'V_fem':V_fem,'V_mal':V_mal,
                'V_fem_mar':V_fem_mar,'V_mal_mar':V_mal_mar}
        
    elif match_type == 'Unplanned pregnancy':
        # expected things are here
        p_access = setup.pars['p_abortion_access']
        u_costs  = setup.pars['abortion_costs']
        
        i_abortion = (V['Female, single']['V'] - u_costs) >= V['Female and child']['V']
        
        V_fem = (1-p_access)*V['Female and child']['V'] + \
                   p_access *(np.maximum( V['Female, single']['V'] - u_costs,V['Female and child']['V'] )) \
                   - setup.pars['disutil_shotgun']
                   
        V_mal = V['Male, single']['V'] \
                - setup.pars['disutil_shotgun'] # child support enforcement here 
                
        V_fem_mar = V['Couple and child']['VF']
        V_mal_mar = V['Couple and child']['VM']
        
        return {'V_fem':V_fem,'V_mal':V_mal,'i_abortion':i_abortion,
                'V_fem_mar':V_fem_mar,'V_mal_mar':V_mal_mar}
    
    elif match_type == 'Single mother':
        
        # expected things are here
        V_fem = V['Female and child']['V'] 
        V_mal = V['Male, single']['V'] 
        
        
        V_fem_mar = V['Couple and child']['VF'] - setup.pars['disutil_marry_sm_fem']
        V_mal_mar = V['Couple and child']['VM'] - setup.pars['disutil_marry_sm_mal']
        
        return {'V_fem':V_fem,'V_mal':V_mal,
                'V_fem_mar':V_fem_mar,'V_mal_mar':V_mal_mar}
    
    else:
        raise KeyError('unknown match type')
    
    
        