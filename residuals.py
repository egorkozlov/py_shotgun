#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 10:58:43 2019

@author: Egor
"""


# this defines model residuals
import numpy as np
import pickle, dill
import os


lb = np.array(   [ 0.0,  1e-4,   0.5,  0.1,  -0.2, 0.0,  0.01, 0.05,  0.05, -0.2, 0.0])
ub = np.array(   [ 2.0,  0.5,  10.0,  1.0,   0.0, 1.0,   3.0,  3.0,  0.9,    0.0, 1.0])
#xdef = np.array(  [0.5,  0.05,   2.0,  0.4, -0.05, 0.8,   0.5,  0.6,  0.3 ])
xdef = np.array([ 1.49701401,0.23225228,0.86106072,0.1669372,-0.01156311,0.10068043,0.86490734,0.23337081,0.89917949,0.0,1/3])
# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience
def mdl_resid(x=xdef,save_to=None,load_from=None,return_format=['distance'],
              store_path = None,verbose=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents
 
    mshift=x[0]
    sigma_psi = x[1] 
    sigma_psi_init = x[1]*x[2]
    pmeet = x[3]
    pmeet_t = x[4]
    pls = x[5]
    util_alp = x[6]
    util_kap = x[7]
    preg_a0 = x[8]
    preg_at = x[9]
    poutsm = x[10]
    
    
    
    # this is for the default model
    dc_k  = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00,eq_split=0.0)
    dc_nk = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00,eq_split=0.0)
    
    
    
    iter_name = 'default' if not verbose else 'default-timed'
    
    
    def join_path(name,path):
        return os.path.join(path,name)
    
    
    
    if load_from is not None:
        if type(load_from) is not list:
            load_from = [load_from]
        if store_path is not None:
            load_from = [join_path(n,store_path) for n in load_from]
    
    
    if save_to is not None:
        if type(save_to) is not list:
            save_to = [save_to]
        if store_path is not None:
            save_to = [join_path(n,store_path) for n in save_to]
    
    
                
    if load_from is None:
        
        kwords = dict(sigma_psi=sigma_psi,
                        sigma_psi_init=sigma_psi_init,
                        pmeet=pmeet,util_alp=util_alp,util_kap=util_kap,
                        pls=pls,u_shift_mar=mshift,preg_a0=preg_a0,
                        pmeet_t=pmeet_t,preg_at=preg_at,poutsm=poutsm)
    
        
        mdl = Model(iterator_name=iter_name,divorce_costs_k=dc_k,
                    divorce_costs_nk=dc_nk,**kwords)
        mdl_list = [mdl]
        
    else:       
        mdl_list = [dill.load(open(l,'rb+')) for l in load_from]
        mdl = mdl_list[0]
        
        
    if save_to is not None:
        
        if len(save_to) > 1:
            print('warning: too much stuff is save_to')
        dill.dump(mdl,open(save_to[0],'wb+'))            
            
     
    agents = Agents( mdl_list, verbose=verbose)
    
    
    n_mark = agents.state_codes['Couple and child']
    n_marnk = agents.state_codes['Couple, no children']
    n_single = agents.state_codes['Female, single']
    n_singlek = agents.state_codes['Female and child']
    
    
    is_mar = (agents.state == n_mark) | (agents.state == n_marnk)
    
    ever_mar = (np.cumsum(is_mar,axis=1) > 0)
    div_now =  (ever_mar) & ((agents.state==n_single) | (agents.state==n_singlek))
    ever_kid = ( np.cumsum( (agents.state == n_mark) | (agents.state == n_singlek),axis=1) > 0)
    
    
    
    nmar_25 = 1-ever_mar[:,4].mean()
    nmar_30 = 1-ever_mar[:,9].mean()
    nmar_35 = 1-ever_mar[:,14].mean()
    
    div_25 = div_now[:,4].mean()
    div_30 = div_now[:,9].mean()
    div_35 = div_now[:,14].mean()
    div_40 = div_now[:,19].mean()
    
    
    nkid_25 = 1-ever_kid[:,4].mean()
    nkid_30 = 1-ever_kid[:,9].mean()
    nkid_35 = 1-ever_kid[:,14].mean()
    
    nkid_25_mar = 1-ever_kid[is_mar[:,4],4].mean() if np.any(is_mar[:,4]) else 0.0
    nkid_30_mar = 1-ever_kid[is_mar[:,9],9].mean() if np.any(is_mar[:,9]) else 0.0
    nkid_35_mar = 1-ever_kid[is_mar[:,14],14].mean() if np.any(is_mar[:,14]) else 0.0
    
    
    #mkids_0_mar = (agents.state[:,1:] == n_mark)[ ~is_mar[:,0:-1] & is_mar[:,1:]].mean()
    mkids_1_mar = (agents.state[:,2:] == n_mark)[ ~is_mar[:,0:-2] & is_mar[:,2:]].mean()
    mkids_2_mar = (agents.state[:,3:] == n_mark)[ ~is_mar[:,0:-3] & is_mar[:,3:]].mean()
    mkids_3_mar = (agents.state[:,4:] == n_mark)[ ~is_mar[:,0:-4] & is_mar[:,4:]].mean()
    
    
    sim = np.array([nmar_25,nmar_30,nmar_35,
                    div_25,div_30,div_35,div_40,
                    nkid_25,nkid_30,nkid_35,
                    nkid_25_mar,nkid_30_mar,nkid_35_mar,
                    mkids_1_mar,mkids_2_mar,mkids_3_mar])
    dat = np.array([0.69,0.50,0.26,
                    0.11,0.13,0.16,0.19,
                    0.68,0.5,0.25,
                    0.43,0.27,0.14,
                    0.45,0.55,0.64])
    
    
    
    
    W = np.eye(sim.size)/(sim.size)
    
    
    
    if len(dat) != len(sim):
        sim = np.full_like(dat,1.0e6)
        
    if not rel_diff:    
        res_all=(dat-sim)
    else:
        res_all = (dat-sim)/(dat)
    
    
    if verbose:
        print('data moments are {}'.format(dat))
        print('simulated moments are {}'.format(sim))
    
    resid_all = np.array([x if (not np.isnan(x) and not np.isinf(x)) else 1e6 for x in res_all])
    
    resid_sc = resid_all*np.sqrt(np.diag(W)) # all residuals scaled
    
    dist =  np.dot(np.dot(resid_all,W),resid_all)


    print('Distance is {}'.format(dist))
    
    
    
    out_dict = {'distance':dist,'all residuals':resid_all,
                'scaled residuals':resid_sc,'models':mdl_list,'agents':agents}
    out = [out_dict[key] for key in return_format]
    
    if len(out) == 1: out = out[0]
  
    return out
