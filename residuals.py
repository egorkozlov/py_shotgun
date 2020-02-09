#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 10:58:43 2019

@author: Egor
"""


# this defines model residuals
import numpy as np
import dill
import os


# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience
def mdl_resid(x=None,save_to=None,load_from=None,return_format=['distance'],
              store_path = None,verbose=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents
    from calibration_params import calibration_params


    if type(x) is dict:
        kwords = x
    else:
        lb, ub, xdef, keys, translator = calibration_params()
    
        if x is None:
            x = xdef
        kwords = translator(x)
            
        
    if verbose: print(kwords)
    
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
    is_mark = (agents.state == n_mark)
    
    ever_mar = (np.cumsum(is_mar,axis=1) > 0)
    div_now =  (ever_mar) & ((agents.state==n_single) | (agents.state==n_singlek))
    ever_kid = ( np.cumsum( (agents.state == n_mark) | (agents.state == n_singlek),axis=1) > 0)
    
    have_kid = (agents.state == n_mark) | (agents.state == n_singlek)
    num_mar = np.cumsum( agents.agreed, axis = 1 )
    one_mar = (num_mar == 1)
    
    share_x = agents.x / np.maximum(1e-3, agents.x  + agents.c)
    mean_x = share_x[:,0:20][is_mark[:,0:20]].mean()
    
    nmar_25 = 1-ever_mar[:,4].mean()
    nmar_30 = 1-ever_mar[:,9].mean()
    nmar_35 = 1-ever_mar[:,14].mean()
    nmar_40 = 1-ever_mar[:,19].mean()
    
    div_25 = div_now[ever_mar[:,4],4].mean()
    div_30 = div_now[ever_mar[:,9],9].mean()
    div_35 = div_now[ever_mar[:,14],14].mean()
    div_40 = div_now[ever_mar[:,19],19].mean()
    
    
    nkid_25 = 1-ever_kid[:,4].mean()
    nkid_30 = 1-ever_kid[:,9].mean()
    nkid_35 = 1-ever_kid[:,14].mean()
    
    nkid_25_mar = 1-ever_kid[is_mar[:,4],4].mean() if np.any(is_mar[:,4]) else 0.0
    nkid_30_mar = 1-ever_kid[is_mar[:,9],9].mean() if np.any(is_mar[:,9]) else 0.0
    nkid_35_mar = 1-ever_kid[is_mar[:,14],14].mean() if np.any(is_mar[:,14]) else 0.0
    
    
    #mkids_0_mar = (agents.state[:,1:] == n_mark)[ ~is_mar[:,0:-1] & is_mar[:,1:]].mean()
    no_kids_1_mar = 1 - ( have_kid[:,2:][ ~is_mar[:,0:-2] & is_mar[:,2:] & one_mar[:,2:]] ).mean()
    no_kids_2_mar = 1 - ( have_kid[:,3:][ ~is_mar[:,0:-3] & is_mar[:,3:] & one_mar[:,3:]] ).mean()
    no_kids_3_mar = 1 - ( have_kid[:,4:][ ~is_mar[:,0:-4] & is_mar[:,4:] & one_mar[:,4:]] ).mean()
    
    
    in_sample = (agents.k_m) | (agents.m_k)
    
    km_25 = agents.k_m[in_sample[:,4],4].mean()
    km_30 = agents.k_m[in_sample[:,9],9].mean()
    km_35 = agents.k_m[in_sample[:,14],14].mean()
    
    just_km_25 = agents.agreed_k[:,4].mean()
    just_km_30 = agents.agreed_k[:,9].mean()
    just_km_35 = agents.agreed_k[:,14].mean()
    
    share_planned = agents.planned_preg[(agents.planned_preg) | (agents.unplanned_preg)].mean()
    
    
    
    sim = np.array([nmar_25,nmar_30,nmar_35,nmar_40,
                    div_25,div_30,div_35,div_40,
                    nkid_25,nkid_30,nkid_35,
                    nkid_25_mar,nkid_30_mar,nkid_35_mar,
                    no_kids_1_mar,no_kids_2_mar,no_kids_3_mar,
                    mean_x,
                    km_25,km_30,km_35,
                    just_km_25,just_km_30,just_km_35])
    dat = np.array([0.75,0.38,0.21,0.15,
                    0.057,0.084,0.11,0.15,
                    0.90,0.60,0.34,
                    0.71,0.39,0.17,
                    0.81,0.66,0.51,
                    0.4,
                    0.21,0.12,0.10,
                    0.0058,0.0081,0.0046])
    
    
    
    
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
    
    del(out_dict)
    
    
    if 'models' not in return_format:
        for m in mdl_list:
            del(m)
        del mdl_list
        
    if 'agents' not in return_format:
        del(agents)
        
    
    if len(out) == 1: out = out[0]
  
    return out
