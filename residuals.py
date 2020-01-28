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

xdef = np.array([0.05,0.01,0.02,0.7,0.25,0.0001,0.5,0.5])


# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience
def mdl_resid(x=xdef,save_to=None,load_from=None,return_format=['distance'],
              store_path = None,verbose=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents
 
    ulost = x[0]
    mshift=x[5]
    sigma_psi = x[1] 
    sigma_psi_init = x[1]*x[2]
    pmeet = x[3]
    pls = x[6]
    util_alp = x[4]
    util_kap = x[7]
    
    
    
    # this is for the default model
    dc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=ulost,u_lost_f=ulost,eq_split=0.0)
    sc = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00)
    
    
    
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
                        pls=pls,u_shift_mar=mshift)
    
        
        mdl = Model(iterator_name=iter_name,divorce_costs=dc,
                    separation_costs=sc,**kwords)
        mdl_list = [mdl]
        
    else:       
        mdl_list = [dill.load(open(l,'rb+')) for l in load_from]
        mdl = mdl_list[0]
        
        
    if save_to is not None:
        
        if len(save_to) > 1:
            print('warning: too much stuff is save_to')
        dill.dump(mdl,open(save_to[0],'wb+'))            
            
     
    agents = Agents( mdl_list, verbose=verbose)
    
    sim = np.array([0.0,0.2,0.4])
    dat = np.array([0.1,0.1,0.1])
    W = np.eye(sim.size)
    
    
    
    if len(dat) != len(sim):
        sim = np.full_like(dat,1.0e6)
        
    if not rel_diff:    
        res_all=(dat-sim)
    else:
        res_all = 100*(dat-sim)/(dat)
    
    
    if verbose:
        print('data moments are {}'.format(dat))
        print('simulated moments are {}'.format(sim))
    
    resid_all = np.array([x if (not np.isnan(x) and not np.isinf(x)) else 1e6 for x in res_all])
    
    resid_sc = resid_all*np.sqrt(np.diag(W)) # all residuals scaled
    
    dist = np.dot(np.dot(resid_all,W),resid_all)


    print('Distance is {}'.format(dist))
    
    
    
    out_dict = {'distance':dist,'all residuals':resid_all,
                'scaled residuals':resid_sc,'models':mdl_list,'agents':agents}
    out = [out_dict[key] for key in return_format]
    
  
            
    return out
