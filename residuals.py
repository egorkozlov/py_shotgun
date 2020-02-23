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

w = {'divorced if k then m and one marriage':10,
     'divorced if m then k and one marriage':10}

def mdl_resid(x=None,targets=None,weights=w,
              save_to=None,load_from=None,return_format=['distance'],
              store_path = None,verbose=False,draw=False,graphs=False,
              rel_diff=True):
    
    
    
    
    
    from model import Model
    from setup import DivorceCosts
    from simulations import Agents
    from calibration_params import calibration_params


    if type(x) is dict:
        kwords = x
        
        if 'targets' in x:
            targets = x.pop('targets')
        
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
    mom = agents.compute_moments()
    
    if targets is None:
        from targets import target_values
        tar = target_values()
    elif type(targets) is str:
        from targets import target_values
        tar = target_values(mode=targets)
    else:
        tar = targets
    
    resid_all, resid_sc, dist = distance_to_targets(mom,tar,weights=weights,
                                                    report=verbose)
    
    #if verbose:
    #    print('data moments are {}'.format(dat))
    #    print('simulated moments are {}'.format(sim))
    
    


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


def distance_to_targets(moments,targets,weights={},relative=True,report=False):
    num = len(targets)
    
    
    resid = np.zeros((num,),dtype=np.float32)
    W = np.eye(num,dtype=np.float32)
    for i, (name, targ) in enumerate(targets.items()):
        try:
            mom = moments[name]
        except KeyError:
            print('Warning: cannot find {} in moments'.format(name))
            raise(Exception('cannot compute moments'))
        resid[i] = (mom - targ)/targ if relative else (mom-targ)
        
        if report:
            print('{} is {:01.2g} (target {:01.2g})'.format(name,mom,targ))
        
        try:
            W[i,i] = weights[name]
        except:
            W[i,i] = 1/num
            
    resid_scaled = resid*np.sqrt(np.diag(W))
    dist = np.dot(np.dot(resid,W),resid)    
    
    return resid, resid_scaled, dist