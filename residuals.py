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
from tiktak import filer


# return format is any combination of 'distance', 'all_residuals' and 'models'
# we can add more things too for convenience

'''
w = {'divorced if k then m and one marriage':1.0,
     'divorced if m then k and one marriage':1.0}
'''
w = dict()

def mdl_resid(x=None,targets=None,weights=w,
              save_to=None,load_from=None,return_format=['distance'],
              store_path = None,verbose=False,draw=False,graphs=False,
              rel_diff=False,cs_moments=False,moments_repeat=5,
              moments_save_name=None):
    
    
    
    
    
    from model import Model
    from setup import DivorceCosts
    from nopar_probs import AgentsEst as Agents
    from crosssection import CrossSection
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
    #dc_k  = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00,eq_split=1.0)
    #dc_nk = DivorceCosts(unilateral_divorce=True,assets_kept = 1.0,u_lost_m=0.00,u_lost_f=0.00,eq_split=1.0)
    
    
    
    
    
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
        
        mdl = Model(verbose=verbose,**kwords)
        mdl_list = [mdl]
        
    else:       
        mdl_list = [dill.load(open(l,'rb+')) for l in load_from]
        mdl = mdl_list[0]
        
        
    if save_to is not None:
        
        if len(save_to) > 1:
            print('warning: too much stuff is save_to')
        dill.dump(mdl,open(save_to[0],'wb+'))            
            
    np.random.seed(18)
    agents = Agents( mdl_list, verbose=verbose, fix_seed=False)
    
    if not cs_moments:
        moments_list = [agents.compute_moments()] + [Agents( mdl_list, verbose=False, fix_seed=False).compute_moments() for _ in range(moments_repeat-1)]
    else:
        moments_list = [CrossSection(mdl_list, verbose=False, N_total=30000, fix_seed=False).compute_moments() for _ in range(moments_repeat)]
    
    
    mom = {key : np.mean([m[key] for m in moments_list],axis=0) for key in moments_list[0].keys()}
    
    
    #mom_join = Agents( mdl_list, N=10000, T=18, female=False, verbose=False).aux_moments()
    #mom_men = agents_extra.compute_moments()
    #mom.update(mom_join)
    
    
    
    if moments_save_name: # is not None can be ommited
        filer('{}.pkl'.format(moments_save_name),mom,True)
        
        
    
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
    
    

    tt = mdl_list[0].get_total_time()
    print('Distance {}, time {}'.format(dist,tt))
    
    
    
    out_dict = {'distance':dist,'all residuals':resid_all,
                'scaled residuals':resid_sc,'models':mdl_list,'agents':agents,
                'moments':mom}
    
    
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
    
    
    resid = np.zeros((num,),dtype=np.float64)
    W = np.eye(num,dtype=np.float64)
    for i, (name, targ) in enumerate(targets.items()):
        try:
            mom = moments[name]
            mom = mom if not np.isnan(mom) else 1000.0
        except KeyError:
            print('Warning: cannot find {} in moments'.format(name))
            raise(Exception('cannot compute moments'))
            
        if type(targ) is tuple:
            targ0, weights[name] = targ[0], 1/targ[1]**2 # should be inverse variances on the diagonal
            relative_here = False
        else:
            targ0 = targ
            relative_here = relative
        
        resid[i] = (mom - targ0)/targ0 if relative_here else (mom-targ0)
        
        
        try:
            W[i,i] = weights[name]
        except:
            W[i,i] = 1/num
        
        
        if report:
            print('{} is {:02.3g} (target {:02.3g})'.format(name,mom,targ0))
        
        
            
    # normalize W
    #W = W/np.sum(W)
    
    W = W/1e4
    
            
    resid_scaled = resid*np.sqrt(np.diag(W))
    dist = np.dot(np.dot(resid,W),resid)    
    
    return resid, resid_scaled, dist



        
    
