#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:42:52 2020

@author: egorkozlov
"""

import numpy as np
import pandas as pd
import scipy.stats as sps
import scipy.optimize as spo
from tiktak import filer

def nonpar_distribution(z,data,nbins):
    vec = z
    vm = np.concatenate( ([-np.inf],vec[:-1]) )
    vp = np.concatenate( (vec[1:],[np.inf]) )
    z_up   = 0.5*(vec + vp)
    z_down = 0.5*(vec+vm)
    nz = z.size
    
    z_probs = np.zeros((z.size,))
    a_probs_by_z = np.zeros((z.size,nbins))
    a_vals_by_z = np.zeros((z.size,nbins))
    
    for iz in range(nz):
        p_groups = np.zeros((nbins))
        a_vals = np.zeros((nbins))
        i_group = (data['z']>=z_down[iz]) & (data['z']<=z_up[iz])
        p_z = 0.0
        
        if np.any(i_group):
            
            p_z = i_group.mean()
            i_a0 = (data['a'] < 1e-2) & i_group
            p_a0 = (i_a0[i_group]*data['w'][i_group]).sum()/(data['w'][i_group].sum()) if np.any(i_a0) else 0.0
            i_ap = (~i_a0) & i_group
            a_val_pos = data['a'][i_ap]
            wa_val_pos = data['w'][i_ap]
            data['a_group'] = np.zeros_like(data['a'])
            p_groups[0] = p_a0
            
            a_unique = np.unique(a_val_pos)
            n_unique = a_unique.size
            
            if n_unique < nbins:
                for i, a in enumerate(a_unique):
                    a_vals[i+1] = a
                    p_groups[i+1] = (1-p_a0)*(wa_val_pos[a_val_pos==a].sum()/wa_val_pos.sum())
            else:
                # do histogram
                
                n_q_groups = nbins - 1 # 0 is always a group
                
                _, be = np.histogram(a_val_pos,bins=n_q_groups,weights=wa_val_pos)
                
                
                # this includes 0th quantile and 100th quantile
                for i in range(n_q_groups):
                    a_low = be[i]
                    a_high = be[i+1]
                    pick = (i_ap) & (data['a'] >= a_low) & ((data['a'] < a_high) if i<n_q_groups-1 else True)                    
                    data['a_group'][pick] = i+1
                    p_groups[i+1] = (pick[i_group]*data['w'][i_group]).sum()/(data['w'][i_group].sum())
                    a_vals[i+1] = (a_low+a_high)/2
        else: # if no one has the right z
            p_groups[0] = 1.0
            
        z_probs[iz] = p_z
        a_probs_by_z[iz,:] = p_groups
        a_vals_by_z[iz,:] = a_vals
        #assert np.all(np.diff(a_vals)>=0)
        assert np.allclose(np.sum(p_groups),1.0)
        
    return z_probs, a_probs_by_z, a_vals_by_z



def get_estimates(fname='income_assets_distribution.csv',
                      save='za_dist.pkl',
                      age_start=23,age_stop=45,weighted=True,
                      age_bw=1,
                      zlist=None):
    print('obtaining estimates from file {}'.format(fname))
    df = pd.read_csv(fname)
        
    ages_array = np.arange(age_start,age_stop+1)
    out_list = []
    for i, age in enumerate(ages_array):
        
        age_pick = (df['age']>=age-age_bw) &  (df['age']<=age+age_bw)
        data = df[age_pick][['z','a','w']]
        
        if not zlist:
            z_std = np.std(data['z'])
            zval = z_std*np.array([-2,-1,0,1,2])
        else:
            zval = zlist[i]
        
        data_npar = pd.DataFrame({'z':data['z'].copy(),
                                  'a':data['a'].copy(),
                                  'w':data['w'].copy()})
        
        out = nonpar_distribution(zval,data_npar,5)
        out_list.append(out)
        
    result = {'age':ages_array,
              'prob_z':[o[0] for o in out_list],
              'prob_a_by_z':[o[1] for o in out_list],
              'val_a_by_z':[o[2] for o in out_list]}
    #filer('za_dist.pkl',result,True)
    return result
        
        
        
        
        
if __name__ == '__main__':
    get_estimates()
        
        