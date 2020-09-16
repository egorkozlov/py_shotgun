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
from mc_tools import int_prob





# male, college: sd: exp(1.056854 -.0106599*t)
# male, college: mean: 1.066965 + .2004596*t + 1.309432*z + .0649392*zt

# female, college: sd: exp(1.153647 -.013831*t)
# female, college: mean:  .496482 + .2205861*t + .888435*z + .0874937*z*t

# male, hs: sd: exp(1.173178  -.0113634*t)
# male, hs: mean: -.5651976 + .1630382*t + 2.479825*z - .0588704*z*t

# female, hs: sd: exp(.8055752 + .0050402*t)
# female, hs: mean: .5254092 + .0645352*t + 1.717792*z + .0047126*z*t


def nonpar_distribution(z,t,data,nbins,*,female,college):
    
    if female and college:
        sd = lambda t : np.exp(1.153647 -.013831*t)
        mu = lambda z, t : .496482 + .2205861*t + .888435*z + .0874937*z*t
        
        
    if (not female) and college:
        sd = lambda t : np.exp(1.056854 -.0106599*t)
        mu = lambda z, t : 1.066965 + .2004596*t + 1.309432*z + .0649392*z*t
    
    if (female) and (not college):
        sd = lambda t : np.exp(.8055752 + .0050402*t)
        mu = lambda z, t : .5254092 + .0645352*t + 1.717792*z + .0047126*z*t
        
    if (not female) and (not college):
        sd = lambda t : np.exp(1.173178  -.0113634*t)
        mu = lambda z, t : .5254092 + .0645352*t + 1.717792*z + .0047126*z*t
        
    
    
    vec = z
    vm = np.concatenate( ([-np.inf],vec[:-1]) )
    vp = np.concatenate( (vec[1:],[np.inf]) )
    z_up   = 0.5*(vec + vp)
    z_down = 0.5*(vec+vm)
    nz = z.size
    
    z_probs = np.zeros((z.size,))
    a_probs_by_z = np.zeros((z.size,nbins))
    a_vals_by_z = np.zeros((z.size,nbins))
    
    
    if nbins == 5:
        sa = np.array([-1.5,-0.75,0.0,0.75,1.5])
    else:
        sa = sps.norm.ppf([(i+1)/(nbins+1) for i in range(nbins)])
    
    
    
    
    muz = np.average(data['z'],weights=data['w'])
    sdz = np.sqrt(np.average((data['z']-muz)**2,weights=data['w']))
    
    
    pz_all = int_prob( vec, mu=muz, sig=sdz )
    
    for iz in range(nz):
        p_groups = np.zeros((nbins))
        a_vals = np.zeros((nbins))
        
        
            
            
        
        p_z = pz_all[iz]
        
    
        
        
        vsd = sd(t)
        vmu = mu(z[iz], t)
        
        la_vals = vmu + sa*vsd
        p_groups = int_prob( la_vals,mu=vmu,sig=vsd,trim=False) 
        
        
        a_vals = np.exp(la_vals)*(la_vals >= 0)
            
            
            
        z_probs[iz] = p_z
        a_probs_by_z[iz,:] = p_groups
        a_vals_by_z[iz,:] = a_vals
        #assert np.all(np.diff(a_vals)>=0)
        assert np.allclose(np.sum(p_groups),1.0)
    
    assert np.allclose(z_probs.sum(),1.0)
    return z_probs, a_probs_by_z, a_vals_by_z

    
    



def get_estimates(fname='income_assets_distribution.csv',
                      save='za_dist.pkl',
                      weighted=True,                      
                      zlist=None,*,
                      age_start,age_stop,
                      female,college
                      ):
    print('obtaining estimates from file {}'.format(fname))
    df = pd.read_csv(fname)
    assert not np.any(np.isnan(df))
        
    ages_array = np.arange(age_start,age_stop+1)
    out_list = []
    
    
    amin = df['age'].min()
    
    for i, age in enumerate(ages_array):
        t = age - max(age_start,amin)
        
        
        age_pick = (df['age']==max(age,amin))
        data = df[age_pick][['z','a','w']]
        
        if not zlist:
            z_std = np.std(data['z'])
            zval = z_std*np.array([-2,-1,0,1,2])
        else:
            zval = zlist[i]
        
        data_npar = pd.DataFrame({'z':data['z'].copy(),
                                  'a':data['a'].copy(),
                                  'w':data['w'].copy()})
        
        out = nonpar_distribution(zval,t,data_npar,5,female=female,college=college)
        out_list.append(out)
        
    result = {'age':ages_array,
              'prob_z':[o[0] for o in out_list],
              'prob_a_by_z':[o[1] for o in out_list],
              'val_a_by_z':[o[2] for o in out_list]}
    #filer('za_dist.pkl',result,True)
    return result
        
        
        
        
        
if __name__ == '__main__':
    get_estimates()
        
        