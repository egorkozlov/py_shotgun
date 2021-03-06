#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:25:40 2020

@author: egorkozlov
"""

import numpy as np
from tiktak import filer
import matplotlib.pyplot as plt

if __name__ == '__main__':
    V_base = filer('v_save_base.pkl',0,0)
    V_comp = filer('v_save_counterfactual.pkl',0,0)
    
    
    
    v_base = V_base[0]['Female, singlle']['V']
    v_compare = V_comp[0]['Female, singlle']['V']
    
    he = True
    v_compare(v_base,v_compare,he)
    
    

def v_compare(v_base,v_compare,he,plot=True,verbose=False):
    ia0 = 0
    t = 0
    a0 = 0.0
    
    
    from setup import ModelSetup
    d = {'high education':he}
    su_def = ModelSetup(**d)
    agrid = su_def.agrid_s
    
    zfgrid = su_def.exogrid.zf_t[t]
    zftrend = su_def.pars['f_wage_trend'][t]
    
   
    
    names = ['extra assets in base to reach 0 in compare',
             'extra assets in compare to reach 0 in base']
    
    i = 0
    
    wall = list()
    aall = list()
    
    for (v0, v1), name in zip([(v_base, v_compare),(v_compare, v_base)],
                                names):
        
        if verbose:
            print('')
            print(name)
            print('')
            print('v0 is larger in {} cases'.format(np.sum(v0>v1)))
            print('v0 is smaller in {} cases'.format(np.sum(v0<v1)))
        
        wlist = list()
        alist = list()
        
        for iz in range(v0.shape[1]):
            if verbose: print('iz is {}'.format(iz))
            if verbose: print('wage is {}'.format(np.exp(zfgrid[iz] + zftrend)))
            usd_mult = 35
            
            ia_find = np.searchsorted(v0[:,iz],v1[ia0,iz])
            
            if verbose: print('v1 is {}, v0 is {}'.format(v1[ia0,iz],v0[ia_find,iz]))
            
            wg = usd_mult*np.exp(zfgrid[iz] + zftrend)
            wlist.append(wg)
            if verbose: print('the USD wage is {}'.format(wg))
            
            a = usd_mult*(agrid[ia_find]-agrid[ia0])
            alist.append(a)
            if verbose: print('the USD equivalent level at the point is {}'.format(a))
            
            
        if plot:
            fig, ax = plt.subplots()
            plt.plot(wlist,alist,'o--k')
            plt.title('Welfare costs of unplanned pregnancy, \n single women at {}, median income = {:3.1f}'.format(21+t,wlist[3])) 
            plt.xlabel('yearly income, $1000 (2015)')
            plt.ylabel('assets equivalent (EV), $1000 (2015)')
            ax.grid(True)
        #plt.savefig('assets_eq{}.pdf'.format(i))
        
        i += 1
        
        wall.append(np.array(wlist))
        aall.append(np.array(alist))
    return wall, aall
    