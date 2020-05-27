#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:25:40 2020

@author: egorkozlov
"""

import numpy as np
from tiktak import filer

V_base = filer('v_save_base.pkl',0,0)
V_comp = filer('v_save_counterfactual.pkl',0,0)
from setup import ModelSetup
import matplotlib.pyplot as plt

d = {'high education':True}
su_def = ModelSetup(**d)
agrid = su_def.agrid_s


t = 4
a0 = 0.0
ia0 = 0#np.searchsorted(agrid,a0)
sname = 'Female, single'
oname = 'V'

zfgrid = su_def.exogrid.zf_t[t]
zftrend = su_def.pars['f_wage_trend'][t]

names = ['extra assets in base to reach 0 in compare',
         'extra assets in compare to reach 0 in base']

i = 0
for (v0, v1), name in zip([(V_base[t][sname][oname], V_comp[t][sname][oname]),(V_comp[t][sname][oname], V_base[t][sname][oname])],
                            names):
    
    print('')
    print(name)
    print('')
    print('v0 is larger in {} cases'.format(np.sum(v0>v1)))
    print('v0 is smaller in {} cases'.format(np.sum(v0<v1)))
    
    wlist = list()
    alist = list()
    
    for iz in range(v0.shape[1]):
        print('iz is {}'.format(iz))
        print('wage is {}'.format(np.exp(zfgrid[iz] + zftrend)))
        usd_mult = 35
        
        ia_find = np.searchsorted(v0[:,iz],v1[ia0,iz])
        
        print('v1 is {}, v0 is {}'.format(v1[ia0,iz],v0[ia_find,iz]))
        
        wg = usd_mult*np.exp(zfgrid[iz] + zftrend)
        wlist.append(wg)
        print('the USD wage is {}'.format(wg))
        
        a = usd_mult*(agrid[ia_find]-agrid[ia0])
        alist.append(a)
        print('the USD equivalent level at the point is {}'.format(a))
        
    fig, ax = plt.subplots()
    plt.plot(wlist,alist,'o--k')
    plt.title('Welfare costs of unplanned pregnancy, \n single women at {}, median income = {:3.1f}'.format(21+t,wlist[3])) 
    plt.xlabel('yearly income, $1000 (2015)')
    plt.ylabel('assets equivalent (EV), $1000 (2015)')
    ax.grid(True)
    plt.savefig('assets_eq{}.pdf'.format(i))
    
    i += 1
    