#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 14:26:27 2020

@author: egorkozlov
"""

import jax.numpy as np
from jax import grad, jit
#from jax.config import config
#config.update("jax_enable_x64", True)


agrid = np.linspace(0.0,20.0,100)
au = agrid[1:]
ad = agrid[:-1]


da = (au-ad)



t0 = 10
t1 = 15

tend = t0+t1

w1 = 0.2
w0 = 1.0
sig = 1.5


R = 1.0

vlast = (w1 + R*agrid)**(1-sig)/(1-sig)

dvlast = np.fmax((vlast[1:] - vlast[:-1])/da,1e-16)

ts = 0.2

@jit
def opt_sav(money,v,au,ad,bet,R,sig):
    
    def ufun(x):
        return x**(1-sig)/(1-sig)

    dv = (v[1:] - v[:-1])/(au-ad)
    co = (bet*dv)**(-1/sig) 
    so = money[:,None] - co[None,:]
    s_best_interval = np.fmax(np.fmin(au,so),ad)
    #print(s_best_interval)
    v_best_interval = v[:-1][None,:] + dv[None,:]*(s_best_interval - ad)
    mleft = np.fmax(money[:,None] - s_best_interval,1e-10)
    Umat = ufun(mleft) + bet*v_best_interval
    vbest = Umat.max(axis=1)
    U_excess_norm = (Umat - vbest[:,None])/ts
    enorm = np.exp(U_excess_norm)
    sumexp = enorm.sum(axis=1)
    ccp = enorm/(sumexp[:,None])
    vopt = vbest + ts*np.log(sumexp)
    return ccp, so, vopt 


ag = agrid
def s0_wrap(w):
    mon = R*ag + w
    vin = vlast
    for t in range(5):
        ccp, so, vout = opt_sav(mon,vin,au,ad,1.0,1.0,1.5)
        vin = vout
    return np.sum(ccp[0]*so[0])
for pt in [0.5,1.0,1.5,2.5,3.5,10.0]:
    print(s0_wrap(pt))
    print(grad(s0_wrap)(pt))

'''

def opt_exact(money,w1,bet,R,sig):
    mult = (bet*R)**(-1/sig)
    s = np.maximum( (money - mult*w1)/(1+R*mult), 0.0)
    return s

def s0_exact(w):
    return (opt_exact(w + 0.0,w1,1.0,1.0,1.5))


print(s0_exact(0.5))
print(grad(s0_exact)(0.5))
'''





def interp_np(grid,xnew,return_wnext=True,trim=False):    
    # this finds grid positions and weights for performing linear interpolation
    # this implementation uses numpy
    
    if trim: xnew = np.minimum(grid[-1], np.maximum(grid[0],xnew) )
    
    j = np.minimum( np.searchsorted(grid,xnew,side='left')-1, grid.size-2 )
    wnext = (xnew - grid[j])/(grid[j+1] - grid[j])
    
    return j, (wnext if return_wnext else 1-wnext) 



