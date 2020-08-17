#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 15:00:07 2020

@author: egorkozlov
"""
import numpy as np
from numba import cuda, f4, f8, b1, i2
from math import ceil

try:
    import cupy as cp
except:
    pass


use_f32 = False

if use_f32:
    gpu_type = f4
    cpu_type = np.float32
else:
    gpu_type = f8
    cpu_type = np.float64


def v_mar_gpu(vfy,vmy,vfn,vmn,ithtgrid,wnthtgrid):
                    
    na, ne, nt_coarse = vfy.shape
    nt = ithtgrid.size
    
    assert nt < 500 
    
    
    
    
    v_m = cuda.device_array((na,ne),dtype=cpu_type)
    v_f = cuda.device_array((na,ne),dtype=cpu_type)
    agree = cuda.device_array((na,ne),dtype=np.bool_)
    nbs = cuda.device_array((na,ne),dtype=cpu_type)
    itheta = cuda.device_array((na,ne),dtype=np.int16)
    
    
    
    #ithtgrid = cuda.to_device(ithtgrid)
    #wnthtgrid = cuda.to_device(wnthtgrid)
    

    threadsperblock = (16, 64)
        
    b_a = ceil(na/threadsperblock[0])
    b_exo = ceil(ne/threadsperblock[1])
    
    blockspergrid = (b_a, b_exo)
    
    '''
    vfy_, vmy_ = [cuda.to_device(
                            np.ascontiguousarray(x)
                                        ) for x in (vfy, vmy)]
    
    vfn_, vmn_ = [cuda.to_device(
                                    np.ascontiguousarray(x)
                                  ) for x in (vfn,vmn)]
    '''
    vfy_, vmy_, vfn_, vmn_ = (vfy, vmy, vfn, vmn)
    
    #itht, wntht = (cuda.const.array_like(x) for x in (itht, wntht))
    
    
    cuda_ker[blockspergrid, threadsperblock](vfy_,vmy_,vfn_,vmn_,
                                                ithtgrid,wnthtgrid,
                                                    v_f,v_m,agree,nbs,itheta)
    
    '''
    v_f, v_m, agree, nbs, itheta = (x.copy_to_host() 
                            for x in (v_f, v_m, agree, nbs, itheta))
    '''
    
    v_f, v_m, agree, nbs, itheta = (cp.asarray(x) #.copy_to_host() 
                            for x in (v_f, v_m, agree, nbs, itheta))
    
    return v_f, v_m, agree, nbs, itheta

@cuda.jit
def cuda_ker(vfy,vmy,vfn,vmn,ithtgrid,wnthtgrid,v_f,v_m,agree,nbs,itheta):
    # this is the core function that does bargaining
    
    def nbs_fun(x,y):
        if x > 0.0 and y > 0.0:
            return x*y
        else:
            return 0.0
    
    na = vfy.shape[0]
    ne = vfy.shape[1]
    #nt_crude = vfy.shape[2]
    nt = ithtgrid.size
    
    ia, ie  = cuda.grid(2)   
    f1 = gpu_type(1.0)
    
    if ia < na and ie < ne:    
        
        vfy_store = cuda.local.array((500,),gpu_type)
        vmy_store = cuda.local.array((500,),gpu_type)
        
        vfn_store = vfn[ia,ie]
        vmn_store = vmn[ia,ie]
        
        # interpolate
        for it in range(nt):
            it_c = ithtgrid[it]
            it_cp = it_c+1
            wn_c = wnthtgrid[it]
            wt_c = f1 - wn_c
            
            vfy_store[it] = vfy[ia,ie,it_c]*wt_c + vfy[ia,ie,it_cp]*wn_c
            vmy_store[it] = vmy[ia,ie,it_c]*wt_c + vmy[ia,ie,it_cp]*wn_c
            
        # optimize  
        nbs_save = 0.0
        itheta_save = -1
        yes = False
        
        for it in range(nt):
            a1 = vfy_store[it] - vfn_store
            a2 = vmy_store[it] - vmn_store
            
            
            if a1 > 0 and a2 > 0:
                nbs_cand = a1*a2
            else:
                nbs_cand = 0.0
                
            if nbs_cand > nbs_save:
                nbs_save = nbs_cand
                itheta_save = it
                yes = True
                
             
        itheta[ia,ie] = itheta_save
        nbs[ia,ie] = nbs_save            
        agree[ia,ie] = yes
        
        v_f[ia,ie] = vfy_store[it] if yes else vfn_store
        v_m[ia,ie] = vmy_store[it] if yes else vmn_store
    
    
    #return v_f, v_m, agree, nbs, itheta



'''

@cuda.jit
def get_marriage_values_core(vfy,vmy,vfn,vmn,ithtgrid,wnthtgrid):
    # this is the core function that does bargaining
    
    def nbs_fun(x,y):
        if x > 0.0 and y > 0.0:
            return x*y
        else:
            return 0.0
    
    na = vfy.shape[0]
    ne = vfy.shape[1]
    #nt_crude = vfy.shape[2]
    nt = ithtgrid.size
    
    v_f = np.empty((na,ne),dtype=vfy.dtype)
    v_m = np.empty((na,ne),dtype=vfy.dtype)
    
    agree = np.zeros((na,ne),dtype=np.bool_)
    nbs = np.zeros((na,ne),dtype=vfy.dtype)
    itheta = np.empty((na,ne),dtype=np.int16)
    
    
    for ia in range(na):
        for ie in range(ne):
            
            vfy_store = np.empty((nt,),dtype=vfy.dtype)
            vmy_store = np.empty((nt,),dtype=vmy.dtype)
            
            vfn_store = vfn[ia,ie]
            vmn_store = vmn[ia,ie]
            
            # interpolate
            for it in range(nt):
                it_c = ithtgrid[it]
                it_cp = it_c+1
                wn_c = wnthtgrid[it]
                wt_c = 1.0 - wn_c
                
                vfy_store[it] = vfy[ia,ie,it_c]*wt_c + vfy[ia,ie,it_cp]*wn_c
                vmy_store[it] = vmy[ia,ie,it_c]*wt_c + vmy[ia,ie,it_cp]*wn_c
                
            # optimize  
            nbs_save = 0.0
            itheta_save = -1
            yes = False
            
            for it in range(nt):
                
                nbs_cand = nbs_fun( vfy_store[it] - vfn_store, vmy_store[it] - vmn_store )
                if nbs_cand > nbs_save:
                    nbs_save = nbs_cand
                    itheta_save = it
                    yes = True
                 
            itheta[ia,ie] = itheta_save
            nbs[ia,ie] = nbs_save            
            agree[ia,ie] = yes
            
            v_f[ia,ie] = vfy_store[it] if yes else vfn_store
            v_m[ia,ie] = vmy_store[it] if yes else vmn_store
    
    
    return v_f, v_m, agree, nbs, itheta
'''