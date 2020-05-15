#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 19:56:40 2020

@author: egorkozlov
"""

# this is to be used with renegotiation_unilateral

from numba import cuda, f4, f8, i2, b1
import numpy as np

use_f32 = False

if use_f32:
    gpu_type = f4
    cpu_type = np.float32
else:
    gpu_type = f8
    cpu_type = np.float64

from math import ceil

def v_ren_gpu_oneopt(v_y_ni, vf_y_ni, vm_y_ni, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, sig):
    
                    
    na, ne, nt_coarse = v_y_ni.shape
    nt = thtgrid.size
    
    assert nt < 500 
    
    
    
    
    v_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    vm_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    vf_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    itheta_out = cuda.device_array((na,ne,nt),dtype=np.int16)
    
    
    
    
    thtgrid = cuda.to_device(thtgrid)
    

    threadsperblock = (16, 64)
        
    b_a = ceil(na/threadsperblock[0])
    b_exo = ceil(ne/threadsperblock[1])
    
    blockspergrid = (b_a, b_exo)
    
    v_y, vf_y, vm_y = [cuda.to_device(
                            np.ascontiguousarray(x)
                                        ) for x in (v_y_ni, vf_y_ni, vm_y_ni)]
    
    vf_n, vm_n = [cuda.to_device(
                                    np.ascontiguousarray(x)
                                  ) for x in (vf_n_ni,vm_n_ni)]
    
    
    #itht, wntht = (cuda.const.array_like(x) for x in (itht, wntht))
    
    
    cuda_ker_one_opt[blockspergrid, threadsperblock](v_y, vf_y, vm_y, vf_n, vm_n, 
                                    itht, wntht, thtgrid,  
                                    v_out, vm_out, vf_out, itheta_out)
    
    v_out, vm_out, vf_out, itheta_out = (x.copy_to_host() 
                            for x in (v_out, vm_out, vf_out, itheta_out))
    
    return v_out, vf_out, vm_out, itheta_out
    



@cuda.jit   
def cuda_ker_one_opt(v_y_ni, vf_y_ni, vm_y_ni, vf_n, vm_n, itht, wntht, thtgrid, v_out, vm_out, vf_out, itheta_out):
    # this assumes block is for the same a and theta
    ia, ie  = cuda.grid(2)
    
    
    
    na = v_y_ni.shape[0]
    ne = v_y_ni.shape[1]
    nt_crude = v_y_ni.shape[2]
    nt = thtgrid.size
    
    
    f1 = gpu_type(1.0)
    
    if ia < na and ie < ne:
        
        vf_no = vf_n[ia,ie]
        vm_no = vm_n[ia,ie]
        
        v_in_store  = cuda.local.array((500,),gpu_type)
        vf_in_store = cuda.local.array((500,),gpu_type)
        vm_in_store = cuda.local.array((500,),gpu_type)
        
        
        is_good_store = cuda.local.array((500,),b1)
        it_left_store = cuda.local.array((500,),i2)
        it_right_store = cuda.local.array((500,),i2)
        it_best_store = cuda.local.array((500,),i2)
        
        
        ittc = 0
        any_good = False
        for it in range(nt):
            
            it_left_store[it] = -1
            it_right_store[it] = -1
            it_best_store[it] = -1
            
            
            it_int = itht[it]
            for ittc in range(ittc,nt_crude):
                if ittc==it_int: break
            
            
            ittp = ittc + 1
            wttp = wntht[it]
            wttc = f1 - wttp
            
            
            
            v_in_store[it]  = wttc*v_y_ni[ia,ie,ittc]  + wttp*v_y_ni[ia,ie,ittp]
            vf_in_store[it] = wttc*vf_y_ni[ia,ie,ittc] + wttp*vf_y_ni[ia,ie,ittp]
            vm_in_store[it] = wttc*vm_y_ni[ia,ie,ittc] + wttp*vm_y_ni[ia,ie,ittp]
            
            
            if vf_in_store[it] >= vf_no and vm_in_store[it] >= vm_no:
                is_good_store[it] = True
                any_good = True
            else:
                is_good_store[it] = False
            
        # everything is in local mem now
        
        # go from the right
        
        
        if not any_good:            
            for it in range(nt):
                tht = thtgrid[it]
                v_out[ia,ie,it] = tht*vf_no + (f1-tht)*vm_no
                vf_out[ia,ie,it] = vf_no
                vm_out[ia,ie,it] = vm_no
                itheta_out[ia,ie,it] = -1
            return
        
        assert any_good
        
        if is_good_store[0]: it_left_store[0] = 0        
        for it in range(1,nt):
            it_left_store[it] = it if is_good_store[it] else it_left_store[it-1]
        
        if is_good_store[nt-1]: it_right_store[nt-1] = nt-1
        for it in range(nt-2,-1,-1):
            it_right_store[it] = it if is_good_store[it] else it_right_store[it+1]
        
        # find the best number
        for it in range(nt):
            if is_good_store[it]:
                it_best_store[it] = it
            else:
                if it_right_store[it] >= 0 and it_left_store[it] >= 0:
                    dist_right = it_right_store[it] - it
                    dist_left = it - it_left_store[it]
                    assert dist_right>0
                    assert dist_left>0
                    
                    if dist_right < dist_left:
                        it_best_store[it] = it_right_store[it]
                    elif dist_right > dist_left:
                        it_best_store[it] = it_left_store[it]
                    else:                                
                        # tie breaker
                        drc = 2*it_right_store[it] - nt
                        if drc<0: drc = -drc
                        dlc = 2*it_left_store[it] - nt
                        if dlc<0: dlc = -dlc                                
                        it_best_store[it] = it_left_store[it] if \
                            dlc <= drc else it_right_store[it]
                elif it_right_store[it] >= 0:
                    it_best_store[it] = it_right_store[it]
                elif it_left_store[it] >= 0:
                    it_best_store[it] = it_left_store[it]
                else:
                    assert False, 'this should not happen'
            
            itb = it_best_store[it]
            
            
            v_out[ia,ie,it] = v_in_store[itb]
            vf_out[ia,ie,it] = vf_in_store[itb]
            vm_out[ia,ie,it] = vm_in_store[itb]
            itheta_out[ia,ie,it] = itb
            
            assert vf_out[ia,ie,it] >= vf_no
            assert vm_out[ia,ie,it] >= vm_no
            
            

def v_ren_gpu_twoopt(v_y_ni0, v_y_ni1, vf_y_ni0, vf_y_ni1, vm_y_ni0, vm_y_ni1, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, sig, 
                          rescale = True):
    
    
    na, ne, nt_coarse = v_y_ni0.shape
    nt = thtgrid.size
    assert rescale, 'no rescale is not implemented'
    
    assert nt < 500 
    
    
    
    
    v_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    vm_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    vf_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    itheta_out = cuda.device_array((na,ne,nt),dtype=np.int16)
    switch_out = cuda.device_array((na,ne,nt),dtype=np.bool_)
    pswitch_out = cuda.device_array((na,ne,nt),dtype=cpu_type)
    
    
    thtgrid = cuda.to_device(thtgrid)
    

    threadsperblock = (16, 32)
        
    b_a = ceil(na/threadsperblock[0])
    b_exo = ceil(ne/threadsperblock[1])
    
    blockspergrid = (b_a, b_exo)
    
    v_y0, vf_y0, vm_y0 = [cuda.to_device(
                            np.ascontiguousarray(x)
                                        ) for x in (v_y_ni0, vf_y_ni0, vm_y_ni0)]

    v_y1, vf_y1, vm_y1 = [cuda.to_device(
                            np.ascontiguousarray(x)
                                        ) for x in (v_y_ni1, vf_y_ni1, vm_y_ni1)]
    
    vf_n, vm_n = [cuda.to_device(
                                    np.ascontiguousarray(x)
                                  ) for x in (vf_n_ni,vm_n_ni)]
    
    
    
    
    
    
    #itht, wntht = (cuda.const.array_like(x) for x in (itht, wntht))
    
    
    cuda_ker_two_opt[blockspergrid, threadsperblock](v_y0, v_y1, vf_y0, vf_y1, vm_y0, vm_y1, vf_n, vm_n, 
                                    itht, wntht, thtgrid, sig,
                                    v_out, vm_out, vf_out, itheta_out, switch_out, pswitch_out)
    
    v_out, vm_out, vf_out, itheta_out, switch_out, pswitch_out = (x.copy_to_host() 
                            for x in (v_out, vm_out, vf_out, itheta_out, switch_out, pswitch_out))
    
    return v_out, vf_out, vm_out, itheta_out, switch_out, pswitch_out  
            
            


from math import exp, log
ofval = 14.0

@cuda.jit   
def cuda_ker_two_opt(v_y_ni0, v_y_ni1, vf_y_ni0, vf_y_ni1, vm_y_ni0, vm_y_ni1, vf_n, vm_n, itht, wntht, thtgrid, sig, v_out, vm_out, vf_out, itheta_out, switch_out, pswitch_out):
    # this assumes block is for the same a and theta
    ia, ie  = cuda.grid(2)
    
    
    
    na = v_y_ni0.shape[0]
    ne = v_y_ni0.shape[1]
    nt_crude = v_y_ni0.shape[2]
    nt = thtgrid.size
    
    
    f1 = gpu_type(1.0)
    nu = log(2.0) #0.5772156649
    correction = sig*nu
    
    if ia < na and ie < ne:
        
        vf_no = vf_n[ia,ie]
        vm_no = vm_n[ia,ie]
        
        v_in_store  = cuda.local.array((500,),gpu_type)
        vf_in_store = cuda.local.array((500,),gpu_type)
        vm_in_store = cuda.local.array((500,),gpu_type)
        
        
        is_good_store = cuda.local.array((500,),b1)
        it_left_store = cuda.local.array((500,),i2)
        it_right_store = cuda.local.array((500,),i2)
        it_best_store = cuda.local.array((500,),i2)
        
        
        ittc = 0
        any_good = False
        
        for it in range(nt):
            
            it_left_store[it] = -1
            it_right_store[it] = -1
            it_best_store[it] = -1
            
            
            it_int = itht[it]
            for ittc in range(ittc,nt_crude):
                if ittc==it_int: break
            
            
            ittp = ittc + 1
            wttp = wntht[it]
            wttc = f1 - wttp
            
            
            
            vy_0 = wttc*v_y_ni0[ia,ie,ittc] + wttp*v_y_ni0[ia,ie,ittp]
            vy_1 = wttc*v_y_ni1[ia,ie,ittc] + wttp*v_y_ni1[ia,ie,ittp]
            
            vfy_0 = wttc*vf_y_ni0[ia,ie,ittc] + wttp*vf_y_ni0[ia,ie,ittp]
            vfy_1 = wttc*vf_y_ni1[ia,ie,ittc] + wttp*vf_y_ni1[ia,ie,ittp]
            
            vmy_0 = wttc*vm_y_ni0[ia,ie,ittc] + wttp*vm_y_ni0[ia,ie,ittp]
            vmy_1 = wttc*vm_y_ni1[ia,ie,ittc] + wttp*vm_y_ni1[ia,ie,ittp]
            
            
            
            v_diff_scaled = (vy_1 - vy_0)/sig
            
            overflow_up = True if v_diff_scaled >= ofval else False 
            overflow_down = True if v_diff_scaled <= -ofval else False
            
            overflow = overflow_up or overflow_down
            
            
            if not overflow:
                emx = exp(-v_diff_scaled)
                p1 = 1.0/(1.0+emx)
                p0 = f1 - p1              
                v_smax = vy_0 + sig*log(1+exp(v_diff_scaled)) - correction                
                v_pure = v_smax - p1*vy_1 - p0*vy_0
                vf_smax = v_pure + p1*vfy_1 + p0*vfy_0
                vm_smax = v_pure + p1*vmy_1 + p0*vmy_0
            elif overflow_up:
                p1 = f1                
                v_smax = vy_1 - correction
                vf_smax = vfy_1 - correction
                vm_smax = vmy_1 - correction
            elif overflow_down:
                p1 = 0.0                
                v_smax = vy_0 - correction
                vf_smax = vfy_0 - correction
                vm_smax = vmy_0 - correction                
            else:
                assert False, 'this should not happen'
                
            
            
            pick1 = (vy_1 > vy_0) 
            pswitch_out[ia,ie,it] = p1
            switch_out[ia,ie,it] = pick1
            
            
            
            v_in_store[it]  = v_smax
            vf_in_store[it] = vf_smax
            vm_in_store[it] = vm_smax
            
            
            
            
            if vf_in_store[it] >= vf_no and vm_in_store[it] >= vm_no:
                is_good_store[it] = True
                any_good = True
            else:
                is_good_store[it] = False
            
        # everything is in local mem now
        
        # go from the right
        
        
        if not any_good:            
            for it in range(nt):
                tht = thtgrid[it]
                v_out[ia,ie,it] = tht*vf_no + (f1-tht)*vm_no
                vf_out[ia,ie,it] = vf_no
                vm_out[ia,ie,it] = vm_no
                itheta_out[ia,ie,it] = -1
            return
        
        assert any_good
        
        if is_good_store[0]: it_left_store[0] = 0        
        for it in range(1,nt):
            it_left_store[it] = it if is_good_store[it] else it_left_store[it-1]
        
        if is_good_store[nt-1]: it_right_store[nt-1] = nt-1
        for it in range(nt-2,-1,-1):
            it_right_store[it] = it if is_good_store[it] else it_right_store[it+1]
        
        # find the best number
        for it in range(nt):
            if is_good_store[it]:
                it_best_store[it] = it
            else:
                if it_right_store[it] >= 0 and it_left_store[it] >= 0:
                    dist_right = it_right_store[it] - it
                    dist_left = it - it_left_store[it]
                    assert dist_right>0
                    assert dist_left>0
                    
                    if dist_right < dist_left:
                        it_best_store[it] = it_right_store[it]
                    elif dist_right > dist_left:
                        it_best_store[it] = it_left_store[it]
                    else:                                
                        # tie breaker
                        drc = 2*it_right_store[it] - nt
                        if drc<0: drc = -drc
                        dlc = 2*it_left_store[it] - nt
                        if dlc<0: dlc = -dlc                                
                        it_best_store[it] = it_left_store[it] if \
                            dlc <= drc else it_right_store[it]
                elif it_right_store[it] >= 0:
                    it_best_store[it] = it_right_store[it]
                elif it_left_store[it] >= 0:
                    it_best_store[it] = it_left_store[it]
                else:
                    assert False, 'this should not happen'
            
            itb = it_best_store[it]
            
        
            
            
            v_out[ia,ie,it] = v_in_store[itb]
            vf_out[ia,ie,it] = vf_in_store[itb]
            vm_out[ia,ie,it] = vm_in_store[itb]
            itheta_out[ia,ie,it] = itb
            
            assert vf_out[ia,ie,it] >= vf_no
            assert vm_out[ia,ie,it] >= vm_no
               
