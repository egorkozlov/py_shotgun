#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for setting up the model.

Created on Tue Oct 29 17:01:07 2019

@author: egorkozlov
"""


#from platform import system

import numpy as np
from timeit import default_timer
#from numba import njit, vectorize
#from memory_profiler import profile
#from IPython import get_ipython
#from asizeof import asizeof
import os
import psutil


#if system() != 'Darwin':
from setup import ModelSetup
from couples import solve_couples
from singles import solve_singles
from single_moms import solve_single_mom




try:
    import cupy as cp
    gpu = True
except:
    gpu = False

class Model(object):
    def __init__(self,verbose=False,**kwargs):
        self.mstart = self.get_mem()
        self.mlast = self.get_mem()
        self.verbose = verbose
        self.start = default_timer()
        self.last = default_timer()
        self.time_dict = dict()
        
        self.setup = ModelSetup(**kwargs)
        self.dtype = self.setup.dtype
        self.time('Creation of setup')
        
        self.solve()

    def get_mem(self):
        return psutil.Process(os.getpid()).memory_info().rss/1e6
        
        
    def time(self,whatisdone,time_lb=0.5,verbose=True):
        verbose = self.verbose
        total_time = default_timer() - self.start
        last_time = default_timer() - self.last

        total_mem = self.get_mem() / 1000
        last_mem  = self.get_mem() - self.mlast

        def r(x): return round(x,2)

        if last_time > time_lb:
            if verbose: print('{} took {}, total {}, mem {}G'.format(whatisdone,r(last_time),r(total_time),r(total_mem)))
        self.last = default_timer()
        self.mlast = self.get_mem()

        if whatisdone in self.time_dict:
            self.time_dict[whatisdone] = self.time_dict[whatisdone] + [last_time]
        else:
            self.time_dict[whatisdone] = [last_time]

    def get_total_time(self):
        return round(default_timer() - self.start,2)
    
    def time_statistics(self,remove_worst=True,remove_single=False):

        print('Total time is {}'.format(default_timer() - self.start))
        for what, timelist in self.time_dict.items():

            if remove_single and len(timelist) == 1: continue

            time_arr = np.array(timelist)

            extra = ''
            if remove_worst and time_arr.size > 1:
                time_worst = time_arr.max()
                time_arr = time_arr[time_arr<time_worst]
                extra = ' (excl the worst)'

            av_time = round(np.mean(time_arr),2) 
            tot_time = round(np.sum(np.array(timelist)),2) 
            print('On average {} took {}, total {} sec{}'.format(what,av_time,tot_time,extra))



    def v_next(self,desc,t,V_next):
        # This comptues next period value function.
        # V_next can be None, in this case this is an initial iteration
        setup = self.setup
        ushift = setup.utility_shifters[desc]



        if desc == 'Female, single' or desc == 'Male, single':



            female = (desc == 'Female, single')

            
            (V, c, s), dec = solve_singles(self,t,V_next,ushift,female)      
            
            if gpu:
                V, c, s = [cp.asnumpy(x) for x in (V, c, s)]
                dec = {k : {d: cp.asnumpy(dec[k][d]) for d in dec[k]} for k in dec}
                self.time('Transfer from gpu')


            return {desc: {'V':V,'c':c,'s':s,'x':0.0*c}}, {desc: dec}
        
        elif desc== 'Couple and child' or desc == 'Couple, no children':
            
            
            haschild = (desc== 'Couple and child')
            
            (V, VF, VM, c, x, s, fls), dec = solve_couples(self,t,V_next,ushift,haschild)
            
            
            if gpu:
                V, VF, VM, c, x, s, fls = [cp.asnumpy(x) for x in (V, VF, VM, c, x, s, fls)]
                dec = {d: cp.asnumpy(dec[d]) for d in dec}
                self.time('Transfer from gpu')
            
            
            return {desc: {'V':V,'VF':VF,'VM':VM,'c':c,'x':x,'s':s,'fls':fls}},\
                   {desc: dec}
        
        elif desc == 'Female and child':
            
                                    
            (V, c, x, s, fls), dec = solve_single_mom(self,t,V_next,ushift)
            
            if gpu:
                (V, c, x, s, fls) = [cp.asnumpy(x) for x in (V, c, x, s, fls)]
                dec = {k : {d: cp.asnumpy(dec[k][d]) for d in dec[k]} for k in dec}
                self.time('Transfer from gpu')
            
            return {desc: {'V':V,'c':c,'x':x,'s':s,'fls':fls}},\
                   {desc: dec}
                
        
        else:
            raise Exception('I do not know this type...')
    
            
    
    def solve(self,save=False):
        
        
        T = self.setup.pars['T']
        self.V = list()
        self.decisions = list()
        
        
        for t in reversed(range(T)):
            Vnow = dict()
            decnow = dict()
            if self.verbose: print('solving for t = {}'.format(t))
            
            Vnext = self.V[0] if t<T-1 else None
            
            if gpu and (Vnext is not None):
                Vnext = {d: {k: cp.array(Vnext[d][k]) for k in Vnext[d]} 
                                                        for d in Vnext}    
                self.time('Transfer to gpu')
            
            for desc in self.setup.state_names:
                V_d, dec = self.v_next(desc,t,Vnext)                   
                Vnow.update(V_d)
                decnow.update(dec)
                
            self.V = [Vnow] + self.V
            self.decisions = [decnow] + self.decisions
            
            
            
            #if show_mem:
            #    print('The size of V is {} giga'.format(asizeof(self.V)/1000000000))
            #    print('The size of decisions is {} giga'.format(asizeof(self.decisions)/1000000000))
        if save:
            import pickle
            pickle.dump(self,open('model_save.pkl','wb+'))
        
        
    def x_reshape(self,x,t):
        # this reshapes couple's values to multidimensional objects
        ss = self.setup
        
        shp = (ss.na,ss.pars['n_zf_t'][t],
                              ss.pars['n_zm_t'][t],
                              ss.pars['n_psi_t'][t],
                              x.shape[-1])
        
        x_reshape = x.reshape(shp)
        return x_reshape
    
        
    #def graph(self,ai,zfi,zmi,psii,ti,thi):        
    #    #Draw some graph of Value and Policy Functions
    #    V=graphs(self,ai,zfi,zmi,psii,ti,thi)        
    #    return V
    
    def get_graph_values(self,fun='V',dec=False,iassets=slice(None),iexo=slice(None),itheta=slice(None),*,state,t):
        Vin = self.V[t][state][fun] if not dec else self.decisions[t][state][fun]
        
        if fun == 'fls': Vin = self.setup.ls_levels[state][Vin]
        
        if fun == 'thetas':
            Vin_old = Vin
            Vin = self.setup.thetagrid_fine[Vin]
            Vin[Vin_old==-1] = None
        
        
        
        has_theta = (Vin.ndim == 3)
        
        if type(iexo) is tuple:
            Vin = self.x_reshape(Vin,t)
            ind = (iassets,) + iexo
        else:
            ind = (iassets,iexo)
        
        agrid = self.setup.agrid_c if has_theta else self.setup.agrid_s
        if has_theta:
            exogrids = (self.setup.exogrid.zf_t[t],self.setup.exogrid.zm_t[t],self.setup.exogrid.psi_t[t])
        else:
            if state == 'Female, single' or state == 'Female and child':
                exogrids = (self.setup.zf_t[t],)
            else:
                exogrids = (self.setup.zm_t[t],)
        thetagrid = ((self.setup.thetagrid,) if not dec else (self.setup.thetagrid_fine,)) \
                            if has_theta else tuple()        
        grids = (agrid,) + exogrids + thetagrid
        
    
        
        if has_theta: ind = ind+(itheta,)
        
        grid_out = None
        
        for num, ii in enumerate(ind):
            if ii == slice(None): grid_out = grids[num]
        
        return grid_out, Vin[ind]
    
    
    def diagnostics(self):
        #from diagnostics import check_value_functions
        #check_value_functions(self)
        
        from graphs import v_graphs
        v_graphs(self)
        
    def mar_graphs(self):
        from mar_graphs import mar_graphs
        mar_graphs(self)
        
      
        
    
