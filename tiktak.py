# -*- coding: utf-8 -*-
"""
Implementation of TikTak code as described in:
    
    'Benchmarking Global Optimizers'
     by Antoine Arnoud,Fatih Guvenen and Tatjana Kleineberg'

@author: Fabio
"""
import sobol_seq
import numpy as np
from scipy.optimize import minimize
from p_client import compute_for_values
from time import sleep
import pickle
from calibration_params import calibration_params

def tiktak(*,N,N_st,xfix=None,skip_global=False,skip_local=False,resume_global=False,resume_local=False):
    
    xl, xu, x0, keys, translator = calibration_params(xfix=xfix)
    #Initial cheks
    assert len(xl)==len(xu)
    
    assert N>=N_st
    
    x0 = np.array(x0)
    
    
    ############################################
    #1 INITIALIZATION
    ###########################################
    
    
    
    if not skip_global:
        #First Create a Sobol Sequence
        init = sobol_seq.i4_sobol_generate(len(xl),N) # generate many draws from uniform
        #init=init[:,0]   
        
        #Get point on the grid
        x_init=xl*(1-init)+xu*init
        x_init=x_init.T
        x_init=x_init.squeeze()
    
        #Get fitness of initial points
        
        pts = [ ('compute',translator(x_init[:,j])) for j in range(N)]
        fx_init = compute_for_values(pts,resume=resume_global)
        
        fx_init = (np.array(fx_init)).squeeze()
         # !! not the optimizer returns squared value of mdl_resid
        
        #Sort in ascending order of fitness
        order=np.argsort(fx_init)
        fx_init=fx_init[order]
        x_init=x_init[:,order]
        
        filer('sobol_results.pkl',(fx_init,x_init),True)
        print('saved the results succesfully')
    else:        
        (fx_init,x_init0) = filer('sobol_results.pkl',None,False)
        
        # this block appends variables if too little are specified in sobol_results
        x_init = np.zeros((x0.size,x_init0.shape[1]),dtype=x_init0.dtype)
        x_init[:,:] = x0[:,None]
        x_init[0:(x_init0.shape[0]),:] = x_init0        
        
        
        print('loaded the results from the file')
        
        
    #Take only the first N_st realization
    fx_init=fx_init[0:N_st]
    x_init=x_init[:,0:N_st]
   
    if skip_local:
        print('local minimizers are skipped')
        return x_init[:,0]
    
    
    #Create a file with sobol sequence points
    filer('sobol.pkl',x_init,True)    
    
    
    if not resume_local:
        #List containing parameters and save them in file
        param=list([ (fx_init[0], x_init[:,0])])
        filer('wisdom.pkl',param,True)
        i_start = 0
    else:
        param = filer('wisdom.pkl',None,False)
        i_start = len(param)
    
    vals = [('minimize',(i,N_st,xfix)) for i in range(i_start,N_st)]
    
    compute_for_values(vals,timeout=7200.0)
    
    param = filer('wisdom.pkl',None,write=False)
    
    ############################################
    #3 TOPPING RULE
    ###########################################
    #print(999,ite)
    #Final Refinement
    
    
    return param[0]
    
##########################################
#Functions
#########################################
    
#Write on Functionsbtach worker_run.sh
def filer(filename,array,write=True,repeat=True):
    
    while True:
        try:
            if write:
                with open(filename, 'wb+') as file:
                    pickle.dump(array,file)
            else:
                with open(filename, 'rb') as file:
                    array=pickle.load(file)
                return array
                
            break
        except KeyboardInterrupt:
            raise KeyboardInterrupt()
        except:
            print('Problems opening the file {}'.format(filename))
            if not repeat: raise Exception('could not open the file')
            #sleep(0.5)
    
