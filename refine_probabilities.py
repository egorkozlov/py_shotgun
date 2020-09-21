#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27
 
@author: Egor Kozlov
"""


#import warnings
#warnings.filterwarnings("error")
 

from platform import system
     
if system() != 'Darwin' and system() != 'Windows':  
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
import gc
import dfols


 
from residuals import mdl_resid
from targets import target_values
from estimates import get_point
import numpy as np

from nopar_probs import AgentsEst    
from calibration_params import calibration_params


print('Hi!')
 
import os
os.environ['MKL_CBWR']='AUTO'

 
if __name__ == '__main__':
    
    
    high_e = True
    
    xinit, targ_mode = get_point(high_e,read_wisdom=False)
    tar = target_values(targ_mode)
    
    
    out, mdl, agents, res, mom = mdl_resid(x=xinit,targets=tar,
                                      load_from='mdl.pkl',
                                      return_format=['distance','models','agents','scaled residuals','moments'],                                      
                                      verbose=False)
    
    print('initial distance is {}'.format(out))
    
    prob_meet_init = np.array(mdl[0].setup.pars['pmeet_t'][:mdl[0].setup.pars['Tmeet']])
    prob_preg_init = np.array([mdl[0].setup.upp_precomputed_fem[t][3] for t in range(mdl[0].setup.pars['Tmeet'])])
    
    
    
    
    
    nopt = 10
    
    for iopt in range(nopt):
        
        print('running esimation round {}'.format(iopt))
        
        print('estimating probabilities:')
        
        prob_meet_est = 0.0
        prob_preg_est = 0.0
        nrep = 4 if nopt > 0 else 1
        np.random.seed(12)
        
        
        for rep in range(nrep):
            o = AgentsEst(mdl,T=30,verbose=False,fix_seed=False)
        
            prob_meet_est += (1/nrep)*o.pmeet_exo.copy()
            prob_preg_est += (1/nrep)*o.ppreg_exo.copy()
            
        print('estimated pmeet = {}'.format(prob_meet_est))
        print('estimated ppreg = {}'.format(prob_preg_est))
        
        
        
        
        
        # this does binary search
        
        w = 1.0
        
        factor = 0.5
        
        ne = prob_meet_est.size
        nw = 10
        print('reference value is {}'.format(out))
        
        y_previous = out
        
        
        
        
        yfactor = 1.5
    
        for i in range(nw):
            prob_meet_w = w*prob_meet_est + (1-w)*prob_meet_init[:ne]
            prob_preg_w = w*prob_preg_est + (1-w)*prob_preg_init[:ne]
            
            xsearch = xinit.copy()
            xsearch.update({'pmeet_exo':prob_meet_w,
                                    'ppreg_exo':prob_preg_w})
        
            out_w = mdl_resid(x=xsearch,targets=tar,
                                          return_format=['distance'],                                      
                                          verbose=False)
            
            print('with weight = {}, distance is {}'.format(w,out_w))
            
            if out_w < yfactor*out: 
                print('found a potentially imporving weight for yfactor = {}'.format(yfactor))
                break
            else:
                w = factor*w
                if i < nw-1:
                    print('trying new weight = {}'.format(w))
                else:
                    print('no luck...')
                        
                
        
        
        
        xfix = {k: xinit[k] for k in ['pmeet_21','pmeet_30','pmeet_40',
                                        'preg_21','preg_28','preg_35']}
        
        
        lb, ub, _, keys, translator = calibration_params(xfix=xfix)
        
        
        
        
        def tr(x):
            xx = translator(x)
            xx.update({'pmeet_exo':prob_meet_w,'ppreg_exo':prob_preg_w})
            return xx
        
        
        
        
        x0 = [xinit[key] for key in keys]
        
        x0, lb, ub = np.array(x0), np.array(lb), np.array(ub)
        
        print('starting from {}'.format(tr(x0)))
        
        tar = target_values('high education')
        def q(pt):
                #print('computing at point {}'.format(translator(pt)))
                try:
                    ans = mdl_resid(tr(pt),return_format=['scaled residuals'])
                except BaseException as a:
                    print('During optimization function evaluation failed at {}'.format(pt))
                    print(a)
                    ans = np.array([1e6])
                finally:
                    gc.collect()
                return ans
                
                
        res=dfols.solve(q, x0, rhobeg = 0.02, rhoend=1e-5, maxfun=60, bounds=(lb,ub),
                        scaling_within_bounds=True, objfun_has_noise=False,
                        npt = len(x0)+5,
                        user_params={'restarts.use_restarts':True,
                                     'restarts.rhoend_scale':0.5,
                                     'restarts.increase_npt':True})
        
        
        
        print(res)
        print('Result is {}'.format(tr(res.x)))
        from tiktak import filer
        filer('wisdom_refined.pkl',tr(res.x),True)
        print('wrote to the file!')
        
        
        xinit = tr(res.x)
        out, mdl, agents, res, mom = mdl_resid(x=xinit,return_format=['distance','models','agents','scaled residuals','moments'])
        if out > y_previous:
            print('no reduction in function value obtained')
            yfactor = 0.5*yfactor + 0.5
            
        y_previous = out