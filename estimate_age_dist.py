#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:43:21 2020

@author: egorkozlov
"""

import numpy as np
import pandas as pd
import scipy.stats as sps
import scipy.optimize as spo
from tiktak import filer
from mc_tools import int_prob


def get_new_estimates(fname='income_assets_distribution.csv',
                      save='za_dist.pkl',
                      age_start=23,age_stop=45,weighted=True):
    print('obtaining estimates from file {}'.format(fname))
    df = pd.read_csv(fname)
    
    
    income_for_normalization = np.median(np.exp(df[df['tage']==27]['log_income']))
    
    ages_array = np.arange(age_start,age_stop+1)
    estimates = np.zeros((ages_array.size,6))
    
    for i, age in enumerate(ages_array):
        
        df_age = df[df['tage']==age][['log_income','log_assets','wpfinwgt']]
        
        
        if weighted:
            weights = df_age['wpfinwgt'] / df_age['wpfinwgt'].sum()
        else:
            weights = np.ones_like(df_age['wpfinwgt'])
            weights = weights / weights.sum()
        
        mean_income = (df_age['log_income']*weights).sum()
        df_age['n_log_income'] = df_age['log_income'] - mean_income
        data = df_age[['n_log_income','log_assets']].copy()
        has_assets = data['log_assets']>0
        data['log_assets'] = sps.mstats.winsorize(data['log_assets'],(None,0.1)) - np.log(income_for_normalization)
        data['n_log_income'] = sps.mstats.winsorize(data['n_log_income'],(0.05,0.05))
        
        data_npar = pd.DataFrame({'z':data['n_log_income'].copy(),'a':data['log_assets'].copy()})
        
        
        def log_l(theta):
            mu_a = theta[0]
            sig_a = theta[1]
            sig_i = theta[2]
            rho = theta[3]
            beta_cons = theta[4]
            beta_i = theta[5]
            
            saa = sig_a**2
            sii = sig_i**2
            sia = sig_a*sig_i*rho
            cov = [[sii,sia],[sia,saa]]
            mu = [0,mu_a]
            
            li_pos = sps.multivariate_normal.logpdf(data,mean=mu,cov=cov).mean()
            li_zero = sps.norm.logpdf(data['n_log_income'],0.0,sig_i)
            p_0 = sps.norm.cdf(beta_cons + beta_i*data['n_log_income'])
            p_pos = 1-p_0
            is_positive = has_assets
            li = weights*( (is_positive)*(np.log(p_pos) + li_pos) + (~is_positive)*(np.log(p_0) + li_zero) )
            return li.sum()
        
        
        bnds_z = (None,None) if has_assets.mean() > 0 else (0.0,0.0)
        
        
        print('\n mean log assets {}, std {}'.format(data['log_assets'].mean(),np.std(data['log_assets']) ))
        print('max assets {}'.format(np.exp(data['log_assets']).max()) )
        
        
        res = spo.minimize(lambda x : -log_l(x),[3.0,0.4,0.4,0.0,0.0,0.0],
                           method = 'L-BFGS-B',
                           bounds=[(None,None),(1e-4,None),(1e-4,None),(-0.95,0.95),(None,None),bnds_z])
        
        #print('for age = {} result is {}, percentage of non-0 is {}'.format(age,res.x,has_assets.mean()))
        estimates[i,:] = res.x
        print('estimated mean log assets {}, estimated std {}\n'.format(res.x[0],res.x[1]) )
        
        
    out = {'Ages':ages_array,'mu_a':estimates[:,0],
           'sig_a':estimates[:,1],
           'sig_i':estimates[:,2],
           'rho':estimates[:,3],
           'beta0':estimates[:,4],
           'beta1':estimates[:,5]}
    
    if save: filer(save,out,True)
        
    return out


def conditional_normal_parameters(zval,mu_a,mu_z,sig_a,sig_z,rho):
    mean = mu_a + (sig_a/sig_z)*rho*(zval-mu_z)
    std = (1-rho**2)**(0.5) * sig_a  
    return mean, std


        
        
    


if __name__ == '__main__':
    dist = get_new_estimates()
    #print(dist)
    dist = filer('za_dist.pkl',0,0)
    #print(dist)
    tval = np.array([-2,-1,0,1,2],dtype=np.float64)
    x, w = np.polynomial.hermite.hermgauss(5)
    x = x*np.sqrt(2)
    w = w/w.sum()
    print((w*(x**2)).sum())
    print((w*(x**3)).sum())
    print((w*(x**4)).sum())
    
    for t in range(dist['Ages'].size):
        zval = dist['sig_i'][t]*tval
        inputs = (zval,dist['mu_a'][t],0.0,dist['sig_a'][t],dist['sig_i'][t],dist['rho'][t])
        print('inputs are {}'.format(inputs))
        mu, sig = conditional_normal_parameters(*inputs)
        share_0 = sps.norm.cdf(dist['beta0'][t] + dist['beta1'][t]*zval)
        
        # predicted assets
        a_nodes = np.exp(mu + x*sig)
        
        print((a_nodes))
        

