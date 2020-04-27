#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:32:21 2020

@author: egorkozlov
"""



import matplotlib.pyplot as plt  
import numpy as np


def make_fit_plots(setup,moments,targets):
    plot_estimates(setup)
    plot_hazards(moments,targets,setup)
    


def plot_estimates(setup):
    tval = np.arange(21,35) 
    pmeet = np.zeros_like(tval,dtype=np.float64)
    ppreg = np.zeros_like(tval,dtype=np.float64)
    
    for i, t in enumerate(tval):
        pmeet[i] = setup.pars['pmeet_t'][i]
        ppreg[i] = setup.upp_precomputed[i][0]
        
    plt.figure()
    #plt.plot(tval,pmeet,label='meeting a partner')
    plt.plot(tval,ppreg*pmeet*100,label='meet and become pregnant')
    plt.legend()
    plt.title('estimated probabilities')
    plt.xlabel('age')
    plt.ylabel('probability (%)')
        
    




def plot_hazards(moments,targets,setup,ci=False):
    # graph 1: hazard of any marriage
    
    
    tval = np.arange(22,36)  
    
    for name in ['hazard of marriage','hazard of marriage & having a child','hazard of new child']:
        haz_model = np.zeros_like(tval,dtype=np.float64)
        haz_data = np.zeros_like(tval,dtype=np.float64)
        haz_data_lci = np.zeros_like(tval,dtype=np.float64)
        haz_data_uci = np.zeros_like(tval,dtype=np.float64)
        
        aux = np.zeros_like(tval,dtype=np.float64)
        name_aux = None
        
        
        for i,t in enumerate(tval):
            haz_model[i] = moments['{} at {}'.format(name,t)]
            haz_data[i] = targets['{} at {}'.format(name,t)][0]
            haz_data_lci[i] = haz_data[i] - 1.96*targets['{} at {}'.format(name,t)][1]
            haz_data_uci[i] = haz_data[i] + 1.96*targets['{} at {}'.format(name,t)][1]
            
            if name == 'hazard of marriage':
                aux[i] = setup.pars['pmeet_t'][i+1]
                name_aux = 'meeting probability'
            elif name == 'hazard of marriage & having a child':
                aux[i] = setup.pars['pmeet_t'][i+1]*setup.upp_precomputed[i+1][0]
                name_aux = 'meeting + pregnancy probability'
           
        
        fig, ax = plt.subplots()
        plt.plot(tval,haz_model*100,label='Model')
        plt.plot(tval,haz_data*100,label='Data')
        if ci: plt.plot(tval,haz_data_lci*100,label='lower 95%')
        if ci: plt.plot(tval,haz_data_uci*100,label='upper 95%')
        #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
        plt.legend()
        plt.title(name) 
        plt.xlabel('age')
        plt.ylabel('hazard (%)')
        if name != 'hazard of marriage & having a child':
            ax.set_aspect(0.5)
        else:
            ax.set_aspect(1.2)