#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:32:21 2020

@author: egorkozlov
"""



import matplotlib.pyplot as plt  
import numpy as np


def make_fit_plots(agents,targets):
    setup = agents.setup
    moments = agents.compute_moments()
    #plot_estimates(setup)
    plot_hazards(moments,targets,setup)
    plot_cumulative(moments,targets,setup)
    plot_by_years_after_marriage(moments,targets,setup)


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
    
    
    tval = np.arange(23,36)  
    
    names = ['hazard of marriage','hazard of new child']
    captions = ["hazard of marriage:\n (% new marriages at [age]) / (% single at [age-1])",
                "hazard of new child:\n (% new births at [age]) / (% childless at [age-1])"]
    
    for name, cap in zip(names,captions):
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
        plt.title(cap) 
        plt.xlabel('age')
        plt.ylabel('hazard (%)')
        

def plot_cumulative(moments,targets,setup):
    # graph 1: hazard of any marriage
    
    import matplotlib.gridspec as gridspec
    
    tval = np.arange(23,36)  
    
    names = ['k then m in population','m then k in population','k then m in sample']
    captions = ["Kids First",
                "Marriage First",
                "Relative % of Kids First by age:\n (Kids First) / (Kids First + Marriage First)"]
    
    probs_model = []
    probs_data = []
    
    for name, cap in zip(names,captions):
        prob_model = np.zeros_like(tval,dtype=np.float64)
        prob_data  = np.zeros_like(tval,dtype=np.float64)
        
        
        
        for i,t in enumerate(tval):
            prob_model[i] = moments['{} at {}'.format(name,t)]
            prob_data[i] = targets['{} at {}'.format(name,t)][0]
        
        probs_model = probs_model + [prob_model*100]
        probs_data = probs_data + [prob_data*100]
    
    fig = plt.figure()#tight_layout=True)
    gs = gridspec.GridSpec(1,2)
    fig.suptitle('% in female population by age:')
    
    for i in range(2):
        #if i < 2:
        ax = fig.add_subplot(gs[0,i])
        #else:
        #    ax = fig.add_subplot(gs[:,1])
        ax.plot(tval,probs_model[i],label='Model')
        ax.plot(tval,probs_data[i],label='Data')
        ax.set_xlabel('age')
        if i==0: ax.set_ylabel('share (%)')
        ax.set_title(captions[i])
        ax.legend()
    
    fig, ax = plt.subplots()
    ax.plot(tval,probs_model[2],label='Model')
    ax.plot(tval,probs_data[2],label='Data')
    ax.set_xlabel('age')
    ax.set_ylabel('share (%)')
    ax.set_title(captions[2])
    ax.legend()
    
    

def plot_by_years_after_marriage(moments,targets,setup):
    # graph 1: hazard of any marriage
    
    
    yval = np.arange(1,11) 
    
    names = ['ever kids by years after marriage','divorced by years after marriage']
    captions = ['% with kids by years after marriage\n (if married)','% divorced by years after marriage \n (excluding remarried)']
    
    for name, cap in zip(names,captions):
        prob_model = np.zeros_like(yval,dtype=np.float64)
        prob_data  = np.zeros_like(yval,dtype=np.float64)
        
        for i,t in enumerate(yval):
            try:
                prob_model[i] = moments['{}, {}'.format(name,t)]
            except:
                prob_model[i] = None
                print('{}, {} not found in moments'.format(name,t))
            try:
                prob_data[i] = targets['{}, {}'.format(name,t)][0]
            except:
                prob_data[i] = None
                print('{}, {} not found in targets'.format(name,t))
            
        
        fig, ax = plt.subplots()
        
        i_data = ~np.isnan(prob_data)
        i_model = ~np.isnan(prob_model)
        plt.plot(yval[i_model],prob_model[i_model]*100,label='Model')
        plt.plot(yval[i_data],prob_data[i_data]*100,label='Data')
        #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
        plt.legend()
        plt.title(cap) 
        plt.xlabel('years after marriage')
        plt.ylabel('share (%)')
        
    