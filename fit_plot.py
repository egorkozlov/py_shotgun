#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:32:21 2020

@author: egorkozlov
"""



import matplotlib.pyplot as plt


import numpy as np


def make_fit_plots(agents,targets,moments_aux=None):
    setup = agents.setup
    moments = agents.compute_moments()
    try:
        plot_estimates(setup)
    except:
        print('failed to plot estimates')
    
    try:
        plot_hazards(moments,targets,setup)
    except:
        print('failed to plot hazards')
        
    try:
        plot_cumulative(moments,targets,setup)
    except:
        print('failed to plot cumulative')
        
    try:
        plot_by_years_after_marriage(moments,targets,setup)
    except:
        print('failed to plot by years after marriage')
        
    try:
        plot_kfmf(moments,targets,setup)
    except:
        print('failed to plot kfmf')
        
    if moments_aux is not None: 
        try:
            plot_men(moments_aux,targets)
        except:
            print('failed to plot men')
    try:
        plot_kfmf_ref(moments,targets,setup)
    except:
        print('failed to plot ref')
    


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
        plt.plot(tval,haz_model*100,'o-k',label='Model')
        plt.plot(tval,haz_data*100,'o--k',label='Data')
        if ci: plt.plot(tval,haz_data_lci*100,label='lower 95%')
        if ci: plt.plot(tval,haz_data_uci*100,label='upper 95%')
        #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
        ax.grid(True)
        xticks = [i for i in range(22,36)]
        ax.set_xticks(xticks)
        plt.legend()
        plt.title(cap) 
        plt.xlabel('age')
        plt.ylabel('hazard (%)')
        plt.savefig('{}.pdf'.format(name))
        

def plot_cumulative(moments,targets,setup):
    # graph 1: hazard of any marriage
    
    import matplotlib.gridspec as gridspec
    
    tval = np.arange(23,36)  
    
    names = ['k then m in population','m then k in population','k then m in sample','ever married']
    captions = ["Kids First",
                "Marriage First",
                "Relative % of Kids First by age:\n (Kids First) / (Kids First + Marriage First)",
                "ever married"]
    
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
    gs = gridspec.GridSpec(1,1)
    fig.suptitle('% in female population by age:')
    
    ch = ['o','x']
    nm = ['KF','MF']
    for i in range(2):
        #if i < 2:
        ax = fig.add_subplot(gs[0,0])
        #else:
        #    ax = fig.add_subplot(gs[:,1])
        ax.plot(tval,probs_model[i],'{}-k'.format(ch[i]),label='{}: model'.format(nm[i]))
        ax.plot(tval,probs_data[i],'{}--k'.format(ch[i]),label='{}: data'.format(nm[i]))
    ax.set_xlabel('age')
    ax.set_ylabel('share (%)')
    ax.set_title('Kids First and Marriage First')
    ax.legend()
    #yticks = [i*5 for i in range(12)]
    #ax.set_yticks(yticks)
    xticks = [i for i in range(22,36)]
    ax.set_xticks(xticks)
    ax.grid(True)
    plt.savefig('popshares.pdf')
    
    fig, ax = plt.subplots()
    ax.plot(tval,probs_model[2],'o-k',label='Model')
    ax.plot(tval,probs_data[2],'o--k',label='Data')
    ax.set_xlabel('age')
    ax.set_ylabel('share (%)')
    ax.set_title(captions[2])
    ax.legend()
    ax.set_xticks(xticks)
    ax.grid(True)
    plt.savefig('relshares.pdf')
    
    
    fig, ax = plt.subplots()
    ax.plot(tval,probs_model[3],'o-k',label='Model')
    ax.plot(tval,probs_data[3],'o--k',label='Data')
    ax.set_xlabel('age')
    ax.set_ylabel('share (%)')
    ax.set_title(captions[3])
    ax.legend()
    ax.set_xticks(xticks)
    ax.grid(True)
    plt.savefig('evermar.pdf')
    

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
        plt.plot(yval[i_model],prob_model[i_model]*100,'o-k',label='Model')
        plt.plot(yval[i_data],prob_data[i_data]*100,'o--k',label='Data')
        #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
        plt.legend()
        plt.title(cap) 
        plt.xlabel('years after marriage')
        plt.ylabel('share (%)')
        ax.grid(True)
        ax.set_xticks(yval)
        plt.savefig('{}.pdf'.format(name))
        
        

def plot_kfmf(moments,targets,setup):
    # graph 1: hazard of any marriage
    
    
    yval = np.arange(1,11) 
    
    names = ['divorced by years after marriage if kids first','divorced by years after marriage if marriage first']
    captions = ['KF','MF']
    mrkrs = ['x','o']
    
    
    fig, ax = plt.subplots()
    
    for name, cap, mrkr in zip(names,captions,mrkrs):
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
        
        print((cap,prob_model))
            
        i_data = ~np.isnan(prob_data)
        i_model = ~np.isnan(prob_model)
        plt.plot(yval[i_model],prob_model[i_model]*100,'{}-k'.format(mrkr),label='{}: model'.format(cap))
        plt.plot(yval[i_data],prob_data[i_data]*100,'{}--k'.format(mrkr),label='{}: data'.format(cap))
        #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
    plt.legend()
    plt.title('Divorced by years after marriage by groups') 
    plt.xlabel('years after marriage')
    plt.ylabel('share (%)')
    ax.grid(True)
    ax.set_xticks(yval)
    plt.savefig('div_kfmf.pdf')
    
    # ref_kf = [0.08114558, 0.10186335, 0.11817027, 0.11984021, 0.11823362, 0.11746032, 0.13162393, 0.15057915, 0.14123007, 0.14550265]
    # ref_mf = [0.        , 0.01083856, 0.02031424, 0.02864134, 0.03768433, 0.04414536, 0.05390435, 0.05976757, 0.0625543 , 0.06301748]
    
def plot_kfmf_ref(moments,targets,setup):
    yval = np.arange(1,11) 
    #ref_kf = np.array([0.08114558, 0.10186335, 0.11817027, 0.11984021, 0.11823362, 0.11746032, 0.13162393, 0.15057915, 0.14123007, 0.14550265])
    #ref_mf = np.array([0.        , 0.01083856, 0.02031424, 0.02864134, 0.03768433, 0.04414536, 0.05390435, 0.05976757, 0.0625543 , 0.06301748])
    
    ref_kf = np.array([0.03219107, 0.06490649, 0.08013937, 0.08991495, 0.09724047,
       0.10518732, 0.11538462, 0.12629758, 0.13806706, 0.15437788])
    
    ref_mf = np.array([0.        , 0.00940496, 0.01853282, 0.02416011, 0.02771235,
       0.03738832, 0.04405594, 0.05178366, 0.05938639, 0.06609948])
    
    here_kf = np.zeros_like(yval,dtype=np.float64)
    here_mf = np.zeros_like(yval,dtype=np.float64)
    
    names = ['divorced by years after marriage if kids first','divorced by years after marriage if marriage first']
    conts = [here_kf,here_mf]
    
    
    
    
    for name, cont in zip(names,conts):
        for i,t in enumerate(yval):
            try:
                cont[i] = moments['{}, {}'.format(name,t)]
            except:
                cont[i] = None
                print('{}, {} not found in moments'.format(name,t))
    
    
    fig, ax = plt.subplots()
    plt.plot(yval,ref_kf*100,'x-k',label='KF: baseline')
    plt.plot(yval,here_kf*100,'o-b',label='KF: $\phi_s = 0$')
    plt.plot(yval,ref_mf*100,'x--k',label='MF: baseline')
    plt.plot(yval,here_mf*100,'o--b',label=r'MF: $\phi_s = 0$')
    plt.legend()
    plt.title('Divorced by years after marriage: social stigma matters') 
    plt.xlabel('years after marriage')
    plt.ylabel('share (%)')
    ax.grid(True)
    ax.set_xticks(yval)
    plt.savefig('div_kfmf_ref.pdf')
    
    
    
def plot_men(moments,targets):
    
    tval = np.arange(24,36)  
    
    
    names = ['men, relative income just married / single','men, relative income with kids / no kids']
    captions = ['men, relative income just married / single','married men, relative income with kids / no kids']
    fnames = ['men_mar_ratio','men_kids_ratio']
    
    for name, cap, fname in zip(names,captions,fnames):
        haz_model = np.zeros_like(tval,dtype=np.float64)
        haz_data = np.zeros_like(tval,dtype=np.float64)
        haz_data_lci = np.zeros_like(tval,dtype=np.float64)
        haz_data_uci = np.zeros_like(tval,dtype=np.float64)
        
        
        
        for i,t in enumerate(tval):
            haz_model[i] = moments['{} at {}'.format(name,t)]
            haz_data[i] = targets['{} at {}'.format(name,t)][0]
            haz_data_lci[i] = haz_data[i] - 1.96*targets['{} at {}'.format(name,t)][1]
            haz_data_uci[i] = haz_data[i] + 1.96*targets['{} at {}'.format(name,t)][1]
            
        
        fig, ax = plt.subplots()
        plt.plot(tval,haz_model*100,'o-k',label='Model')
        plt.plot(tval,haz_data*100,'o--k',label='Data')
        ax.grid(True)
        xticks = tval#[i for i in range(24,36)]
        ax.set_xticks(xticks)
        plt.legend()
        plt.title(cap) 
        plt.xlabel('age')
        plt.ylabel('ratio (%)')
        plt.savefig('{}.pdf'.format(fname))
        
    