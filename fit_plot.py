#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:32:21 2020

@author: egorkozlov
"""



import matplotlib.pyplot as plt
from targets import all_targets
from tiktak import filer
from setup import ModelSetup


import numpy as np

class FitPlots(object):
    
    def __init__(self,*,base,targ_mode,compare=None,
                 setup=None,
                 base_name='Model',
                 compare_name='Data',
                 graphs_title_add=None,
                 moments_aux=None,
                 noplot=False):
        
        
        if setup is None: self.setup = ModelSetup()
        
        
        if type(base) is str:
            self.moments = filer(base,0,0,repeat=False)
        else:
            self.moments = base
            
            
                    
            
        self.targ_mode = targ_mode
        self.base_name = base_name
        self.compare_name = compare_name
        
        if graphs_title_add:
            self.graphs_title_add = '\n ' + graphs_title_add
        else:
            self.graphs_title_add = ''
        
            
        if type(compare) is str:
            targ_load = filer(compare,0,0,repeat=False)
            # emulating targets
            self.targets = {key: (targ_load[key],0.0) for key in targ_load}
        else:
            self.targets = all_targets(targ_mode)
        
        try:
            self.print_things()
        except:
            print('failed to print')
            
        self.moments_aux = moments_aux
        
        if not noplot: self.plot()
    
    def plot(self):
        '''
        try:
            self.plot_estimates()
        except:
            print('failed to plot estimates')
        
        try:
            self.plot_hazards()
        except:
            print('failed to plot hazards')            
        '''
        
        
        self.plot_trends()
        
        try:
            self.plot_cumulative()
        except:
            print('failed to plot cumulative')
            
        try:
            self.plot_by_years_after_marriage()
        except:
            print('failed to plot by years after marriage')
            
        try:
            self.plot_kfmf()
        except:
            print('failed to plot kfmf')
            
        if self.moments_aux is not None: 
            try:
                self.plot_men()
            except:
                print('failed to plot men')
        
        try:
            self.plot_single_moms()
        except:
            print('failed to plot single moms')
            
            
        try:
            self.plot_single_moms_new()
        except:
            print('failed to plot single moms (new)')
                    
        
        try:
            self.plot_welfare()
        except:
            print('failed to plot welfare')
        
        try:
            self.plot_kfmf_ref()
        except:
            print('failed to plot ref')
            
        plt.show()
        
    
    def print_things(self):
        labels = ['divorced if k then m and one marriage',
                  'divorced by years after marriage if kids first, 10',
                  'divorced if m then k and one marriage',
                  'divorced by years after marriage if marriage first, 10']
        for l in labels:
            print('{}: base {}, compare {}'.format(l,self.moments[l],self.targets[l][0]))
    
    
    
    def plot_trends(self):
        setup = self.setup
        tval = np.arange(21,45)
        
        trend_f = np.zeros_like(tval,dtype=np.float32)
        trend_m = np.zeros_like(tval,dtype=np.float32)
        for i, t in enumerate(tval):
            trend_f[i] = np.exp( setup.pars['f_wage_trend'][i] )
            trend_m[i] = np.exp( setup.pars['m_wage_trend'][i-2] ) if i >= 2 else None
        
        
        fig, ax = plt.subplots()
        ax.plot(tval,trend_f,'--b',label='female')
        ax.plot(tval,trend_m,'-k',label='male')

        ax.set_xlabel('age')
        ax.set_ylabel('income (USD 1000s)')
        ax.set_title('estimated earnings trends:\n{}'.format(self.base_name))
        ax.legend()
        ax.grid(True)
        plt.savefig('earnings_trend.pdf')

        
        
        
        
        
        
    
    def plot_estimates(self):
        setup = self.setup
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
        plt.title('estimated probabilities' + self.graphs_title_add)
        plt.xlabel('age')
        plt.ylabel('probability (%)')
            
        
        
    
        
    
    
    
    
    def plot_hazards(self,ci=False):
        moments = self.moments
        targets = self.targets
        setup = self.setup
        # graph 1: hazard of any marriage
        
        
        tval = np.arange(23,36)  
        
        names = ['hazard of marriage','hazard of new child']
        captions = ["hazard of marriage:\n (% new marriages at [age]) / (% single at [age-1])"  + self.graphs_title_add,
                    "hazard of new child:\n (% new births at [age]) / (% childless at [age-1])"  + self.graphs_title_add]
        
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
            plt.plot(tval,haz_model*100,'o-b',label=self.base_name)
            plt.plot(tval,haz_data*100,'o-k',label=self.compare_name)
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
            
    
    def plot_cumulative(self):
        # graph 1: hazard of any marriage
        moments,targets,setup = self.moments, self.targets, self.setup
        import matplotlib.gridspec as gridspec
        
        tval = np.arange(23,36)  
        
        names = ['k then m in population','m then k in population','k then m in sample','ever married']
        captions = ["Kids First" ,
                    "Marriage First",
                    "Relative % of Kids First by age:\n (Kids First) / (Kids First + Marriage First)" + self.graphs_title_add,
                    "ever married" + self.graphs_title_add]
        
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
        
        
        ch = ['o','x']
        nm = ['KF','MF']
        for i in range(2):
            #if i < 2:
            ax = fig.add_subplot(gs[0,0])
            #else:
            #    ax = fig.add_subplot(gs[:,1])
            ax.plot(tval,probs_model[i],'{}-b'.format(ch[i]),label='{}: {}'.format(nm[i],self.base_name.lower()))
            ax.plot(tval,probs_data[i],'{}-k'.format(ch[i]),label='{}: {}'.format(nm[i],self.compare_name.lower()))
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('% in female population by age:\nKids First and Marriage First' + self.graphs_title_add)
        ax.legend()
        #yticks = [i*5 for i in range(12)]
        #ax.set_yticks(yticks)
        xticks = [i for i in range(22,36)]
        ax.set_xticks(xticks)
        ax.grid(True)
        plt.ylim(bottom=-1.0)
        plt.savefig('popshares.pdf')
        
        fig, ax = plt.subplots()
        ax.plot(tval,probs_model[2],'o-b',label=self.base_name)
        ax.plot(tval,probs_data[2],'o-k',label=self.compare_name)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title(captions[2])
        ax.legend()
        ax.set_xticks(xticks)
        ax.grid(True)
        plt.ylim(bottom=-1.0)
        plt.savefig('relshares.pdf')
        
        
        fig, ax = plt.subplots()
        ax.plot(tval,probs_model[3],'o-b',label=self.base_name)
        ax.plot(tval,probs_data[3],'o-k',label=self.compare_name)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title(captions[3])
        ax.legend()
        ax.set_xticks(xticks)
        ax.grid(True)
        plt.savefig('evermar.pdf')
        
    
    def plot_by_years_after_marriage(self):
        # graph 1: hazard of any marriage
        
        moments, targets = self.moments, self.targets
        
        yval = np.arange(1,11) 
        
        names = ['ever kids by years after marriage','divorced by years after marriage']
        captions = ['% with kids by years after marriage\n (if married)' + self.graphs_title_add,'% divorced by years after marriage \n (excluding remarried)'+ self.graphs_title_add]
        fnames = ['kids_by_duration','divorced_by_duration']
        
        
        for name, cap, fname in zip(names,captions, fnames):
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
            plt.plot(yval[i_model],prob_model[i_model]*100,'o-b',label=self.base_name)
            plt.plot(yval[i_data],prob_data[i_data]*100,'d-k',label=self.compare_name)
            #if name == 'ever kids by years after marriage' and self.compare_name == 'data':
            #    plt.plot(yval[i_data][3:],prob_data[i_data][3:]*100,'+',label='targeted',linewidth=6,markersize=14)
            #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
            plt.legend()
            plt.title(cap) 
            plt.xlabel('years after marriage')
            plt.ylabel('share (%)')
            ax.grid(True)
            ax.set_xticks(yval)
            plt.ylim(bottom=-1.0)
            plt.savefig('{}.pdf'.format(fname))
            
            
    
    def plot_kfmf(self):
        # graph 1: hazard of any marriage
        
        moments,targets,setup = self.moments,self.targets,self.setup
        
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
            plt.plot(yval[i_model],prob_model[i_model]*100,'{}-b'.format(mrkr),label='{}: {}'.format(cap,self.base_name.lower()))
            plt.plot(yval[i_data],prob_data[i_data]*100,'{}-k'.format(mrkr),label='{}: {}'.format(cap,self.compare_name.lower()))
            #if name_aux is not None: plt.plot(tval,aux,label=name_aux)
        plt.legend()
        plt.title('Divorced by years after marriage by groups' + self.graphs_title_add) 
        plt.xlabel('years after marriage')
        plt.ylabel('share (%)')
        ax.grid(True)
        ax.set_xticks(yval)
        plt.savefig('div_kfmf.pdf')
        
        # ref_kf = [0.08114558, 0.10186335, 0.11817027, 0.11984021, 0.11823362, 0.11746032, 0.13162393, 0.15057915, 0.14123007, 0.14550265]
        # ref_mf = [0.        , 0.01083856, 0.02031424, 0.02864134, 0.03768433, 0.04414536, 0.05390435, 0.05976757, 0.0625543 , 0.06301748]
        
    def plot_kfmf_ref(self):
        
        assert False
        
        moments,targets,setup,targ_mode = self.moments, self.targets, self.setup, self.targ_mode
        
        yval = np.arange(1,11) 
        #ref_kf = np.array([0.08114558, 0.10186335, 0.11817027, 0.11984021, 0.11823362, 0.11746032, 0.13162393, 0.15057915, 0.14123007, 0.14550265])
        #ref_mf = np.array([0.        , 0.01083856, 0.02031424, 0.02864134, 0.03768433, 0.04414536, 0.05390435, 0.05976757, 0.0625543 , 0.06301748])
        
        if targ_mode == 'high education':
            ref_kf = np.array([0.03697876, 0.07397737, 0.09510086, 0.12460401, 0.1457346 ,
           0.16467463, 0.18292683, 0.18918919, 0.19642857, 0.21967963])
            
            ref_mf = np.array([0.        , 0.00997299, 0.01585933, 0.02122016, 0.0277903 ,
           0.03765505, 0.04619769, 0.05458448, 0.06238917, 0.06841011])
        else:
            ref_kf = np.array([0.02801336, 0.05332585, 0.08019741, 0.10971467, 0.13539367,
           0.1536    , 0.17518248, 0.19001832, 0.20896249, 0.23069054])
            
            ref_mf = np.array([0.        , 0.01145147, 0.0232318 , 0.0351265 , 0.04966411,
           0.06725495, 0.08349328, 0.10104955, 0.11824745, 0.13481709])
        
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
        
        lbl_cmp = '$\phi_s = 0$'
        #lbl_cmp = 'equal pay'
        plt.plot(yval,ref_kf*100,'x-b',label='KF: {}'.format(self.base_name.lower()))
        plt.plot(yval,here_kf*100,'o-b',label='KF: {}'.format(lbl_cmp))
        plt.plot(yval,ref_mf*100,'x-k',label='MF: {}'.format(self.base_name.lower()))
        plt.plot(yval,here_mf*100,'o--b',label='MF: {}'.format(lbl_cmp))
        plt.legend()
        plt.title('Divorced by years after marriage: social stigma matters') 
        #plt.title('Divorced by years after marriage: \n removing gender pay gap') 
        plt.xlabel('years after marriage')
        plt.ylabel('share (%)')
        ax.grid(True)
        ax.set_xticks(yval)
        plt.savefig('div_kfmf_ref.pdf')
        
    
    def plot_single_moms(self):
        # graph 1: hazard of any marriage
        
        moments,targets,setup = self.moments,self.targets,self.setup
        
        tval = np.arange(23,36)
        
        
        nm_sm_model = np.zeros_like(tval,dtype=np.float64)
        div_sm_model = np.zeros_like(tval,dtype=np.float64)
        all_sm_model = np.zeros_like(tval,dtype=np.float64)
        
        nm_sm_data = np.zeros_like(tval,dtype=np.float64)
        div_sm_data = np.zeros_like(tval,dtype=np.float64)
        all_sm_data = np.zeros_like(tval,dtype=np.float64)
        
        div_childless_data =  np.zeros_like(tval,dtype=np.float64)
        div_childless_model = np.zeros_like(tval,dtype=np.float64)
        
        nm_childless_data = np.zeros_like(tval,dtype=np.float64)
        nm_childless_model = np.zeros_like(tval,dtype=np.float64)
        
        for i,t in enumerate(tval):
            
            nm_sm_model[i] = moments['never married and kids in population at {}'.format(t)]
            div_sm_model[i] = moments['divorced and kids in population at {}'.format(t)]
            all_sm_model[i] = nm_sm_model[i] + div_sm_model[i]
            
            nm_sm_data[i] = targets['never married and kids in population at {}'.format(t)][0]
            div_sm_data[i] = targets['divorced and kids in population at {}'.format(t)][0]
            all_sm_data[i] = nm_sm_data[i] + div_sm_data[i]
            
            nm_childless_data[i] = targets['never married and no kids in population at {}'.format(t)][0]
            nm_childless_model[i] = moments['never married and no kids in population at {}'.format(t)]
            
            div_childless_data[i] = targets['divorced and no kids in population at {}'.format(t)][0]
            div_childless_model[i] = moments['divorced and no kids in population at {}'.format(t)]
            
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*all_sm_model,'o-b',label=self.base_name)
        ax.plot(tval,100*all_sm_data,'o-k',label=self.compare_name)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Single mothers in population')
        ax.legend()
        ax.set_xticks(tval)
        plt.ylim(bottom=-1.0)
        ax.grid(True)
        plt.savefig('fit_single_mothers_all.pdf')
        
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*div_sm_model,'o-b',label=self.base_name)
        ax.plot(tval,100*div_sm_data,'d-k',label=self.compare_name)
        if self.compare_name == 'data': ax.plot(tval[5:],100*div_sm_data[5:],'+k',label='targeted',linewidth=6,markersize=14)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Divorced with kids in population')
        ax.legend()
        ax.set_xticks(tval)
        plt.ylim(bottom=-0.2)
        ax.grid(True)
        plt.savefig('fit_single_mothers_div.pdf')
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*nm_sm_model,'o-b',label=self.base_name)
        ax.plot(tval,100*nm_sm_data,'d-k',label=self.compare_name)
        if self.compare_name == 'data': ax.plot(tval[5:],100*nm_sm_data[5:],'+k',label='targeted',linewidth=6,markersize=14)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        plt.ylim(bottom=-0.2)
        ax.set_title('Never married with kids in population')
        ax.legend()
        ax.set_xticks(tval)
        ax.grid(True)
        plt.savefig('fit_single_mothers_nm.pdf')
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*div_childless_model,'o-b',label=self.base_name)
        ax.plot(tval,100*div_childless_data,'d-k',label=self.compare_name)
        if self.compare_name == 'data': ax.plot(tval[5:],100*div_childless_data[5:],'+k',label='targeted',linewidth=6,markersize=14)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Divorced witout kids in population')
        ax.legend()
        ax.set_xticks(tval)
        plt.ylim(bottom=-0.2)
        ax.grid(True)
        plt.savefig('fit_single_nonmothers_div.pdf')
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*nm_childless_model,'o-b',label=self.base_name)
        ax.plot(tval,100*nm_childless_data,'d-k',label=self.compare_name)
        if self.compare_name == 'data': ax.plot(tval[0:],100*nm_childless_data[0:],'+k',label='targeted',linewidth=6,markersize=14)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Never married without kids in population')
        ax.legend()
        ax.set_xticks(tval)
        plt.ylim(bottom=-1.0)
        ax.grid(True)
        plt.savefig('fit_single_nonmothers_nm.pdf')
        
        
        
    def plot_single_moms_new(self):
        # graph 1: hazard of any marriage
        
        moments,targets,setup = self.moments,self.targets,self.setup
        
        tval = np.arange(25,36)
        
        
        
        all_sm_model = np.zeros_like(tval,dtype=np.float64)
        rel_sm_model = np.zeros_like(tval,dtype=np.float64)
        abortions_model = np.zeros_like(tval,dtype=np.float64)
        
        
        
        all_sm_data = np.zeros_like(tval,dtype=np.float64)
        rel_sm_data = np.zeros_like(tval,dtype=np.float64)
        abortions_data = np.zeros_like(tval,dtype=np.float64)
        
        
        for i,t in enumerate(tval):
            
            all_sm_model[i] = moments['single mothers in total at {}'.format(t)]
            all_sm_data[i] = targets['single mothers in total at {}'.format(t)][0]
            
            rel_sm_model[i] = moments['single mothers among mothers at {}'.format(t)]
            rel_sm_data[i] = targets['single mothers among mothers at {}'.format(t)][0]
            
            all_sm_model[i] = moments['single mothers in total at {}'.format(t)]
            all_sm_data[i] = targets['single mothers in total at {}'.format(t)][0]
            
            abortions_model[i] = moments['aborted unplanned pregnancies at {}'.format(t)]
            abortions_data[i] = targets['aborted unplanned pregnancies at {}'.format(t)][0]
            
            
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*all_sm_model,'o-b',label=self.base_name)
        ax.plot(tval,100*all_sm_data,'o-k',label=self.compare_name)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Single mothers in population')
        ax.legend()
        ax.set_xticks(tval)
        ax.grid(True)
        
        
        
        
        fig, ax = plt.subplots()
        ax.plot(tval,100*rel_sm_model,'o-b',label=self.base_name+', abortions {:02.3g}%'.format(np.nanmean(abortions_model)))
        ax.plot(tval,100*rel_sm_data,'o-k',label=self.compare_name+', abortions {:02.3g}%'.format(np.nanmean(abortions_data)))
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Single mothers among mothers + abortions')
        ax.legend()
        ax.set_xticks(tval)
        ax.grid(True) 
        plt.savefig('sm_among_mothers.pdf')

        
        fig, ax = plt.subplots()
        ax.plot(tval,abortions_model,'o-b',label=self.base_name)
        ax.plot(tval,abortions_data,'o-k',label=self.compare_name)
        ax.set_xlabel('age')
        ax.set_ylabel('share (%)')
        ax.set_title('Share of abortions')
        ax.legend()
        ax.set_xticks(tval)
        ax.grid(True)
        
        
    
    def plot_welfare(self,a_mult=35):
        # singles
        z_fem = self.setup.exogrid.zf_t[0]
        w_fem = np.exp(z_fem + self.setup.pars['f_wage_trend'][0])
        
        z_mal = self.setup.exogrid.zm_t[0]
        w_mal = np.exp(z_mal + self.setup.pars['f_wage_trend'][0])
        
        moments,targets = self.moments, self.targets
        
        
        
        # males and females
        a_female = \
                  v_compare(self.setup.agrid_s,
                            moments['value function: female, single, all assets'],
                            targets['value function: female, single, all assets'][0],
                            a_mult=a_mult)
        
        
        a_male = \
                  v_compare(self.setup.agrid_s,
                  moments['value function: male, single, all assets'],
                  targets['value function: male, single, all assets'][0],
                  a_mult=a_mult)
                  
        
        
        fig, ax = plt.subplots()
        plt.plot(z_fem,a_female,'o-b',label='{} - {}'.format(self.base_name,self.compare_name))
        ax.grid(True)
        xticks = z_fem
        ax.set_xticks(xticks)
        plt.legend()
        plt.title('Welfare differences: female, 21, no kids, no assets') 
        plt.xlabel('female productivity')
        plt.ylabel('asset equivalent variation (2016 USD, 1000s)') 
        plt.savefig('welfare_female.pdf')
        
        fig, ax = plt.subplots()        
        plt.plot(z_mal,a_male,'o-b',label='{} - {}'.format(self.base_name,self.compare_name))
        ax.grid(True)
        xticks = z_mal
        ax.set_xticks(xticks)
        plt.legend()
        plt.title('Welfare differences: male, 23, no kids, no assets') 
        plt.xlabel('male productivity')
        plt.ylabel('asset equivalent variation (2016 USD, 1000s)') 
        plt.savefig('welfare_male.pdf')
        
        # couples 
        # transofrm the state
        ipsi = np.arange(self.setup.pars['n_psi'])
        imal = 2*np.ones_like(ipsi)
        ifem = 3*np.ones_like(ipsi)
        
        iexo, _, _, _ = self.setup.all_indices(0,(ifem,imal,ipsi))
        
        
        a_couple = \
                  v_compare(self.setup.agrid_c,
                            moments['value function: couple, no children, all assets'][:,iexo],
                            targets['value function: couple, no children, all assets'][0][:,iexo],
                            a_mult=a_mult)
        
        
        
        psi = self.setup.exogrid.psi_t[0]
        
        
        fig, ax = plt.subplots()
        #plt.plot(psi,v_couple_base_val[iexo],'o-b',label=self.base_name)
        #if v_couple_compare_val is not None: plt.plot(psi,v_couple_compare_val[iexo],'o-k',label=self.compare_name)
        plt.plot(psi,a_couple,'o-b',label='{} - {}'.format(self.base_name,self.compare_name))
        ax.grid(True)
        xticks = np.linspace(psi.min(),psi.max(),7)
        ax.set_xticks(xticks)
        plt.legend()
        plt.title('Welfare comparison: couple, 21/23, no kids, no assets') 
        plt.xlabel('love shock')
        plt.ylabel('asset equivalent variation (2016 USD, 1000s)')  
        
    
        
    def plot_men(self):
        
        moments,targets = self.moments, self.targets
        
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
            plt.plot(tval,haz_model*100,'o-b',label=self.base_name)
            plt.plot(tval,haz_data*100,'o-k',label=self.compare_name)
            ax.grid(True)
            xticks = tval#[i for i in range(24,36)]
            ax.set_xticks(xticks)
            plt.legend()
            plt.title(cap) 
            plt.xlabel('age')
            plt.ylabel('ratio (%)')
            plt.savefig('{}.pdf'.format(fname))
    
    
    
from interp_np import interp
def v_compare(agrid,v_base,v_compare,*,a_mult):
    
    ia0 = 0
    
    names = ['extra assets in base to reach 0 in compare',
             'extra assets in compare to reach 0 in base']
    
    i = 0
    
    aall = list()
    
    amin = agrid.min()
    amax = agrid.max()
    
    for (v0, v1), name in zip([(v_base, v_compare),(v_compare, v_base)],
                                names):
        
        print('')
        print(name)
        print('')
        print('v0 is larger in {} cases'.format(np.sum(v0>v1)))
        print('v0 is smaller in {} cases'.format(np.sum(v0<v1)))
        
        alist = list()
        
        for iz in range(v0.shape[1]):
            
            
            j, wn = interp(v0[:,iz],v1[ia0,iz],return_wnext=True)
            
            if j <= -1:
                a = 0.0
            else:
                a = agrid[j]*(1-wn) + agrid[j+1]*wn
                
            a = np.clip(a,amin,amax)
            
            alist.append(a_mult*a)
            
        
        i += 1
        
        aall.append(np.array(alist))
    return aall[1] - aall[0]
    