#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:19:14 2020

@author: egorkozlov
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma

def mar_graphs(mdl,t=1):
    setup = mdl.setup
    V = mdl.V
    from marriage import v_mar_igrid
    from gridvec import VecOnGrid
    agrid_c = mdl.setup.agrid_c
    agrid_s = mdl.setup.agrid_s
    icouple = VecOnGrid(agrid_c,2*agrid_s).i
    
    npsi = setup.pars['n_psi_t'][t]
    psig = setup.exogrid.psi_t[t]
    nzf = setup.pars['n_zf_t'][t]
    zfg = setup.exogrid.zf_t[t]
    zmg = setup.exogrid.zm_t[t]
    
    
    
    izm = 2
    wm = setup.pars['m_wage_trend'][t] + zmg[izm]
    wzf = np.exp( setup.pars['f_wage_trend'][t] + zfg) / np.exp(wm)
    
    '''
    
    print(psig)
    izf = 3
    izm = 2
    inds_ = (izf*np.ones(npsi,dtype=np.int16),izm*np.ones(npsi,dtype=np.int16),np.arange(npsi,dtype=np.int16))
    inds = mdl.setup.all_indices(t,inds_)[0]
    
    '''
    inds = mdl.setup.all_indices(t)[0]
    iac = np.searchsorted(agrid_c,2.0*wm)
    ias = np.searchsorted(agrid_s,1.0*wm)
    results = list()
    for upp in [False,True]:
        uloss = 0.0 #setup.pars['disutil_shotgun'] if upp else 0.0
        res = v_mar_igrid(setup,t,V[t],icouple,inds,female=True,giveabirth=upp,
                                          unplanned_pregnancy=upp,
                                          uloss_fem=0.0,uloss_mal=0.0,
                                          uloss_fem_single=uloss,uloss_mal_single=uloss,
                                          return_all=True)
    
        res_r = mdl.x_reshape(res['theta'][...,None],t).squeeze(axis=-1)
        tht_r = setup.thetagrid_fine[res_r]
        #tht_v[res['theta'] < 0] = None
        tht_all = ma.masked_where(res_r<0,tht_r)
        tht_pick = tht_all[iac,:,izm,:]
        print(tht_pick.shape)
        fig, ax = plt.subplots()
        cs = ax.contourf(zfg,psig,tht_pick.T,cmap='Blues',vmin=0.0,vmax=0.8) 
        cb = fig.colorbar(cs)
        cb.set_label(r'Resulting female bargaining power ($\theta$)')
        plt.xlabel('Female productivity')
        plt.ylabel(r'Love shock at match ($\psi$)',labelpad=-3.0)
        plt.title('Bargaining: {}'.format('Unplanned Pregnancy (No Stigma)' if upp else 'Regular Match')) 
        results.append(res)
        #ax.grid(True)
        ax.set_xticks(zfg)
        ax.set_yticks(psig)#np.arange(-4,5))
        plt.savefig(('barg_upp.pdf' if upp else 'barg_reg_match.pdf'))
        
        
    result_noupp, result_upp = results
    
    
    izf = 3
    izm = 3
    ipsi = 9
    
    fig, axs = plt.subplots(1,2)
    fig2, axs2 = plt.subplots()
    
    
    
    
    
    for n, upp, res in zip([0,1],[False,True],[result_noupp,result_upp]):
        Vfm,Vmm,Vfs,Vms,gamma = res['ins']
        Vfm0, Vmm0 = [mdl.x_reshape(x,t) for x in [Vfm,Vmm]]
        Vfs0,Vms0 = [mdl.x_reshape(x[...,None],t) for x in [Vfs,Vms]]
        
        
        
        if not upp: 
            norm_f = Vfs0[iac,izf,izm,ipsi]
            norm_m = Vms0[iac,izf,izm,ipsi]
            
        vfm_tht = Vfm0[iac,izf,izm,ipsi,:] - norm_f
        vmm_tht = Vmm0[iac,izf,izm,ipsi,:] - norm_m
        
        
        
        
        vfs_tht = Vfs0[iac,izf,izm,ipsi]*np.ones(vfm_tht.size) - norm_f
        vms_tht = Vms0[iac,izf,izm,ipsi]*np.ones(vmm_tht.size) - norm_m
        
        
        
        tht_fine = setup.thetagrid_fine
        
        l = 'unplanned pregnancy' if upp else 'regular match'
        c = 'o' if upp else 'x'

        axs[0].plot(tht_fine,vfm_tht,'{}-k'.format(c),label='Agree, {}'.format(l),markevery=25)
        axs[0].plot(tht_fine,vfs_tht,'{}--k'.format(c),label='Disagree, {}'.format(l),markevery=25)
        axs[1].plot(tht_fine,vmm_tht,'{}-k'.format(c),label='Agree, {}'.format(l),markevery=25)
        axs[1].plot(tht_fine,vms_tht,'{}--k'.format(c),label='Disagree, {}'.format(l),markevery=25)
        axs[1].set_title('Male')
        axs[0].set_title('Female')
        
        axs2.plot(tht_fine,vfm_tht-vfs_tht,'{}-k'.format(c),label='F, {}'.format(l),markevery=25)
        axs2.plot(tht_fine,vmm_tht-vms_tht,'{}--k'.format(c),label='M, {}'.format(l),markevery=25)
        if not upp: axs2.plot(tht_fine,np.zeros_like(tht_fine),':k'.format(c),markevery=25)
        
    axs[0].set_ylabel(r'Value functions (normalized)')
    [ax.set_xlabel(r'Bargaining power ($\theta$)') for ax in axs]
    axs[1].legend(bbox_to_anchor=(-1.4, -0.175),loc='upper left',ncol=2)
    fig.suptitle('Change in values b/c of unplanned pregnancy',fontsize=16)
    fig.subplots_adjust(bottom=+0.25)
    fig.savefig('change_upp.pdf',bbox='tight',pad_inches=0.5)
    fig2.suptitle('Surplus over disagreement',fontsize=16)
    axs2.legend(ncol=1)
    axs2.set_xlim(0.0,1.0)
    #axs[0].set_ylim(-50,10)
    #axs[1].set_ylim(-50,10)
    
        
    

        
    
    
    
    
        
        