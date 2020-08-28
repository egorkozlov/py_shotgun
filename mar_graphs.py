#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:19:14 2020

@author: egorkozlov
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma

def mar_graphs(mdl,t=2):
    setup = mdl.setup
    V = mdl.V
    from marriage import v_mar
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
    
    
    inds = mdl.setup.all_indices(t)[0] # all possible outcomes for couples
    iac = np.searchsorted(agrid_c,5.0*wm)
    
    
    icouple = np.broadcast_to(icouple[:,None],(icouple.size,inds.size))
    results = list()
    for upp in [False,True]:
        mt = 'Regular' if not upp else 'Unplanned pregnancy'
        res = v_mar(setup,V[t],t,icouple,inds,match_type=mt,female=True,
                                                        return_insides=True)
    
        res_r = mdl.x_reshape(res['itheta'][...,None],t).squeeze(axis=-1)
        tht_r = setup.thetagrid_fine[res_r]
        tht_all = ma.masked_where(res_r<0,tht_r)
        tht_pick = tht_all[iac,:,izm,:]
        print(tht_pick.shape)
        fig, ax = plt.subplots()
        cs = ax.contourf(zfg,psig,tht_pick.T,cmap='Blues',vmin=0.0,vmax=0.95) 
        cb = fig.colorbar(cs)
        cb.set_label(r'Resulting female bargaining power ($\theta$)')
        plt.xlabel('Female productivity')
        plt.ylabel(r'Love shock at match ($\psi$)',labelpad=-3.0)
        plt.title('Bargaining: {}'.format('Unplanned Pregnancy' if upp else 'Regular Match')) 
        results.append(res)
        #ax.grid(True)
        ax.set_xticks(zfg)
        ax.set_yticks(np.linspace(psig.min(),psig.max(),7))
        plt.savefig(('barg_upp.pdf' if upp else 'barg_reg_match.pdf'))
        
        
    result_noupp, result_upp = results
    
    
    izf = 1
    izm = 3
    ipsi = 5
    
    fig, axs = plt.subplots(1,2)
    
    
    
    fig2, axs2 = plt.subplots()
    
    
    
    
    
    for n, upp, res in zip([0,1],[False,True],[result_noupp,result_upp]):
        Vfm,Vmm,Vfs,Vms = res['insides']
        Vfm0, Vmm0 = [mdl.x_reshape(x,t) for x in [Vfm,Vmm]]
        Vfs0,Vms0 = [mdl.x_reshape(x[...,None],t) for x in [Vfs,Vms]]
        
        
        
        if not upp: 
            norm_f = Vfs0[iac,izf,izm,ipsi]
            norm_m = Vms0[iac,izf,izm,ipsi]
            
        vfm_tht = Vfm0[iac,izf,izm,ipsi,:] - norm_f
        vmm_tht = Vmm0[iac,izf,izm,ipsi,:] - norm_m
        
        
        
        
        vfs_tht = Vfs0[iac,izf,izm,ipsi]*np.ones(vfm_tht.size) - norm_f
        vms_tht = Vms0[iac,izf,izm,ipsi]*np.ones(vmm_tht.size) - norm_m
        
        print((iac,izf,izm,ipsi))
        print(Vfs0.shape)
        print(Vfm0.shape)
        tht_fine = setup.thetagrid
        
        l = 'unplanned pregnancy' if upp else 'regular match'
        c = 'o' if upp else 'x'
        k = 'b' if upp else 'k'

        axs[0].plot(tht_fine,vfm_tht,'{}-{}'.format(c,k),label='Agree, {}'.format(l),markevery=1)
        axs[0].plot(tht_fine,vfs_tht,'{}--{}'.format(c,k),label='Disagree, {}'.format(l),markevery=1)
        axs[1].plot(tht_fine,vmm_tht,'{}-{}'.format(c,k),label='Agree, {}'.format(l),markevery=1)
        axs[1].plot(tht_fine,vms_tht,'{}--{}'.format(c,k),label='Disagree, {}'.format(l),markevery=1)
        axs[1].set_title('Male')
        axs[0].set_title('Female')
        
        axs2.plot(tht_fine,vfm_tht-vfs_tht,'{}-{}'.format(c,k),label='F, {}'.format(l),markevery=1)
        axs2.plot(tht_fine,vmm_tht-vms_tht,'{}--{}'.format(c,k),label='M, {}'.format(l),markevery=1)
        if not upp: axs2.plot(tht_fine,np.zeros_like(tht_fine),'{}:k'.format(c),markevery=1)
        
    axs[0].set_ylabel(r'Value functions (normalized)')
    [ax.set_xlabel(r'Bargaining power ($\theta$)') for ax in axs]
    axs2.set_xlabel(r'Bargaining power ($\theta$)')
    axs[1].legend(bbox_to_anchor=(-1.4, -0.175),loc='upper left',ncol=2)
    fig.suptitle('Change in values b/c of unplanned pregnancy',fontsize=16)
    fig.subplots_adjust(bottom=+0.25)
    fig.savefig('change_upp.pdf',bbox='tight',pad_inches=0.5)
    fig2.suptitle('Surplus over disagreement',fontsize=16)
    axs2.legend(ncol=1)
    axs2.set_xlim(0.0,1.0)
    fig2.savefig('surp_upp.pdf',bbox='tight',pad_inches=0.5)
    #axs[0].set_ylim(-50,10)
    #axs[1].set_ylim(-50,10)
    
        
    

        
    
    
    
    
        
        