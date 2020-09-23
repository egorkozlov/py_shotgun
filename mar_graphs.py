#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:19:14 2020

@author: egorkozlov
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma

try:
    import cupy as cp
except:
    pass

def mar_graphs(self,t=2):
    mdl = self
    setup = mdl.setup
    V = mdl.V
    from marriage import v_mar
    from gridvec import VecOnGrid
    agrid_c = mdl.setup.agrid_c
    agrid_s = mdl.setup.agrid_s
    
    
    npsi = setup.pars['n_psi_t'][t]
    psig = setup.exogrid.psi_t[t]
    nzf = setup.pars['n_zf_t'][t]
    zfg = setup.exogrid.zf_t[t]
    zmg = setup.exogrid.zm_t[t]
    
    
    
    izm = 2
    
    wm = np.exp(setup.pars['m_wage_trend'][t] + zmg[izm])
    wf = np.exp( setup.pars['f_wage_trend'][t] + zfg) 
    
    
    inds = mdl.setup.all_indices(t)[0] # all possible outcomes for couples
    
    
    af = 0.0*wf[3]
    am = 0.0*wm
    
    iaf = np.searchsorted(agrid_c,af)
    iac = np.searchsorted(agrid_c,af+am)
    
    icouple = iac*np.ones_like(agrid_c,dtype=np.int16)
    icouple = np.broadcast_to(icouple[:,None],(icouple.size,inds.size))
    results = list()
    
    
    
    
    for upp, ss in zip([False,True,True],[True,True,False]):
        mt = 'Regular' if not upp else 'Unplanned pregnancy'
        Vnext = self.cupyfy_v(V[t]) if self.gpu else V[t]
        
        if not ss:
            f = setup.pars['disutil_shotgun']
            setup.pars['disutil_shotgun'] = 0.0
        
        res = v_mar(setup,Vnext,t,icouple,inds,match_type=mt,female=True,
                                                        return_insides=True)
        
        if not ss: setup.pars['disutil_shotgun'] = f
        
        if self.gpu: res = {k:cp.asnumpy(res[k]) for k in res}
        results.append(res)

    
    
    for upp in [False,True]:
        mt = 'Regular' if not upp else 'Unplanned pregnancy'
        
        res = results[0] if not upp else results[1]
        
        
        res_r = mdl.x_reshape(res['itheta'][...,None],t).squeeze(axis=-1)
        tht_r = setup.thetagrid_fine[res_r]
        tht_all = ma.masked_where(res_r<0,tht_r)
        tht_pick = tht_all[iaf,:,izm,:].astype(np.float32)
        fig, ax = plt.subplots()
        cs = ax.contourf(zfg,psig,tht_pick.T,cmap='Blues',vmin=0.05,vmax=0.95) 
        cb = fig.colorbar(cs)
        cb.set_label(r'Resulting female bargaining power ($\theta$)')
        plt.xlabel('Female productivity')
        plt.ylabel(r'Love shock at match ($\psi$)',labelpad=-3.0)
        plt.title('Bargaining: {}'.format('Unplanned Pregnancy' if upp else 'Regular Match')) 
        #ax.grid(True)
        ax.set_xticks(zfg)
        ax.set_yticks(np.linspace(psig.min(),psig.max(),7))
        plt.savefig(('barg_upp.pdf' if upp else 'barg_reg_match.pdf'))
        
        
    result_noupp, result_upp, result_upp_noss = results
    
    
    izf = 3
    izm = 2
    ipsi = 7
    
    fig, axs = plt.subplots(1,2)
    
    
    
    fig2, axs2 = plt.subplots()
    
    fig3, axs3 = plt.subplots()
    
    
    
    
    
    
    
    for n, upp, res in zip([0,2,1],[False,True,True],[result_noupp,result_upp_noss,result_upp]):
    #for n, upp, res in zip([0,1],[False,True],[result_noupp,result_upp]):
        Vfm,Vmm,Vfs,Vms = res['V_f_yes'], res['V_m_yes'], res['V_f_no'], res['V_m_no']        
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
        
        if n == 2:
            l = 'UP, but no stigma'
            k, c = 'g', 'd'

        if n < 2: axs[0].plot(tht_fine,vfm_tht,'{}-{}'.format(c,k),label='Agree, {}'.format(l),markevery=1)
        axs[0].plot(tht_fine,vfs_tht,'{}--{}'.format(c,k),label='Disagree, {}'.format(l),markevery=1)
        if n < 2: axs[1].plot(tht_fine,vmm_tht,'{}-{}'.format(c,k),label='Agree, {}'.format(l),markevery=1)
        axs[1].plot(tht_fine,vms_tht,'{}--{}'.format(c,k),label='Disagree, {}'.format(l),markevery=1)
        axs[1].set_title('Male')
        axs[0].set_title('Female')
        
        axs2.plot(tht_fine,vfm_tht-vfs_tht,'{}-{}'.format(c,k),label='F, {}'.format(l),markevery=1)
        axs2.plot(tht_fine,vmm_tht-vms_tht,'{}--{}'.format(c,k),label='M, {}'.format(l),markevery=1)
        if not upp: axs2.plot(tht_fine,np.zeros_like(tht_fine),'{}:k'.format(c),markevery=1)
        
        
        surp_f = Vfm0[iac,izf,izm,:,:] - Vfs0[iac,izf,izm,:]
        surp_m = Vfm0[iac,izf,izm,:,:] - Vfs0[iac,izf,izm,:]
        
        surp = np.minimum(surp_f,surp_m).max(axis=1) 
        
        axs3.plot(psig,surp,'{}-{}'.format(c,k),label='Surplus, {}'.format(l))
        if not upp: axs3.plot(psig,np.zeros_like(psig),'--')
        axs3.set_ylim(-20,20)
        
        
        
        
        
    axs[0].set_ylabel(r'Value functions (normalized)')
    [ax.set_xlabel(r'Bargaining power ($\theta$)') for ax in axs]
    axs2.set_xlabel(r'Bargaining power ($\theta$)')
    lgd = axs[1].legend(bbox_to_anchor=(-1.4, -0.175),loc='upper left',ncol=2)
    txt = fig.suptitle('Change in values b/c of unplanned pregnancy',fontsize=16)
    fig.subplots_adjust(bottom=+0.25)
    fig.savefig('change_upp.pdf',bbox_extra_artists=(txt,lgd), bbox_inches='tight') #bbox='tight',pad_inches=0.5)
    fig2.suptitle('Surplus over disagreement',fontsize=16)
    lgd = axs2.legend(ncol=1)
    axs2.set_xlim(0.0,1.0)
    fig2.savefig('surp_upp.pdf',bbox_extra_artists=(lgd,), bbox_inches='tight') #bbox='tight',pad_inches=0.5)
    
    axs3.set_title('Defining the agreement threshold:\n point where the surplus hits zero',fontsize=16)
    axs3.legend(ncol=1)
    fig3.savefig('surp_upp.pdf',bbox='tight',pad_inches=0.5)
    axs3.set_ylabel(r'Egalitarian Bargaining Surplus $M(\psi)$')
    axs3.set_xlabel(r'Match quality ($\psi$)')
    
    
    
    
    
    
    
    
    
    fig, axs = plt.subplots()
    
    for n, upp, res in zip([0,1,2],[False,True,True],[result_noupp,result_upp,result_upp_noss]):
        
        l = 'unplanned pregnancy' if upp else 'regular match'
        if n == 2: l = 'UP, but no stigma'
        c = 'o' if upp else 'x'
        k = 'b' if upp else 'k'
        if n == 2: k, c = 'g', 'd'
        
        Vfm,Vmm,Vfs,Vms = res['V_f_yes'], res['V_m_yes'], res['V_f_no'], res['V_m_no']        
        Vfm0, Vmm0 = [mdl.x_reshape(x,t) for x in [Vfm,Vmm]]
        Vfs0,Vms0 = [mdl.x_reshape(x[...,None],t) for x in [Vfs,Vms]]
        surplus = np.minimum(Vfm0-Vfs0,Vmm0-Vms0).max(axis=4)
        ind_last_neg = np.sum((surplus < 0),axis=3,keepdims=True)-1
        last_neg_v = np.take_along_axis(surplus,ind_last_neg,3)
        first_pos_v = np.take_along_axis(surplus,ind_last_neg+1,3)
        w_pos = (-last_neg_v) / (first_pos_v - last_neg_v)
        assert np.all(w_pos<=1)
        assert np.all(w_pos>=0)
        
        tholds = psig[ind_last_neg+1]*w_pos + psig[ind_last_neg]*(1-w_pos)
        # this is an array of thresholds, dimensionality is (a,zf,zm), but 
        # a is quirky
        
        axs.plot(zfg,tholds[iac,:,izm],'{}-{}'.format(c,k),label='{}'.format(l))
        
    axs.set_ylabel(r'Agreement threshold for relationship quality $\psi_m$')
    axs.set_xlabel(r'Female productivity $z^f$')
    txt = axs.set_title('Agreement thresholds for marriage \n by unplanned pregnancy status',fontsize=16)
    axs.legend()
    fig.savefig('tholds_upp.pdf',bbox_extra_artists=(txt,),bbox='tight',pad_inches=2.0)
    plt.show()
        
    
    
    
    
    

        
    
    
    
    
        
        