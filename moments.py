#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:09:59 2020

@author: egorkozlov
"""
import numpy as np

def compute_moments(self):
    moments = dict()
    
    n_mark = self.state_codes['Couple and child']
    n_marnk = self.state_codes['Couple, no children']
    n_single = self.state_codes['Female, single'] if self.female else self.state_codes['Male, single']
    n_singlek = self.state_codes['Female and child']
    
    
    is_mar = (self.state == n_mark) | (self.state == n_marnk)
    is_mark = (self.state == n_mark)
    
    ever_mar = (np.cumsum(is_mar,axis=1) > 0)
    div_now =  (ever_mar) & ((self.state==n_single) | (self.state==n_singlek))
    ever_div = (np.cumsum(div_now,axis=1) > 0)
    ever_kid = ( np.cumsum( (self.state == n_mark) | (self.state == n_singlek),axis=1) > 0)
    ever_upp = (np.cumsum(self.unplanned_preg,axis=1)>0)
    ever_km  = (np.cumsum(self.k_m,axis=1)>0)
    
    
    
    nmar_cum = np.zeros_like(self.agreed,dtype=np.uint8)
    nmar_cum[:,1:] = np.cumsum(self.agreed[:,:-1],axis=1)
    
    assert np.all(nmar_cum == self.nmar)
    
    one_mar = (nmar_cum == 1)
    
    obs_k_m = self.k_m & one_mar
    obs_m_k = self.m_k & one_mar
    
    
    have_kid = (self.state == n_mark) | (self.state == n_singlek)
    num_mar = np.cumsum( self.agreed, axis = 1 )
    one_mar = (num_mar == 1)
    
    #share_x = self.x / np.maximum(1e-3, self.x  + self.c)  
    #mean_x = share_x[:,0:20][is_mark[:,0:20]].mean()
    
    share_x = self.x[:,:20][is_mark[:,:20]] / self.couple_earnings[:,:20][is_mark[:,:20]]
    mean_x = np.mean(share_x)
    
    ls_min = min(self.setup.ls_levels['Couple and child'])
    
    ls_fem_30 = (self.labor_supply[is_mark[:,9],9]>ls_min).mean()
    
    
    
    age = np.broadcast_to((21+np.arange(self.T)[None,:]),(self.N,self.T))
    
    age_m30 = age - 30
    
    
    
    moments['mean x share'] = mean_x
    moments['in labor force at 30 if kids'] = ls_fem_30
    
    age_pick = ((age>=21) & (age<=40))
    
    never_mar = 1 - ever_mar
    reg_y = never_mar[age_pick]
    reg_x = age_m30[age_pick]
    try:
        pol = np.polyfit(reg_x,reg_y,2)
    except:
        pol = [0,0,0]
    moments['never married by age, b0'] = pol[2]
    moments['never married by age, b1'] = pol[1]
    moments['never married by age, b2'] = pol[0]
    
    
    
    
    
    
    reg_y = ever_kid[age_pick]
    reg_x = age_m30[age_pick]
    try:
        pol = np.polyfit(reg_x,reg_y,2)
    except:
        pol = [0,0,0]
    moments['ever kids by age, b0'] = pol[2]
    moments['ever kids by age, b1'] = pol[1]
    moments['ever kids by age, b2'] = pol[0]
    
    
    yaftmar_pick = (self.yaftmar>=1) & (self.yaftmar<=10)
    pick = yaftmar_pick & age_pick & one_mar & is_mar
    
    if np.any(pick):
        reg_y = ever_kid[pick]
        reg_x = self.yaftmar[pick]
        pol = np.polyfit(reg_x,reg_y,2)
    else:
        pol = [1,1,1]
        
    moments['ever kids by years after marriage, b0'] = pol[2]
    moments['ever kids by years after marriage, b1'] = pol[1]
    moments['ever kids by years after marriage, b2'] = pol[0]
    
    
    for t in range(1,11):
        moments['ever kids by years after marriage, {}'.format(t)] = ever_kid[pick & (self.yaftmar==t)].mean()
    
    
    
    pick = yaftmar_pick & age_pick & one_mar
    if np.any(pick):
        reg_y = div_now[pick]
        reg_x = self.yaftmar[pick]
        pol = np.polyfit(reg_x,reg_y,2)
    else:
        pol = [1,1,1]
        
    moments['divorced by years after marriage, b0'] = pol[2]
    moments['divorced by years after marriage, b1'] = pol[1]
    moments['divorced by years after marriage, b2'] = pol[0]
    
    for t in range(1,11):
        moments['divorced by years after marriage, {}'.format(t)]  = div_now[pick & (self.yaftmar==t)].mean()
    
    
    
    
    # divorce hazards
    # simulating by years after marriage
    pick = age_pick & one_mar
    
    still_mar = np.zeros((11,),dtype=np.float64)
    just_div = np.zeros((11,),dtype=np.float64)
    haz_div = np.zeros((11,),dtype=np.float64)
    for yam in range(0,11):
        thing = pick & (self.yaftmar==yam)
        if np.any(thing):
            still_mar[yam] = is_mar[thing].mean()
            just_div[yam] = self.just_divorced[thing].mean()
            haz_div[yam] = just_div[yam]/still_mar[yam]
    
    

    
    moments['divorced at 30 if one marriage'] = div_now[one_mar[:,9],9].mean()
    
    t = 9
    pick = (self.state[:,t] == n_single)    
    moments['median savings to earnings at 30, single'] = np.median(self.savings_to_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average savings to earnings at 30, single'] = np.mean(self.savings_to_earnings[pick,t]) if np.any(pick) else 0.0
    moments['median consumption to earnings at 30, single']  =  np.median(self.c[pick,t]/self.female_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average consumption to earnings at 30, single']  =  np.mean(self.c[pick,t]/self.female_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average savings at 30, single']  =  np.mean(self.s[pick,t]) if np.any(pick) else 0.0
    
    
    pick = is_mar[:,t]
    moments['median savings to earnings at 30, couples'] = np.median(self.savings_to_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average savings to earnings at 30, couples'] = np.mean(self.savings_to_earnings[pick,t]) if np.any(pick) else 0.0
    moments['median consumption to earnings at 30, couples']  =  np.median(self.c[pick,t]/self.couple_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average consumption to earnings at 30, couples']  =  np.mean(self.c[pick,t]/self.couple_earnings[pick,t]) if np.any(pick) else 0.0
    moments['average savings at 30, couple']  =  np.mean(self.s[pick,t]) if np.any(pick) else 0.0
    
    
    
    
    
    
    moments['never married at 25'] = 1-ever_mar[:,4].mean()
    moments['never married at 30'] = 1-ever_mar[:,9].mean()
    moments['never married at 35'] = 1-ever_mar[:,14].mean()
    moments['never married at 40'] = 1-ever_mar[:,19].mean()
    
    
    
    
    
    
    moments['divorced right now at 25'] = div_now[ever_mar[:,4],4].mean()
    moments['divorced right now at 30'] = div_now[ever_mar[:,9],9].mean()
    moments['divorced right now at 35'] = div_now[ever_mar[:,14],14].mean()
    moments['divorced right now at 40'] = div_now[ever_mar[:,19],19].mean()
    
    
    
    moments['no kids at 25'] = 1-ever_kid[:,4].mean()
    moments['no kids at 30'] = 1-ever_kid[:,9].mean()
    moments['no kids at 35'] = 1-ever_kid[:,14].mean()
    
    
    moments['more than one mar at 40'] = (num_mar[:,19]>1).mean()
    moments['more than one mar at 30'] = (num_mar[:,9]>1).mean()
    
    if np.any(num_mar[:,9]>1):
        moments['ever kids if remarried at 30'] = ever_kid[(num_mar[:,9]>1),9].mean()
    else:
        moments['ever kids if remarried at 30'] = 0.0
    
    
    
    moments['no kids at 25 if married'] = 1-ever_kid[is_mar[:,4],4].mean() if np.any(is_mar[:,4]) else 0.0
    moments['no kids at 30 if married'] = 1-ever_kid[is_mar[:,9],9].mean() if np.any(is_mar[:,9]) else 0.0
    moments['no kids at 35 if married'] = 1-ever_kid[is_mar[:,14],14].mean() if np.any(is_mar[:,14]) else 0.0
    
    
    in_sample = obs_k_m | obs_m_k
    
    
    for i in range(1,15):
        moments['k then m in sample at {}'.format(21+i)] = obs_k_m[in_sample[:,i],i].mean()
    
    
    for i in range(1,15):
        moments['k then m in population at {}'.format(21+i)] = obs_k_m[:,i].mean()
    
    for i in range(1,15):
        moments['m then k in population at {}'.format(21+i)] = obs_m_k[:,i].mean()
   
    
    
    pick = age_pick & in_sample
    
    try:
        reg_y = obs_k_m[pick]
        reg_x = age_m30[pick]
        pol = np.polyfit(reg_x,reg_y,2)
    except:
        pol = [0,0,0]
        
    moments['k then m by age in sample, b0'] = pol[2]
    moments['k then m by age in sample, b1'] = pol[1]
    moments['k then m by age in sample, b2'] = pol[0]
    
    pick = self.agreed & one_mar & age_pick
    just_mark_t0 = self.agreed_unplanned
    
    try:
        reg_y = just_mark_t0[pick]
        reg_x = age_m30[pick]
        pol = np.polyfit(reg_x,reg_y,1)
    except:
        pol = [0,0,0]
        
    moments['share of kids in new marriages at 25'] = just_mark_t0[self.agreed[:,4],4].mean()   if np.any(self.agreed[:,4]) else 0.0
    moments['share of kids in new marriages at 30'] = just_mark_t0[self.agreed[:,9],9].mean()   if np.any(self.agreed[:,9]) else 0.0
    moments['share of kids in new marriages at 35'] = just_mark_t0[self.agreed[:,14],14].mean() if np.any(self.agreed[:,15]) else 0.0
    
    
    
    moments['just k & m at 25'] = (self.agreed_k & one_mar)[:,4].mean()
    moments['just k & m at 30'] = (self.agreed_k & one_mar)[:,9].mean()
    moments['just k & m at 35'] = (self.agreed_k & one_mar)[:,14].mean()
    
    
    inc = self.female_earnings if self.female else self.male_earnings
    inc_30 = inc[:,9]
    pick = (inc_30 >= np.median(inc_30))
    im_pick = is_mar[:,9]    
    em_pick = ever_mar[:,9]
    moments['divorced at 30, above median'] = div_now[pick & em_pick,9].mean()
    moments['ever married at 30, above median'] = ever_mar[pick,9].mean()
    moments['ever kids at 30, above median'] = ever_kid[pick & im_pick,9].mean()
    pick = (inc_30 <= np.median(inc_30))
    moments['divorced at 30, below median'] = div_now[pick & em_pick,9].mean()
    moments['ever married at 30, below median'] = ever_mar[pick,9].mean()
    moments['ever kids at 30, below median'] = ever_kid[pick & im_pick,9].mean()
    
    
    #ls_fem_30 = self.labor_supply[is_mark[:,9],9].mean()
    me_med = np.median(self.male_earnings[is_mark[:,9],9])
    pick_above = (self.male_earnings[:,9] >= me_med) & (is_mark[:,9])
    ls_fem_30_abovemed = (self.labor_supply[pick_above,9]>ls_min).mean()
    pick_below = (self.male_earnings[:,9] <= me_med) & (is_mark[:,9])
    ls_fem_30_below = (self.labor_supply[pick_below,9]>ls_min).mean()
    ls_fem_30_ratio = ls_fem_30_abovemed/ls_fem_30_below if ls_fem_30_below > 0 else 1.0
    moments['in labor force at 30 if kids ratio'] = ls_fem_30_ratio
    moments['in labor force at 30, k then m'] = (self.labor_supply[:,9]>ls_min)[obs_k_m[:,9]].mean() if np.any(self.k_m[:,9]) else 0.0
    moments['in labor force at 30, m then k'] = (self.labor_supply[:,9]>ls_min)[obs_m_k[:,9]].mean() if np.any(self.m_k[:,9]) else 0.0
    
    try:
        moments['divorced at 30, ratio'] = moments['divorced at 30, above median']/moments['divorced at 30, below median']
    except:
        moments['divorced at 30, ratio'] = 0.0
        
    try:
        moments['ever married at 30, ratio'] = moments['ever married at 30, above median'] / moments['ever married at 30, below median'] 
    except:
        moments['ever married at 30, ratio'] = 0.0
        
    try:
        moments['ever kids at 30, ratio'] = moments['ever kids at 30, above median'] / moments['ever kids at 30, below median']
    except:
        moments['ever kids at 30, ratio'] = 0.0
    
    moments['divorced with kids at 30']      = (div_now[:,9]     &  ever_kid[:,9])[ever_mar[:,9]].mean()
    moments['divorced never kids at 30']     = (div_now[:,9]     & ~ever_kid[:,9])[ever_mar[:,9]].mean()
    moments['never married with kids at 30'] = ((~ever_mar)[:,9] &  ever_kid[:,9]).mean()
    
    moments['share of divorced with kids at 30'] = moments['divorced with kids at 30']/(moments['divorced with kids at 30'] + moments['divorced never kids at 30'])
    
    
    share_planned = self.planned_preg[(self.planned_preg) | (self.unplanned_preg)].mean()
    moments['share of planned pregnancies'] = share_planned
    share_rejected = self.disagreed.sum() / (self.disagreed | self.agreed).sum()
    moments['share of rejected proposals'] = share_rejected
    if self.verbose: print('Rejected: {}, planned preg: {}'.format(share_rejected,share_planned))
    
    
    divorced_km = div_now[:,:20][self.k_m[:,:20]].mean()
    divorced_mk = div_now[:,:20][self.m_k[:,:20]].mean()    
    if self.verbose: print('Anything: divorced k_m = {}, divorced m_k = {}'.format(divorced_km,divorced_mk))
    moments['divorced if km (all)'] = divorced_km
    moments['divorced if mk (all)'] = divorced_mk
    
    
    
    
    divorced_km_1m = div_now[:,:20][obs_k_m[:,:20] & one_mar[:,:20]].mean()
    divorced_mk_1m = div_now[:,:20][obs_m_k[:,:20] & one_mar[:,:20]].mean()    
    if self.verbose: print('One mar: divorced k_m = {}, divorced m_k = {}'.format(divorced_km_1m,divorced_mk_1m))
    
    
    moments['divorced if k then m and one marriage, panel'] = divorced_km_1m
    moments['divorced if m then k and one marriage, panel'] = divorced_mk_1m
    
    from age_dist import get_age_distribution
    age_dist = get_age_distribution()
    age_dist_cumulative = np.cumsum(age_dist)
    
    N = div_now.shape[0]
    t_pick = (np.random.random_sample(N)[:,None] > age_dist_cumulative[None,:]).sum(axis=1)
    divorced_km_1m_cs = div_now[np.arange(N),t_pick][obs_k_m[np.arange(N),t_pick]].mean()
    divorced_mk_1m_cs = div_now[np.arange(N),t_pick][obs_m_k[np.arange(N),t_pick]].mean()
    
    moments['divorced if k then m and one marriage'] = divorced_km_1m_cs
    moments['divorced if m then k and one marriage'] = divorced_mk_1m_cs
    
    
    
    e_divorced_upp  = ever_div[ever_upp[:,20],20].mean()
    e_divorced_nupp = ever_div[~ever_upp[:,20],20].mean()    
    if self.verbose: print('Ever divorced upp = {}, ever divorced nupp = {}'.format(e_divorced_upp,e_divorced_nupp))
    moments['ever divorced if had unplanned pregnancy'] = e_divorced_upp
    moments['ever divorced if no unplanned pregnancy'] = e_divorced_nupp
    
    e_divorced_ekm  = ever_div[ever_km[:,20],20].mean()
    e_divorced_nekm = ever_div[~ever_km[:,20],20].mean()    
    if self.verbose: print('Ever divorced ever km = {}, ever divorced never km = {}'.format(e_divorced_ekm,e_divorced_nekm))
    moments['ever divorced if ever km'] = e_divorced_upp
    moments['ever divorced if never km'] = e_divorced_nupp
    
    
    
    def std_pos(x):
        return np.std(np.log(x[x>0]))
    
    
    pick_24 = (self.labor_supply[:,3]>ls_min)
    sd_f_24 = std_pos(self.female_earnings[pick_24,3]) if np.any(pick_24) else 0.0
    pick_30 = (self.labor_supply[:,9]>ls_min)
    sd_f_30 = std_pos(self.female_earnings[pick_30,9]) if np.any(pick_30) else 0.0
    
    sd_m_24 = std_pos(self.male_earnings[:,3])
    sd_m_30 = std_pos(self.male_earnings[:,9])
    
    moments['std earnings at 24, female'] = sd_f_24
    moments['std earnings at 30, female'] = sd_f_30
    moments['std earnings at 24, male'] = sd_m_24
    moments['std earnings at 30, male'] = sd_m_30
    
    
    
    
    if self.verbose:
        print('std of earnings is {} at 24 and {} at 30 for males'.format(sd_m_24,sd_m_30))
        print('std of earnings is {} at 24 and {} at 30 for females'.format(sd_f_24,sd_f_30))
    
    
    i25 = 4
    p25 = (self.male_wage[:,i25] > 0)
    i30 = 9
    p30 = (self.male_wage[:,i30] > 0)
    i40 = 19
    p40 = (self.male_wage[:,i40] > 0)
    
    
    
    
    med_25 = np.median(self.male_wage[p25,i25])
    med_30 = np.median(self.male_wage[p30,i30])    
    
    above_med_25 = ((self.male_wage[:,i25] >= med_25) & p25)
    below_med_25 = ((self.male_wage[:,i25] <= med_25) & p25)
    above_med_30 = ((self.male_wage[:,i30] >= med_30) & p30)
    below_med_30 = ((self.male_wage[:,i30] <= med_30) & p30)
    
    moments['log earnings coef at 25'] = ever_kid[above_med_25,i25].mean() - ever_kid[below_med_25,i25].mean() 
    moments['log earnings coef at 30'] = ever_kid[above_med_30,i30].mean() - ever_kid[below_med_30,i30].mean()
    

    
    p_1yr = (~is_mar[:,0:-2] & is_mar[:,2:] & one_mar[:,2:] & (self.male_wage[:,2:]>0))
    linc_own = np.log(self.female_wage[:,2:][p_1yr] )
    linc_sp =  np.log(self.male_wage[:,2:][p_1yr])
    
    
    c, o, m = self.marriage_stats()
    
    
    n_childless = (~have_kid).sum(axis=0)
    n_newkids = np.zeros_like(n_childless)
    n_newkids[:-1] = ((have_kid[:,1:]) & (~have_kid[:,:-1])).sum(axis=0)
    
    
    
    haz_m = m['Everyone']/c['Everyone']
    haz_mk = (m['Single, pregnant'])/c['Single']
    haz_nk = n_newkids / n_childless 
    
    for t in range(1,15):
        moments['hazard of marriage at {}'.format(21+t)] = haz_m[t] if not np.isnan(haz_m[t]) else 0.0
    for t in range(1,15):
        moments['hazard of marriage & having a child at {}'.format(21+t)] = haz_mk[t] if not np.isnan(haz_mk[t]) else 0.0
        
    for t in range(1,15):
        moments['hazard of new child at {}'.format(21+t)] = haz_nk[t] if not np.isnan(haz_nk[t]) else 0.0
    
    try:
        moments['spouse log coef 1 year after'] = np.polyfit(linc_sp,linc_own,1)[0]
    except:
        moments['spouse log coef 1 year after'] = 0.0
    
    if self.verbose:
        print('Coefficients are {} at 25 and {} at 30'.format(moments['log earnings coef at 25'],moments['log earnings coef at 30']))
    
    
    if self.verbose:
        print('')
        print('')
        print('Key target: km {}, mk {}, ratio {}'.format(divorced_km_1m,divorced_mk_1m,divorced_km_1m/divorced_mk_1m))


    
    return moments

