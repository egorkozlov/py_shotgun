# -*- coding: utf-8 -*- 
""" 
Created on Wed Dec 18 12:52:29 2019 
 
@author: Fabio 
""" 
 
import pandas as pd 
import numpy as np 
import pickle 
import tabulate 
import statsmodels.api as sm 
import statsmodels.formula.api as smf 
 
################################ 
#Functions 
############################### 
def hazards(dataset,event,duration,end,listh,number,wgt): 
     #Create hazard given some spells in 
     #dataframe 
     
    #Number of unit weights-needed later 
    lgh=np.sum(dataset[wgt])#len(dataset) 
     
    for t in range(number): 
     
        #Get who has the event realized in t+1 
        cond=np.array(dataset[duration])==t+1 
        temp=dataset[[end,wgt]][cond] 
        cond1=temp[end]==event 
        temp1=temp[cond1] 
         
        #Compute the hazard 
        if lgh>0: 
            haz1=np.sum(temp1[wgt])/lgh#len(temp1)/lgh 
            lgh=lgh-np.sum(temp[wgt])#lgh-len(temp) 
        else: 
            haz1=0.0 
         
         
        #If hazard is zero do substitute this 
        #with a random very small number. This 
        #will help later for computing the variance 
        #covariance matrix 
        if (haz1<0.0001): 
            haz1=np.random.uniform(0,0.0001) 
             
        #Add hazard to the list 
        listh=[haz1]+listh 
     
    listh.reverse() 
    listh=np.array(listh).T 
    return listh 
 
 
 
 
##################################### 
#Routine that computes moments 
##################################### 
def compute(hi,period=1): 
    #compute moments, period 
    #says how many years correspond to one 
    #period 
 
    #Get Date at Interview 
    hi.insert(0, 'IDN', range(0,  len(hi))) 
    hi['res']=hi['NUMUNION']+hi['NUMCOHMR'] 
     
    #Get Duration bins 
    bins_d=np.linspace(0,1200,(100/period)+1) 
    bins_d_label=np.linspace(1,len(bins_d)-1,len(bins_d)-1) 
     
    ########################## 
    #Gen cohabitation Dataset 
    ######################### 
     
    #Get date at interview 
    hi['int']=hi['IDATMM']+(hi['IDATYY']-1900)*12 
     
    #Take only if cohabitations 
    coh=hi[(hi['NUMUNION']-hi['NUMMAR']>0) |  (hi['NUMCOHMR']>0)].copy() 
     
     
    #Create number of cohabitations 
    coh['num']=0 
    for i in range(9): 
        coh.loc[coh['HOWBEG0'+str(i+1)]=='coh','num']=coh.loc[coh['HOWBEG0'+str(i+1)]=='coh','num']+1 
             
    #Expand the data     
    cohe=coh.loc[coh.index.repeat(coh.num)] 
     
     
    #Link each cohabitation to relationship number 
    cohe['rell'] = cohe.groupby(['IDN']).cumcount()+1 
    cohe['cou']=1 
    cohe['rel']=None 
    for i in range(9): 
        cohe.loc[(cohe['HOWBEG0'+str(i+1)]=='coh') & (cohe['rell']==cohe['cou']),'rel']=i+1 
        cohe.loc[cohe['HOWBEG0'+str(i+1)]=='coh','cou']= cohe.loc[cohe['HOWBEG0'+str(i+1)]=='coh','cou']+1 
         
    #Get beginning and end of relationhip 
    cohe['beg']=-1 
    cohe['endd']=-1 
    cohe['how']=-1 
    cohe['mar']=-1     
    for i in range(9): 
        cohe.loc[(i+1==cohe['rel']),'beg']=cohe.loc[(i+1==cohe['rel']),'BEGDAT0'+str(i+1)] 
        cohe.loc[(i+1==cohe['rel']),'endd']=cohe.loc[(i+1==cohe['rel']),'ENDDAT0'+str(i+1)] 
        cohe.loc[(i+1==cohe['rel']),'how']=cohe.loc[(i+1==cohe['rel']),'HOWEND0'+str(i+1)] 
        cohe.loc[(i+1==cohe['rel']),'mar']=cohe.loc[(i+1==cohe['rel']),'MARDAT0'+str(i+1)] 
        #add here an indicator of whether it should be unilateral duvorce scenario 
         
    #Get how relationship end 
    cohe['fine']='censored' 
    cohe.loc[cohe['how']=='sep','fine']='sep' 
    cohe.loc[cohe['how']=='div','fine']='mar' 
    cohe.loc[(cohe['how']=='intact') & (cohe['mar']>1),'fine']='mar' 
     
    #Replace censored date if still together 
    cohe['end']=-1 
    cohe.loc[cohe['fine']=='sep','end']=cohe.loc[cohe['fine']=='sep','endd'] 
    cohe.loc[cohe['fine']=='mar','end']=cohe.loc[cohe['fine']=='mar','mar'] 
    cohe.loc[cohe['fine']=='censored','end']=cohe.loc[cohe['fine']=='censored','int'] 
     
    #Duration 
    cohe['dur']=cohe['end']-cohe['beg'] 
     
    #Keep if no error for duration 
    cohe=cohe[(cohe['dur']>0) & (cohe['dur']<2000)] 
     
    #Transform Duration in Years 
    cohe['dury'] = pd.cut(x=cohe['dur'], bins=bins_d,labels=bins_d_label)  
     
    cohe['dury']=cohe['dury'].astype(float)   
     
    #Eliminate non useful things 
    del coh 
     
    ########################## 
    #Gen marriage Dataset 
    ######################### 
     
    #Take only if marriages 
    mar=hi[hi['NUMMAR']>0].copy() 
     
    #Create number of cohabitations 
    mar['num']=0 
    for i in range(9): 
        mar.loc[mar['MARDAT0'+str(i+1)]>0,'num']=mar.loc[mar['MARDAT0'+str(i+1)]>0,'num']+1 
             
    #Expand the data     
    mare=mar.loc[mar.index.repeat(mar.num)] 
     
     
    #Link each marriage to relationship number 
    mare['rell'] = mare.groupby(['IDN']).cumcount()+1 
    mare['cou']=1 
    mare['rel']=None 
    for i in range(9): 
        mare.loc[(mare['MARDAT0'+str(i+1)]>0) & (mare['rell']==mare['cou']),'rel']=i+1 
        mare.loc[mare['MARDAT0'+str(i+1)]>0,'cou']= mare.loc[mare['MARDAT0'+str(i+1)]>0,'cou']+1 
         
    #Get beginning and end of relationhip 
    mare['beg']=-1 
    mare['endd']=-1 
    mare['how']=-1 
    mare['mar']=-1     
    for i in range(9): 
        mare.loc[(i+1==mare['rel']),'beg']=mare.loc[(i+1==mare['rel']),'MARDAT0'+str(i+1)] 
        mare.loc[(i+1==mare['rel']),'endd']=mare.loc[(i+1==mare['rel']),'ENDDAT0'+str(i+1)] 
        mare.loc[(i+1==mare['rel']),'how']=mare.loc[(i+1==mare['rel']),'HOWEND0'+str(i+1)] 
     
         
    #Get how relationship end 
    mare['fine']='censored' 
    mare.loc[mare['how']=='div','fine']='div' 
     
     
    #Replace censored date if still together 
    mare['end']=-1 
    mare.loc[mare['fine']=='div','end']=mare.loc[mare['fine']=='div','endd'] 
    mare.loc[mare['fine']=='censored','end']=mare.loc[mare['fine']=='censored','int'] 
     
    #Duration 
    mare['dur']=mare['end']-mare['beg'] 
     
    #Keep if no error for duration 
    mare=mare[(mare['dur']>0) & (mare['dur']<2000)] 
     
    #Transform Duration in Years 
    mare['dury'] = pd.cut(x=mare['dur'], bins=bins_d,labels=bins_d_label)  
     
    mare['dury']=mare['dury'].astype(float)  
     
    del mar 
     
    ############################# 
    #Build relationship by month 
    ############################## 
     
    #Eliminate observation if info on beg-end not complete 
    #for i in range(9): 
     #   hi=hi[(np.isfinite(hi['BEGDAT0'+str(i+1)])) & (hi['BEGDAT0'+str(i+1)]<3999)] 
         
    #Get date in time at which the guy is 20,25...,50 (9) 
    for j in range(7): 
        hi['time_'+str(20+(j)*5)]=hi['DOBY']*12+hi['DOBM']+(20+(j)*5)*12 
         
    #Get the status 
    for j in range(7): 
         
        #Create the variable of Status 
        hi['status_'+str(20+(j)*5)]='single' 
         
 
         
        for i in range(9): 
             
             
            #Get if in couple 
            hi.loc[(hi['time_'+str(20+(j)*5)]>=hi['BEGDAT0'+str(i+1)]) & (hi['BEGDAT0'+str(i+1)]<3999) & 
                   (((hi['time_'+str(20+(j)*5)]<=hi['ENDDAT0'+str(i+1)]) & (hi['ENDDAT0'+str(i+1)]>0))  | 
                    (hi['ENDDAT0'+str(i+1)]==0) | (hi['WIDDAT0'+str(i+1)]>0) )
                   ,'status_'+str(20+(j)*5)]='mar' 
                    
                   
            #Substitute if actually cohabitation  
            hi.loc[(hi['time_'+str(20+(j)*5)]>=hi['BEGDAT0'+str(i+1)]) & (hi['BEGDAT0'+str(i+1)]<3999) & 
                   (((hi['time_'+str(20+(j)*5)]<=hi['ENDDAT0'+str(i+1)]) & (hi['ENDDAT0'+str(i+1)]>0))  |  
                    (hi['ENDDAT0'+str(i+1)]==0) | (hi['WIDDAT0'+str(i+1)]>0) ) &  
                    (hi['status_'+str(20+(j)*5)]=='mar') &  
                   (hi['HOWBEG0'+str(i+1)]=='coh')    &  
                   ((hi['MARDAT0'+str(i+1)]==0) | (hi['MARDAT0'+str(i+1)]>hi['time_'+str(20+(j)*5)]))      
                   ,'status_'+str(20+(j)*5)]='coh'  
                    
    #Create the variables ever cohabited and ever married 
    for j in range(7): 
         
        #Create the variable of ever married or cohabit 
        hi['everm_'+str(20+(j)*5)]=0.0 
        hi['everc_'+str(20+(j)*5)]=0.0 
         
        for i in range(9): 
                    
            #Get if ever cohabited  
            hi.loc[((hi['everc_'+str(20+(max(j-1,0))*5)]>=0.1) | ((hi['HOWBEG0'+str(i+1)]=='coh') & (hi['time_'+str(20+(j)*5)]>=hi['BEGDAT0'+str(i+1)]))),'everc_'+str(20+(j)*5)]=1.0 
             
            #Get if ever cohabited  
            hi.loc[((hi['everm_'+str(20+(max(j-1,0))*5)]>=0.1) |  (hi['time_'+str(20+(j)*5)]>=hi['MARDAT0'+str(i+1)])),'everm_'+str(20+(j)*5)]=1.0 
             
    ###################################### 
    #Build employment by status in 1986 
    ###################################### 
    empl=hi[(hi['M2DP01']=='FEMALE') & (hi['weeks']<99)].copy() 
    empl['stat']='single' 
    empl['dist']=99999 
    for j in range(7): 
        empl.loc[np.abs(empl['time_'+str(20+(j)*5)]-86*12)<empl['dist'],'stat']=hi['status_'+str(20+(j)*5)] 
             
    ########################## 
    #BUILD HAZARD RATES 
    #########################  
     
    #Hazard of Separation 
    hazs=list() 
    hazs=hazards(cohe,'sep','dury','fine',hazs,int(6/period),'SAMWT') 
     
    #Hazard of Marriage 
    hazm=list() 
    hazm=hazards(cohe,'mar','dury','fine',hazm,int(6/period),'SAMWT') 
     
    #Hazard of Divorce 
    hazd=list() 
    hazd=hazards(mare,'div','dury','fine',hazd,int(6/period),'SAMWT') 
     
    ######################################## 
    #Construct share of each relationship 
    ####################################### 
    mar=np.zeros(7) 
    coh=np.zeros(7) 
    emar=np.zeros(7) 
    ecoh=np.zeros(7) 
     
    for j in range(7): 
        mar[j]=np.average(hi['status_'+str(20+(j)*5)]=='mar', weights=np.array(hi['SAMWT'])) 
        coh[j]=np.average(hi['status_'+str(20+(j)*5)]=='coh', weights=np.array(hi['SAMWT'])) 
        emar[j]=np.average(hi['everm_'+str(20+(j)*5)], weights=np.array(hi['SAMWT'])) 
        ecoh[j]=np.average(hi['everc_'+str(20+(j)*5)], weights=np.array(hi['SAMWT'])) 
         
         
    ######################################### 
    # Femle Labor Supply 
    ######################################### 
    fls_ratio=np.average(empl.loc[empl['stat']=='mar','work'], weights=np.array(empl.loc[empl['stat']=='mar','SAMWT']))/np.average(empl.loc[empl['stat']=='coh','work'], weights=np.array(empl.loc[empl['stat']=='coh','SAMWT'])) 
         
     
    ######################################### 
    #Create the age at unilateral divorce+ 
    #regression on the effect of unilateral divorce 
    ########################################### 
     
    #Number of relationships for the person 
    hi['numerl']=0.0 
     
    #List of variables to keep 
    keep_var=list() 
    keep_var=keep_var+['numerl']+['state']+['SAMWT'] 
     
    for i in range(9): 
         
        #Make sure that some relationship of order i exist 
        if (np.any(hi['BEGDAT0'+str(i+1)])): 
             
            #Add relationship order 
            hi['order'+str(i+1)]=np.nan 
            hi.loc[np.isnan(hi['BEGDAT0'+str(i+1)])==False,'order'+str(i+1)]=i+1 
             
            #Add number of relationships 
            hi.loc[np.isnan(hi['BEGDAT0'+str(i+1)])==False,'numerl']+=1.0 
             
            #Get whether the relationship started in marriage or cohabitation 
            hi['imar'+str(i+1)]=np.nan 
            hi.loc[hi['HOWBEG0'+str(i+1)]=='coh','imar'+str(i+1)]=0.0 
            hi.loc[hi['HOWBEG0'+str(i+1)]=='mar','imar'+str(i+1)]=1.0 
             
            #Get age at relationship 
            hi['iage'+str(i+1)]=np.nan 
            hi.loc[np.isnan(hi['BEGDAT0'+str(i+1)])==False,'iage'+str(i+1)]=round((hi['BEGDAT0'+str(i+1)]-hi['birth_month'])/12) 
             
            #Get if unilateral divorce when relationship started 
            hi['unid'+str(i+1)]=np.nan 
            hi.loc[np.isnan(hi['BEGDAT0'+str(i+1)])==False,'unid'+str(i+1)]=0.0 
            hi.loc[(round(hi['BEGDAT0'+str(i+1)]/12+1900)>=hi['unil']) & (hi['unil']>0.1),'unid'+str(i+1)]=1.0 
             
            #Year Realationship Started 
            hi['year'+str(i+1)]=np.nan 
            hi.loc[np.isnan(hi['BEGDAT0'+str(i+1)])==False,'year'+str(i+1)]=round(hi['BEGDAT0'+str(i+1)]/12+1900) 
                        
            #Keep variables 
            keep_var=keep_var+['year'+str(i+1)]+['unid'+str(i+1)]+['iage'+str(i+1)]+['imar'+str(i+1)]+['order'+str(i+1)] 
     
         
         
    #New Dataset to reshape 
    hi2=hi[keep_var] 
     
    #Reahspe Dataset 
    years = ([col for col in hi2.columns if col.startswith('year')]) 
    unids = ([col for col in hi2.columns if col.startswith('unid')]) 
    iages = ([col for col in hi2.columns if col.startswith('iage')]) 
    imars = ([col for col in hi2.columns if col.startswith('imar')]) 
    order = ([col for col in hi2.columns if col.startswith('order')]) 
     
    hi3 = pd.lreshape(hi2, {'year' : years,'unid' : unids,'iage' : iages,'imar' : imars,'order' : order})  
     
    #Eliminate if missing 
    hi3.replace([np.inf, -np.inf], np.nan) 
    hi3.dropna(subset=['imar','unid']) 
     
    #Regression 
    FE_ols = smf.wls(formula='imar ~ unid+C(order)+C(iage)+C(state)+C(year)',weights=hi3['SAMWT'], data = hi3.dropna()).fit() 
    beta_unid=FE_ols.params['unid'] 
     
    #Get age at which unilateral divorced was introduced 
    hi['age_unid']=0.0 
    hi.loc[hi['unil']==0,'age_unid']=1000.0 
    hi.loc[hi['unil']!=0,'age_unid']=hi['unil']-hi['birth']   
       
     
    #Get age in the second survey 
    date_age=pd.read_csv('age_drop.csv') 
     
    #From hi make '-1' if law changed before the guy starts 
    hi.loc[hi['age_unid']<0,'age_unid']=-1 
     
     
    #Function to get frequencies 
    def CountFrequency(my_list):  
   
        # Creating an empty dictionary   
        freq = {}  
        for item in my_list:  
            if (item in freq):  
                freq[item] += 1 
            else:  
                freq[item] = 1 
       
        #for key, value in freq.items():  
         #   print ("% d : % d"%(key, value))  
         
        return freq 
         
    #Frequencies for age policy change 
    freq_pc = CountFrequency(hi['age_unid'].tolist()) 
     
    #Frequencies for age in the second wave 
    freq_i= CountFrequency(date_age['age'].tolist()) 
     
    #Drop if errors 
     
    #hi3.to_csv(r'D:\Downloads\temp.csv') 
     
    #Create a dictionary for saving simulated moments 
    listofTuples = [("hazs" , hazs), ("hazm" , hazm),("hazd" , hazd),("emar" , emar),
                    ("ecoh" , ecoh), ("fls_ratio" , fls_ratio),("mar" , mar),("coh" , coh),
                    ("freq_pc" , freq_pc), ("freq_i" , freq_i),("beta_unid" , beta_unid)] 
    dic_mom=dict(listofTuples) 
     
    del hi,hi2,hi3 
    return dic_mom 
 
 
 
############################################## 
#Actual moments computation + weighting matrix 
################################################ 
 
def dat_moments(sampling_number=4,weighting=False,covariances=False,period=1): 
     
     
     
         
    #Import Data 
    data=pd.read_csv('histo.csv') 
     
    #Call the routine to compute the moments 
    dic=compute(data.copy(),period=period) 
    hazs=dic['hazs'] 
    hazm=dic['hazm'] 
    hazd=dic['hazd'] 
    emar=dic['emar'] 
    ecoh=dic['ecoh'] 
    fls_ratio=dic['fls_ratio'] 
    mar=dic['mar'] 
    coh=dic['coh'] 
    freq_pc=dic['freq_pc'] 
    freq_i=dic['freq_i'] 
    beta_unid=dic['beta_unid'] 
     
     
    #Use bootstrap samples to compute the weighting matrix 
    n=len(data) 
    boot=sampling_number 
    nn=n*boot 
     
    hazsB=np.zeros((len(hazs),boot)) 
    hazmB=np.zeros((len(hazm),boot)) 
    hazdB=np.zeros((len(hazd),boot)) 
    marB=np.zeros((len(mar),boot)) 
    cohB=np.zeros((len(coh),boot)) 
    emarB=np.zeros((len(emar),boot)) 
    ecohB=np.zeros((len(ecoh),boot)) 
    fls_ratioB=np.zeros((1,boot)) 
    beta_unidB=np.zeros((1,boot)) 
     
    aa=data.sample(n=nn,replace=True,weights='SAMWT',random_state=4) 
     
    #Make weights useless, we already used them for sampling 
    aa['SAMWT']=1 
    for i in range(boot): 
     
        a1=aa[(i*n):((i+1)*n)].copy().reset_index() 
        dicti=compute(a1.copy(),period=period) 
        hazsB[:,i]=dicti['hazs'] 
        hazmB[:,i]=dicti['hazm'] 
        hazdB[:,i]=dicti['hazd'] 
        emarB[:,i]=dicti['emar'] 
        ecohB[:,i]=dicti['ecoh'] 
        fls_ratioB[:,i]=dicti['fls_ratio'] 
        marB[:,i]=dicti['mar'] 
        beta_unidB[:,i]=dicti['beta_unid'] 
     
         
     
    ################################# 
    #Confidence interval of moments 
    ################################ 
    hazmi=np.array((np.percentile(hazmB,2.5,axis=1),np.percentile(hazmB,97.5,axis=1))) 
    hazsi=np.array((np.percentile(hazsB,2.5,axis=1),np.percentile(hazsB,97.5,axis=1))) 
    hazdi=np.array((np.percentile(hazdB,2.5,axis=1),np.percentile(hazdB,97.5,axis=1))) 
    mari=np.array((np.percentile(marB,2.5,axis=1),np.percentile(marB,97.5,axis=1))) 
    cohi=np.array((np.percentile(cohB,2.5,axis=1),np.percentile(cohB,97.5,axis=1))) 
    emari=np.array((np.percentile(emarB,2.5,axis=1),np.percentile(emarB,97.5,axis=1))) 
    ecohi=np.array((np.percentile(ecohB,2.5,axis=1),np.percentile(ecohB,97.5,axis=1))) 
    fls_ratioi=np.array((np.percentile(fls_ratioB,2.5,axis=1),np.percentile(fls_ratioB,97.5,axis=1))) 
    beta_unidi=np.array((np.percentile(beta_unidB,2.5,axis=1),np.percentile(beta_unidB,97.5,axis=1))) 
     
    #Do what is next only if you want the weighting matrix    
    if weighting: 
         
        #Compute optimal Weighting Matrix 
        col=np.concatenate((hazmB,hazsB,hazdB,emarB,ecohB,fls_ratioB,beta_unidB),axis=0)     
        dim=len(col) 
        W_in=np.zeros((dim,dim)) 
        for i in range(dim): 
            for j in range(dim): 
                W_in[i,j]=(1/(boot-1))*np.cov(col[i,:],col[j,:])[0][1] 
               
        if not covariances: 
            W_in = np.diag(np.diag(W_in)) 
         
        #Invert 
        W=np.linalg.inv(W_in) 
         
        # normalize 
        W = W/W.sum() 
         
    else: 
         
        
        #If no weighting, just use sum of squred deviations as the objective function         
        W=np.diag(np.ones(len(hazm)+len(hazs)+len(hazd)+len(emar)+len(ecoh)+2))#two is for fls+beta_unid 
        
            
            
            
    listofTuples = [("hazs" , hazs), ("hazm" , hazm),("hazd" , hazd),("emar" , emar), 
                ("ecoh" , ecoh), ("fls_ratio" , fls_ratio),("mar" , mar),("coh" , coh), 
                ("beta_unid" , beta_unid),
                ("hazsi" , hazsi), ("hazmi" , hazmi),("hazdi" , hazdi),("emari" , emari),
                ("ecohi" , ecohi), ("fls_ratioi" , fls_ratioi),("mari" , mari),("cohi" , cohi),
                ("beta_unidi" , beta_unidi),("W",W)] 
    packed_stuff=dict(listofTuples) 
    #packed_stuff = (hazm,hazs,hazd,emar,ecoh,fls_ratio,W,hazmi,hazsi,hazdi,emari,ecohi,fls_ratioi,mar,coh,mari,cohi) 
     
     
    #Export Moments 
    with open('moments.pkl', 'wb+') as file: 
        pickle.dump(packed_stuff,file)   
         
    #Export Age at Unilateral Divorce 
    with open('age_uni.pkl', 'wb+') as file: 
        pickle.dump(freq_pc,file)  
         
    #Export Age at second  
    with open('age_sw.pkl', 'wb+') as file: 
        pickle.dump(freq_i,file)  
         
    del packed_stuff,freq_i,freq_pc,aa,data 
         
################################################################### 
#If script is run as main, it performs a data comparison with SIPP 
################################################################### 
if __name__ == '__main__': 
     
    import matplotlib.pyplot as plt 
    import matplotlib.backends.backend_pdf 
    import pickle 
 
 
    #Get stuff about moments 
    dat_moments(period=1) 
    print(11111111111111)
    ########################## 
    #Import and work SIPP data 
    ########################## 
     
    #SAmples-variables over year 
    samples=('08','04','01','96') 
     
    for j in samples: 
        print(j) 
        name_sample='sipp'+j 
        vars()[name_sample]=pd.read_stata('D:/blasutto/Data/SIPP raw/sipp'+j+'t.dta') 
         
        #Keep if always observed 
        vars()[name_sample]=vars()[name_sample].dropna(thresh=2) 
 
 
        name_date='date'+j 
        vars()[name_date] = np.ones(16)*np.nan 
        name_mar='mar'+j 
        vars()[name_mar] = np.ones(16)*np.nan 
        name_mar1='mar1'+j 
        vars()[name_mar1] = np.ones(16)*np.nan 
        name_coh='coh'+j 
        vars()[name_coh] = np.ones(16)*np.nan 
         
 
        for i in range(16): 
             
            #Create date 
            vars()[name_date][i]=np.max(vars()[name_sample][['date'+str(i+1)]])[0] 
             
            #Create Marriage and cohabitation rates 
 
            try: 
                vars()[name_mar][i]=np.average(vars()[name_sample].loc[vars()[name_sample]['married_'+str(i+1)]>=0,'married_'+str(i+1)], weights=np.array(vars()[name_sample].loc[vars()[name_sample]['married_'+str(i+1)]>=0,'wpfinwgt']))#*2.0 
                vars()[name_mar1][i]=np.average(vars()[name_sample].loc[vars()[name_sample]['mara_'+str(i+1)]>=0,'mara_'+str(i+1)], weights=np.array(vars()[name_sample].loc[vars()[name_sample]['mara_'+str(i+1)]>=0,'wpfinwgt']))*2.0 
                vars()[name_coh][i]=np.average(vars()[name_sample].loc[vars()[name_sample]['cohab_'+str(i+1)]>=0,'cohab_'+str(i+1)], weights=np.array(vars()[name_sample].loc[vars()[name_sample]['cohab_'+str(i+1)]>=0,'wpfinwgt']))*2.0 
            except: 
                pass 
     
     
    #SAmples-variables over age 
    merged=pd.DataFrame.append(pd.DataFrame.append(pd.DataFrame.append(sipp01,sipp04,sort=False),sipp08,sort=False),sipp96,sort=False) 
    age = np.linspace(20, 60, 41) 
    marr_age=np.ones(41)*np.nan 
    coh_age=np.ones(41)*np.nan 
     
     
     
    def share_age(sample,boot=False): 
        #Function that computes share 
        #of married and cohabiting by age 
        marr_age_s=np.ones(41)*np.nan 
        coh_age_s=np.ones(41)*np.nan  
         
        if boot: 
            sample['wpfinwgt']=1 
     
        for i in age: 
            try: 
                marr_age_s[int(i-21)]=np.average(sample.loc[((abs(sample['age1'])-i)<=0.1),'married_1'],weights=sample.loc[((abs(sample['age1'])-i)<=0.1),'wpfinwgt']) 
                coh_age_s[int(i-21)]=np.average(sample.loc[((abs(sample['age1'])-i)<=0.1),'cohab_1'],weights=sample.loc[((abs(sample['age1'])-i)<=0.1),'wpfinwgt'])*2.0 
            except: 
                pass 
        return marr_age_s,coh_age_s 
     
    #Compute acutal share of cohabiting and married over time 
    marr_age,coh_age=share_age(merged,boot=False) 
     
    ####################### 
    #Compute CI of SIPP data 
    ######################## 
    n=len(merged) 
    boot=100 
    nn=n*boot 
     
    marr_ageB=np.zeros((len(marr_age),boot)) 
    coh_ageB=np.zeros((len(coh_age),boot)) 
    sampling=merged.sample(n=nn,replace=True,weights='wpfinwgt',random_state=4) 
 
    for i in range(boot): 
 
       samp_small=sampling[(i*n):((i+1)*n)].copy().reset_index() 
       marr_ageB[:,i],coh_ageB[:,i]=share_age(samp_small.copy(),boot=True) 
        
        
    marr_agei=np.array((np.percentile(marr_ageB,2.5,axis=1),np.percentile(marr_ageB,97.5,axis=1))) 
    coh_agei=np.array((np.percentile(coh_ageB,2.5,axis=1),np.percentile(coh_ageB,97.5,axis=1))) 
 
    #################### 
    #Get NLSFH data      
    ##################### 
    with open('moments.pkl', 'rb') as file: 
        packed_data=pickle.load(file) 
    #datanlsh=np.array[()] 
      
    #Unpack Moments (see data_moments.py to check if changes) 
    #(hazm,hazs,hazd,mar,coh,fls_ratio,W) 
    hazm_d=packed_data['hazm'] 
    hazs_d=packed_data['hazs'] 
    hazd_d=packed_data['hazd'] 
    mar_d=packed_data['emar'] 
    coh_d=packed_data['ecoh'] 
    fls_d=np.ones(1)*packed_data['fls_ratio'] 
    beta_unid_d=np.ones(1)*packed_data['beta_unid'] 
    hazm_i=packed_data['hazmi'] 
    hazs_i=packed_data['hazsi'] 
    hazd_i=packed_data['hazdi'] 
    mar_i=packed_data['emari'] 
    coh_i=packed_data['ecohi'] 
    fls_i=np.ones(1)*packed_data['fls_ratioi'] 
    beta_unid_i=np.ones(1)*packed_data['beta_unidi'] 
             
    #Create Graph-Marriage 
   
     
    #Create Cohabitation Over Age 
    fig1 = plt.figure() 
    age_d=np.linspace(20,50,7) 
    plt.plot(age_d, coh_d,'r',linewidth=1.5, label='Share Cohabiting NLSFH') 
    plt.fill_between(age_d, coh_i[0,:], coh_i[1,:],alpha=0.2,facecolor='r') 
    plt.plot(age, coh_age,'g',linewidth=1.5, label='Share Cohabiting SIPP') 
    plt.fill_between(age, coh_agei[0,:], coh_agei[1,:],alpha=0.2,facecolor='g') 
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), 
              fancybox=True, shadow=True, ncol=2, fontsize='x-small') 
    plt.ylim(ymax=0.1) 
    plt.ylim(ymin=0.0) 
    plt.xlabel('Age') 
    plt.ylabel('Share') 
     
     
    #Create Marriage Over Age 
    fig2 = plt.figure() 
    age_d=np.linspace(20,50,7) 
    plt.plot(age_d, mar_d,'r',linewidth=1.5, label='Share Married NLSFH') 
    plt.fill_between(age_d, mar_i[0,:], mar_i[1,:],alpha=0.2,facecolor='r') 
    plt.plot(age, marr_age,'g',linewidth=1.5, label='Share Married SIPP') 
    plt.fill_between(age, marr_agei[0,:], marr_agei[1,:],alpha=0.2,facecolor='g') 
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), 
              fancybox=True, shadow=True, ncol=2, fontsize='x-small') 
    plt.ylim(ymax=1.0) 
    plt.ylim(ymin=0.0) 
    plt.xlabel('Age') 
    plt.ylabel('Share') 
     
     
 
 
 
 
 
 
     
 
 
 
         
         
 
 
 
 
 
