#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

        
def target_values(mode='high education'):
    targets = dict()
        
        
    if mode=='high education':
        
        #targets['ever kids by years after marriage, 1'] = (.1901175,.0018521)
        #targets['ever kids by years after marriage, 2'] = (.3416996,.0022054)
        targets['ever kids by years after marriage, 3'] = (.4902692,.0023415)
        targets['ever kids by years after marriage, 4'] = (.6114712,.0023096)
        targets['ever kids by years after marriage, 5'] = (.7039068,.0021802)
        targets['ever kids by years after marriage, 6'] = (.7730147,.0020489)
        targets['ever kids by years after marriage, 7'] = (.8182496,.0019268)
        targets['ever kids by years after marriage, 8'] = (.8568491,.00179)
        #targets['ever kids by years after marriage, 9'] = (.8809587,.001696)
        #targets['ever kids by years after marriage, 10']= (.8978599,.0016307)
        
        
        targets['divorced by years after marriage, 1']  =  (.0084912,.0004315)
        targets['divorced by years after marriage, 2']  =  (.0199536,.0006448)
        targets['divorced by years after marriage, 3']  =  (.0345649,.0008425)        
        targets['divorced by years after marriage, 4']  =  (.0536216,.001042)
        targets['divorced by years after marriage, 5']  =  (.0683989,.0011694)
        targets['divorced by years after marriage, 6']  =  (.0807483,.0012872)
        targets['divorced by years after marriage, 7']  =  (.0879516,.0013595)
        targets['divorced by years after marriage, 8']  =  (.1011097,.0014718*(1/2))
        targets['divorced by years after marriage, 9']  =  (.1027806,.001517*(1/2))
        targets['divorced by years after marriage, 10'] =  (.116314,.0016379*(1/2))
        
        
        
        
        #targets['hazard of marriage at 22'] = (.0506949,.001038)
        targets['hazard of marriage at 23'] = (.0624283,.0011034)
        targets['hazard of marriage at 24'] = (.0787830,.0012318)
        targets['hazard of marriage at 25'] = (.0945611,.0013867)
        targets['hazard of marriage at 26'] = (.1066170,.0015351)
        targets['hazard of marriage at 27'] = (.1105425,.0016506)
        targets['hazard of marriage at 28'] = (.1181793,.0018265)
        targets['hazard of marriage at 29'] = (.1185538,.0019239)
        targets['hazard of marriage at 30'] = (.1116017,.001993)
        targets['hazard of marriage at 31'] = (.1070343,.0020717)
        targets['hazard of marriage at 32'] = (.1068933,.0022018)
        targets['hazard of marriage at 33'] = (.0913318,.0021093)
        targets['hazard of marriage at 34'] = (.0844846,.0020832)
        targets['hazard of marriage at 35'] = (.0780720,.0020774)
        
        
        #targets['k then m in population at 22'] = (.0031368,.0003093)
        targets['k then m in population at 23'] = (.0055518,.0003391)
        targets['k then m in population at 24'] = (.0080609,.000383)
        targets['k then m in population at 25'] = (.0126187,.0004615)
        targets['k then m in population at 26'] = (.0177664,.0005398)
        targets['k then m in population at 27'] = (.0214169,.0005799)
        targets['k then m in population at 28'] = (.0281133,.0006556)
        targets['k then m in population at 29'] = (.0321021,.0006961)
        targets['k then m in population at 30'] = (.0369590,.000731)
        targets['k then m in population at 31'] = (.0386226,.0007545)
        targets['k then m in population at 32'] = (.0424404,.0007881)
        targets['k then m in population at 33'] = (.0436167,.0008042)
        targets['k then m in population at 34'] = (.0483524,.0008466)
        targets['k then m in population at 35'] = (.0535521,.0008891)
        
        

        
        #targets['m then k in population at 22'] = (.0058368,.0004213)
        targets['m then k in population at 23'] = (.0150030,.0005547)
        targets['m then k in population at 24'] = (.0279987,.0007066)
        targets['m then k in population at 25'] = (.0480545,.0008843)
        targets['m then k in population at 26'] = (.0788468,.0011012)
        targets['m then k in population at 27'] = (.1209867,.0013062)
        targets['m then k in population at 28'] = (.1736392,.0015023)
        targets['m then k in population at 29'] = (.2316302,.0016661)
        targets['m then k in population at 30'] = (.2820957,.0017437)
        targets['m then k in population at 31'] = (.3449800,.0018613)
        targets['m then k in population at 32'] = (.3962855,.0019121)
        targets['m then k in population at 33'] = (.4331152,.0019509)
        targets['m then k in population at 34'] = (.4618895,.0019675)
        targets['m then k in population at 35'] = (.4719915,.0019714)
        
        
        #targets['k then m in sample at 22'] = (.3495573,.0264092)
        targets['k then m in sample at 23'] = (.2700996,.0134549)
        targets['k then m in sample at 24'] = (.2235426,.008985)
        targets['k then m in sample at 25'] = (.2079784,.0063642)
        targets['k then m in sample at 26'] = (.1838923,.0048233)
        targets['k then m in sample at 27'] = (.1503958,.0036041)
        targets['k then m in sample at 28'] = (.1393456,.0029425)
        targets['k then m in sample at 29'] = (.1217223,.0024172)
        targets['k then m in sample at 30'] = (.1158390,.0021046)
        targets['k then m in sample at 31'] = (.1006838,.0018461)
        targets['k then m in sample at 32'] = (.0967356,.0016985)
        targets['k then m in sample at 33'] = (.0914910,.001601)
        targets['k then m in sample at 34'] = (.0947637,.0015799)
        targets['k then m in sample at 35'] = (.1018985,.001609)
        
        #targets['hazard of new child at 22'] = (0.0136671,0.0005373)
        targets['hazard of new child at 23'] = (.0188677,.0005959)
        targets['hazard of new child at 24'] = (.0277849,.0007042)
        targets['hazard of new child at 25'] = (.0389933,.0008353)
        targets['hazard of new child at 26'] = (.0519110,.000965)
        targets['hazard of new child at 27'] = (.0661566,.0011131)
        targets['hazard of new child at 28'] = (.0813480,.0012859)
        targets['hazard of new child at 29'] = (.0893035,.0013959)
        targets['hazard of new child at 30'] = (.1023455,.0016000)
        targets['hazard of new child at 31'] = (.1039228,.0017302)
        targets['hazard of new child at 32'] = (.1029690,.0018698)
        targets['hazard of new child at 33'] = (.0970099,.0019398)
        targets['hazard of new child at 34'] = (.0884263,.0019678)
        targets['hazard of new child at 35'] = (.0799924,.0019688)
        
        
        
        
        
        
        targets['mean x share'] = (0.4,0.001)
        targets['divorced at 30 if one marriage'] = (.0686837,.0012484)
        targets['divorced if k then m and one marriage'] = (.1480121,.0024232*(1/2))
        targets['divorced if m then k and one marriage'] = (.0536316,.0004861*(1/2))
        targets['divorced with kids at 30']      = (.0267382,.0007793)
        targets['never married with kids at 30'] = (.0463395,.0008145)      
        targets['more than one mar at 40']       = (.1190139,.0012735)
        targets['in labor force at 30 if kids'] = (.739675,.0028066*(1/4))
        
        #targets['relative income at 30 if childless'] = (1.141977,.0144169*(1/5))
        #targets['men, relative income just married / single at 30'] = (1.197164,.0206279*(1/5))
        #targets['men, relative income with kids / no kids at 30'] = (0.9478894,.011772*(1/5))
        
        
        '''
        targets['men, relative income with kids / no kids at 24'] = (1.005178,.0293669)
        targets['men, relative income with kids / no kids at 25'] = (.9937299,.018525)
        targets['men, relative income with kids / no kids at 26'] = (.9629551,.0138367)
        targets['men, relative income with kids / no kids at 27'] = (.9670526,.0122425)
        targets['men, relative income with kids / no kids at 28'] = (.9550475,.0096132)
        targets['men, relative income with kids / no kids at 29'] = (.9499401,.0091448)
        targets['men, relative income with kids / no kids at 30'] = (.9641474,.008827)
        targets['men, relative income with kids / no kids at 31'] = (1.000548,.0096499)
        targets['men, relative income with kids / no kids at 32'] = (1.012243,.0093645)
        targets['men, relative income with kids / no kids at 33'] = (1.035536,.0105833)
        targets['men, relative income with kids / no kids at 34'] = (1.069198,.0112587)
        targets['men, relative income with kids / no kids at 35'] = (1.097902,.0121228)
        
        
        targets['men, relative income just married / single at 24'] = (1.158090,.0192292)
        targets['men, relative income just married / single at 25'] = (1.091558,.0177764)
        targets['men, relative income just married / single at 26'] = (1.108925,.0149643)
        targets['men, relative income just married / single at 27'] = (1.125526,.015914)
        targets['men, relative income just married / single at 28'] = (1.159013,.015304)
        targets['men, relative income just married / single at 29'] = (1.217950,.018949)
        targets['men, relative income just married / single at 30'] = (1.195528,.0175883)
        targets['men, relative income just married / single at 31'] = (1.208334,.0215594)
        targets['men, relative income just married / single at 31'] = (1.252952,.0225951)
        targets['men, relative income just married / single at 32'] = (1.252952,.0225951)
        targets['men, relative income just married / single at 33'] = (1.285158,.029411)
        targets['men, relative income just married / single at 34'] = (1.260427,.0308268)
        targets['men, relative income just married / single at 35'] = (1.222406,.0344223)
        '''
        
        
        targets['divorced by years after marriage if kids first, 1'] = (.0113142,.0020309)
        targets['divorced by years after marriage if kids first, 2'] = (.0216539,.0027502)
        targets['divorced by years after marriage if kids first, 3'] = (.0492118,.0041846)
        targets['divorced by years after marriage if kids first, 4'] = (.0893572,.0053861)
        targets['divorced by years after marriage if kids first, 5'] = (.1054101,.0058989)
        targets['divorced by years after marriage if kids first, 6'] = (.1386944,.0068164)
        targets['divorced by years after marriage if kids first, 7'] = (.1506841,.0072127)
        targets['divorced by years after marriage if kids first, 8'] = (.1866493,.0080136)
        targets['divorced by years after marriage if kids first, 9'] = (.1786823,.0078907)
        targets['divorced by years after marriage if kids first, 10'] = (.20118,.0084048)
        
        
        targets['divorced by years after marriage if marriage first, 1'] = (.000818,.0005133)
        targets['divorced by years after marriage if marriage first, 2'] = (.0050849,.000689)
        targets['divorced by years after marriage if marriage first, 3'] = (.0102901,.0007538)
        targets['divorced by years after marriage if marriage first, 4'] = (.0165388,.000839)
        targets['divorced by years after marriage if marriage first, 5'] = (.0201483,.000850)
        targets['divorced by years after marriage if marriage first, 6'] = (.0268928,.0009448)
        targets['divorced by years after marriage if marriage first, 7'] = (.0320433,.0010079)
        targets['divorced by years after marriage if marriage first, 8'] = (.0401776,.001116)
        targets['divorced by years after marriage if marriage first, 9'] = (.0479392,.0012191)
        targets['divorced by years after marriage if marriage first, 10'] = (.0602147,.0013719)
        
        
        
        
        
        
    elif mode=='low education':
    
        
        
        targets['ever kids by years after marriage, 1'] = (.6371565,.0030758)
        targets['ever kids by years after marriage, 2'] = (.694856,.0028766)
        targets['ever kids by years after marriage, 3'] = (.7534413,.0027062)
        targets['ever kids by years after marriage, 4'] = (.7929944,.0026026)
        targets['ever kids by years after marriage, 5'] = (.8249145,.0024629)
        targets['ever kids by years after marriage, 6'] = (.8519912,.00235419)
        #targets['ever kids by years after marriage, 7'] = (.8737004,.0022449)
        #targets['ever kids by years after marriage, 8'] = (.881806,.0022094)
        #targets['ever kids by years after marriage, 9'] = (.8990071,.0021148)
        #targets['ever kids by years after marriage, 10']= (.9025669,.0020562)
        
        
        
 
        
        targets['divorced by years after marriage, 1'] =  (.0231847,.0009513)
        targets['divorced by years after marriage, 2'] =  (.0601223,.0014428)
        targets['divorced by years after marriage, 3'] =  (.0943086,.0017544)        
        targets['divorced by years after marriage, 4'] =  (.1245948,.001997)
        targets['divorced by years after marriage, 5'] =  (.153307,.0021667)
        targets['divorced by years after marriage, 6'] =  (.1724024,.0022962)
        targets['divorced by years after marriage, 7'] =  (.1908162,.002408)
        targets['divorced by years after marriage, 8'] =  (.2068432,.002495*(1/2))
        targets['divorced by years after marriage, 9'] =  (.216311,.0025916)
        targets['divorced by years after marriage, 10'] = (.2174107,.0025545*(1/2))
        
        
        
        
  
        
        targets['hazard of marriage at 22'] = (.0550399,.0010074)
        targets['hazard of marriage at 23'] = (.0573569,.001092)
        targets['hazard of marriage at 24'] = (.058633,.0011528)
        targets['hazard of marriage at 25'] = (.0594683,.0011878)
        targets['hazard of marriage at 26'] = (.059702,.0012534)
        targets['hazard of marriage at 27'] = (.0616718,.0013152)
        targets['hazard of marriage at 28'] = (.0620773,.0013581)
        targets['hazard of marriage at 29'] = (.0621857,.0014114)
        targets['hazard of marriage at 30'] = (.0596791,.001373)
        targets['hazard of marriage at 31'] = (.0600205,.0014634)
        targets['hazard of marriage at 32'] = (.0585699,.0014619)
        targets['hazard of marriage at 33'] = (.055078,.0014608)
        targets['hazard of marriage at 34'] = (.0535931,.0014591)
        targets['hazard of marriage at 35'] = (.0505956,.0014132)
        
        
        
        targets['hazard of marriage & having a child at 22'] = (.0250424,.0006871)
        targets['hazard of marriage & having a child at 23'] = (.0261222,.0007443)
        targets['hazard of marriage & having a child at 24'] = (.0286651,.0008126)
        targets['hazard of marriage & having a child at 25'] = (.0328221,.0008876)
        targets['hazard of marriage & having a child at 26'] = (.0336944,.0009461)
        targets['hazard of marriage & having a child at 27'] = (.0360572,.0010091)
        targets['hazard of marriage & having a child at 28'] = (.0381652,.0010673)
        targets['hazard of marriage & having a child at 29'] = (.0399589,.0011331)
        targets['hazard of marriage & having a child at 30'] = (.0389913,.0011104)
        targets['hazard of marriage & having a child at 31'] = (.039547, .0011884)
        targets['hazard of marriage & having a child at 32'] = (.0393054,.0011971)
        targets['hazard of marriage & having a child at 33'] = (.0387907,.0012253)
        targets['hazard of marriage & having a child at 34'] = (.0382023,.0012309)
        targets['hazard of marriage & having a child at 35'] = (.0363907,.0011973)
        
        
        
        targets['hazard of new child at 22'] = (.0706884,.0012251)
        targets['hazard of new child at 23'] = (.0715649,.001351)
        targets['hazard of new child at 24'] = (.0673874,.0014081)
        targets['hazard of new child at 25'] = (.064768,.0014436)
        targets['hazard of new child at 26'] = (.0631538,.0015125)
        targets['hazard of new child at 27'] = (.061789,.0015726)
        targets['hazard of new child at 28'] = (.0617785,.0016493)
        targets['hazard of new child at 29'] = (.0631569,.001762)
        targets['hazard of new child at 30'] = (.0551772,.0016781)
        targets['hazard of new child at 31'] = (.0521211,.0016806)
        targets['hazard of new child at 32'] = (.0496799,.0017484)
        targets['hazard of new child at 33'] = (.0425128,.0016656)
        targets['hazard of new child at 34'] = (.0393337,.0016417)
        targets['hazard of new child at 35'] = (.0389228,.001653)
        
        
        
        
        targets['mean x share'] = (0.4,0.001)
        
        targets['divorced at 30 if one marriage'] = (.1427936,.0019867)
        
        
        
        targets['k then m by age, b0'] = (.3759993,.0015512)
        targets['k then m by age, b1'] = (-.0113458,.000246)
        targets['k then m by age, b2'] = (-.0002164,.0000395)
        
        
        targets['divorced if k then m and one marriage'] = (.1727448,.0014886*(1/4))
        targets['divorced if m then k and one marriage'] = (.1394705,.000983*(1/4))
        
        
        targets['divorced with kids at 30']      = (.0899426,.001527 )
        targets['never married with kids at 30'] = (.2156013,.0017011)      
        targets['more than one mar at 40']       = (.1770624,.0014789)
        
        targets['in labor force at 30 if kids'] = (.5466037,.0033879*(1/4))
        
        
        
    else:
        raise Exception('this mode for targets is not found')
    
    return targets

    
