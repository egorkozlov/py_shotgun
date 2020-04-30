#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

        
def target_values(mode='high education'):
    targets = dict()
        
        
    if mode=='high education':
        
        targets['ever kids by years after marriage, 1'] = (.1901175,.0018521)
        targets['ever kids by years after marriage, 2'] = (.3416996,.0022054)
        targets['ever kids by years after marriage, 3'] = (.4902692,.0023415)
        targets['ever kids by years after marriage, 4'] = (.6114712,.0023096)
        targets['ever kids by years after marriage, 5'] = (.7039068,.0021802)
        targets['ever kids by years after marriage, 6'] = (.7730147,.0020489)
        #targets['ever kids by years after marriage, 7'] = (.8182496,.0019268)
        #targets['ever kids by years after marriage, 8'] = (.8568491,.00179)
        #targets['ever kids by years after marriage, 9'] = (.8809587,.001696)
        #targets['ever kids by years after marriage, 10']= (.8978599,.0016307)
        
        
        targets['divorced by years after marriage, 1'] =  (.0084912,.0004315)
        targets['divorced by years after marriage, 2'] =  (.0199536,.0006448)
        targets['divorced by years after marriage, 3'] =  (.0345649,.0008425)        
        targets['divorced by years after marriage, 4'] =  (.0536216,.001042)
        targets['divorced by years after marriage, 5'] =  (.0683989,.0011694)
        targets['divorced by years after marriage, 6'] =  (.0807483,.0012872)
        targets['divorced by years after marriage, 7'] =  (.0879516,.0013595)
        targets['divorced by years after marriage, 8'] =  (.1011097,.0014718*(1/2))
        targets['divorced by years after marriage, 9'] =  (.1027806,.001517)
        targets['divorced by years after marriage, 10'] = (.116314,.0016379*(1/2))
        
        
        
        
        
        targets['hazard of marriage at 22'] = (0.0414921,0.0011577)
        targets['hazard of marriage at 23'] = (0.0555908,0.0010842)
        targets['hazard of marriage at 24'] = (0.0680466,0.0011492)
        targets['hazard of marriage at 25'] = (0.0843476,0.0012718)
        targets['hazard of marriage at 26'] = (0.1013612,.0014339)
        targets['hazard of marriage at 27'] = (0.1140159,0.0015872)
        targets['hazard of marriage at 28'] = (0.1189503,0.0017143)
        targets['hazard of marriage at 29'] = (0.1275065,0.0019024)
        targets['hazard of marriage at 30'] = (0.1265993,0.0019926)
        targets['hazard of marriage at 31'] = (0.1189422,0.0020686)
        targets['hazard of marriage at 32'] = (0.1137438,0.0021385)
        targets['hazard of marriage at 33'] = (0.1141971,0.002287)
        targets['hazard of marriage at 34'] = (0.09789,0.0021931)
        targets['hazard of marriage at 35'] = (0.0909237,0.0021694)
        
        targets['hazard of marriage & having a child at 22'] = (0.0020043,0.0002575)
        targets['hazard of marriage & having a child at 23'] = (0.0037913,0.0002894)
        targets['hazard of marriage & having a child at 24'] = (0.0050758,0.0003216)
        targets['hazard of marriage & having a child at 25'] = (0.0078502,0.0003979)
        targets['hazard of marriage & having a child at 26'] = (0.0099153,0.0004596)
        targets['hazard of marriage & having a child at 27'] = (0.0128215,0.0005419)
        targets['hazard of marriage & having a child at 28'] = (0.0150063,0.0006154)
        targets['hazard of marriage & having a child at 29'] = (0.0201393,0.0007585)
        targets['hazard of marriage & having a child at 30'] = (0.0222486,0.0008317)
        targets['hazard of marriage & having a child at 31'] = (0.0235701,0.0009146)
        targets['hazard of marriage & having a child at 32'] = (0.0272916,0.0010366)
        targets['hazard of marriage & having a child at 33'] = (0.0330752,0.0012157)
        targets['hazard of marriage & having a child at 34'] = (0.0321491,0.0012423)
        targets['hazard of marriage & having a child at 35'] = (0.0314278,0.0012609)
        
        targets['hazard of new child at 22'] = (0.0069326,0.0004666)
        targets['hazard of new child at 23'] = (0.0138225,0.0005367)
        targets['hazard of new child at 24'] = (0.0188708,0.0005915)
        targets['hazard of new child at 25'] = (0.0277794,0.0006976)
        targets['hazard of new child at 26'] = (0.0386756,0.0008233)
        targets['hazard of new child at 27'] = (0.0515395,0.0009507)
        targets['hazard of new child at 28'] = (0.0652943,0.0010927)
        targets['hazard of new child at 29'] = (0.0799002,0.0012581)
        targets['hazard of new child at 30'] = (0.0876844,0.0013633)
        targets['hazard of new child at 31'] = (0.1004221,0.0015608)
        targets['hazard of new child at 32'] = (0.1016521,0.0016817)
        targets['hazard of new child at 33'] = (0.1003621,0.0018118)
        targets['hazard of new child at 34'] = (0.0947235,0.0018794)
        targets['hazard of new child at 35'] = (0.0858545,0.0018973)
        
        
        
        targets['mean x share'] = (0.4,0.001)
        
        targets['divorced at 30 if one marriage'] = (.0650167,.0011904)
        
        
        
        targets['k then m by age, b0'] = (.1204838,.0009941)
        targets['k then m by age, b1'] = (-.0101522,.0003555)
        targets['k then m by age, b2'] = (.0006733,.0000399)
        
        
        targets['divorced if k then m and one marriage'] = (.1480121,.0024232*(1/3))
        targets['divorced if m then k and one marriage'] = (.0536316,.0004861*(1/3))
        
        
        targets['divorced with kids at 30']      = (.025285,.0007414)
        targets['never married with kids at 30'] = (.0447814,.0007897)      
        targets['more than one mar at 40']       = (.1221077,.0012642)
        
        targets['in labor force at 30 if kids'] = (.739675,.0028066*(1/4))
        

        
    
        
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
