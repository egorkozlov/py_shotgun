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
        
        
        
        '''
        targets['hazard of marriage & having a child at 22'] = (.0006303,.0001446)
        targets['hazard of marriage & having a child at 23'] = (.0015481,.0001853)
        targets['hazard of marriage & having a child at 24'] = (.002005,.0002028)
        targets['hazard of marriage & having a child at 25'] = (.0028903,.0002427)
        targets['hazard of marriage & having a child at 26'] = (.0032498,.0002653)
        targets['hazard of marriage & having a child at 27'] = (.0037572,.000297)
        targets['hazard of marriage & having a child at 28'] = (.004184,.0003306)
        targets['hazard of marriage & having a child at 29'] = (.0055574,.0004085)
        targets['hazard of marriage & having a child at 30'] = (.0048407,.0004016)
        targets['hazard of marriage & having a child at 31'] = (.005135,.0004482)
        targets['hazard of marriage & having a child at 32'] = (.0051018,.0004756)
        targets['hazard of marriage & having a child at 33'] = (.0072451,.0006166)
        targets['hazard of marriage & having a child at 34'] = (.0062195,.0005979)
        targets['hazard of marriage & having a child at 35'] = (.005723,.0006062)
        '''

        #targets['k then m in population at 22'] = (.0031368,.0003093)
        targets['k then m in population at 23'] = (.0054749,.0003343)
        targets['k then m in population at 24'] = (.0079265,.0003768)
        targets['k then m in population at 25'] = (.0123658,.0004524)
        targets['k then m in population at 26'] = (.0173624,.0005278)
        targets['k then m in population at 27'] = (.0208649,.0005656)
        targets['k then m in population at 28'] = (.0273335,.0006384)
        targets['k then m in population at 29'] = (.0311646,.0006774)
        targets['k then m in population at 30'] = (.0357163,.0007086)
        targets['k then m in population at 31'] = (.0374948,.0007337)
        targets['k then m in population at 32'] = (.0409537,.0007628)
        targets['k then m in population at 33'] = (.0421048,.0007787)
        targets['k then m in population at 34'] = (.0467038,.0008199)
        targets['k then m in population at 35'] = (.0516670,.0008606)
        
        #targets['m then k in population at 22'] = (.0058368,.0004213)
        targets['m then k in population at 23'] = (.0147950,.000547)
        targets['m then k in population at 24'] = (.0275320,.0006952)
        targets['m then k in population at 25'] = (.0470913,.0008672)
        targets['m then k in population at 26'] = (.0770536,.0010776)
        targets['m then k in population at 27'] = (.1178683,.0012759)
        targets['m then k in population at 28'] = (.1688230,.0014667)
        targets['m then k in population at 29'] = (.2248659,.0016276)
        targets['m then k in population at 30'] = (.2726108,.0017002)
        targets['m then k in population at 31'] = (.3349063,.0018227)
        targets['m then k in population at 32'] = (.3824028,.0018704)
        targets['m then k in population at 33'] = (.4181024,.0019126)
        targets['m then k in population at 34'] = (.4461409,.0019316)
        targets['m then k in population at 35'] = (.4553771,.0019361)
        
        
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
        
        
        
        
        #targets['hazard of new child at 22'] = (0.0069326,0.0004666)
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
