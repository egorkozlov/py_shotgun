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
        
        
        targets['divorced if k then m and one marriage'] = (.1480121,.0024232*(1/4))
        targets['divorced if m then k and one marriage'] = (.0536316,.0004861*(1/4))
        
        
        targets['divorced with kids at 30']      = (.025285,.0007414)
        targets['never married with kids at 30'] = (.0447814,.0007897)      
        targets['more than one mar at 40']       = (.1221077,.0012642)
        
        targets['in labor force at 30 if kids'] = (.739675,.0028066*(1/4))
        

        
    
        
    elif mode=='low education':
    
        targets['ever kids by years after marriage, 1'] = (.1888399,.0018637)
        targets['ever kids by years after marriage, 2'] = (.3420612,.0022199)
        targets['ever kids by years after marriage, 3'] = (.4909945,.0023519)
        targets['ever kids by years after marriage, 4'] = (.6122718,.0023161)
        targets['ever kids by years after marriage, 5'] = (.7052276,.0021829)
        targets['ever kids by years after marriage, 6'] = (.7739091,.0020495)
        #targets['ever kids by years after marriage, 7'] = (.8188877,.0019268)
        #targets['ever kids by years after marriage, 8'] = (.8575637,.0017885)
        #targets['ever kids by years after marriage, 9'] = (.8813007,.0016954)
        #targets['ever kids by years after marriage, 10']= (.8981363,.0016297)
        
        
        targets['divorced by years after marriage, 1'] =  (.0074631,.0004086)
        targets['divorced by years after marriage, 2'] =  (.0194861,.0006415)
        targets['divorced by years after marriage, 3'] =  (.0340267,.0008401)        
        targets['divorced by years after marriage, 4'] =  (.0529621,.0010397)
        targets['divorced by years after marriage, 5'] =  (.0680016,.0011695)
        targets['divorced by years after marriage, 6'] =  (.0802289,.0012859)
        targets['divorced by years after marriage, 7'] =  (.087411,.001358)
        targets['divorced by years after marriage, 8'] =  (.1002189,.0014686*(1/2))
        targets['divorced by years after marriage, 9'] =  (.1018378,.0015131)
        targets['divorced by years after marriage, 10'] = (.115365,.001633*(1/2))
        
        
        targets['divorced at 30 if one marriage'] = (0.0619966, 0.0005809)
        
        targets['mean x share'] = (0.4,0.001)
        
        
        targets['k then m by age, b0'] = (0.1199747,0.0009933)
        targets['k then m by age, b1'] = (-0.0100934,0.000355)
        targets['k then m by age, b2'] = (0.0006642,0.0000398)
        
        
        targets['divorced if k then m and one marriage'] = (0.1474027, (1/4)*0.00243)
        targets['divorced if m then k and one marriage'] = (0.0531433, (1/4)*0.0004851)
        
        
        targets['divorced with kids at 30']      = (0.0251494,0.0007414)
        targets['never married with kids at 30'] = (0.0443802,0.0007902)      
        targets['more than one mar at 40']       = (0.1225949,0.0012733)
        
        targets['in labor force at 30 if kids'] = (0.7394259,(1/4)*0.0035188)
        targets['in labor force at 30 if kids ratio'] = (0.8641976,0.0085405)
        
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
        
        
        
        
    else:
        raise Exception('this mode for targets is not found')
    
    return targets
