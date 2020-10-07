#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 07:59:29 2020

@author: egorkozlov
"""

from itertools import chain

def target_values(mode='low education'):
        
    
    
    def unroll_key(key):
        assert type(key) is tuple
        if len(key) == 1:
            name, = key
            return [(name,1.0)]
        else:
            try:
                name, (start,finish), factor = key
            except ValueError:
                name, (start,finish) = key
                factor = 1.0
                
            return [((name + str(num)),factor) for num in range(start,finish+1)]
        
        
    
    all_t_pick = all_targets(mode)
    
    
    keys = [
            ('ever kids by years after marriage, ',(4,8),1),
            ('k then m in population at ',(23,35),1),
            ('m then k in population at ',(23,35),1),
            ('k then m in sample at ',(23,35),1),
            ('divorced by years after marriage if kids first, ',(1,10),1),
            ('divorced by years after marriage if marriage first, ',(1,10),1),
            ('mean x share',),
            ('more than one mar at 40',),
            ('in labor force at 30 if kids',),
            ('ever married at ',(23,35),1),
            #('divorced and kids in population at ',(28,35),1),
            ('divorced and no kids in population at ',(28,35),1),
            #('remarriage chance if kids, 26-30',),  
            #('remarriage chance if kids, 31-35',),
            #('remarriage chance if no kids, 26-30',),  
            #('remarriage chance if no kids, 31-35',),
            ('unplanned pregnancies aborted',),
            ('abortion 30s over 20s',),
            ('sorting overall',)            
           ]
    
    
    keys_all = chain(*[unroll_key(key) for key in keys])
    targets = {key: (all_t_pick[key][0], factor*all_t_pick[key][1]) for (key, factor) in keys_all}
        
    return targets



   
def all_targets(pick=None):
    
    all_targets = dict()
    targets = dict()
    
    targets['ever kids by years after marriage, 1'] = (.1901175,.0018521)
    targets['ever kids by years after marriage, 2'] = (.3416996,.0022054)
    targets['ever kids by years after marriage, 3'] = (.4902692,.0023415)
    targets['ever kids by years after marriage, 4'] = (.6114712,.0023096)
    targets['ever kids by years after marriage, 5'] = (.7039068,.0021802)
    targets['ever kids by years after marriage, 6'] = (.7730147,.0020489)
    targets['ever kids by years after marriage, 7'] = (.8182496,.0019268)
    targets['ever kids by years after marriage, 8'] = (.8568491,.00179)
    targets['ever kids by years after marriage, 9'] = (.8809587,.001696)
    targets['ever kids by years after marriage, 10']= (.8978599,.0016307)
    
    
    
    targets['divorced by years after marriage, 1']  =  (.0084912,.0004315)
    targets['divorced by years after marriage, 2']  =  (.0199536,.0006448)
    targets['divorced by years after marriage, 3']  =  (.0345649,.0008425)        
    targets['divorced by years after marriage, 4']  =  (.0536216,.001042)
    targets['divorced by years after marriage, 5']  =  (.0683989,.0011694)
    targets['divorced by years after marriage, 6']  =  (.0807483,.0012872)
    targets['divorced by years after marriage, 7']  =  (.0879516,.0013595)
    targets['divorced by years after marriage, 8']  =  (.1011097,.0014718)
    targets['divorced by years after marriage, 9']  =  (.1027806,.001517)
    targets['divorced by years after marriage, 10'] =  (.116314,.0016379)
    
    
    
    
    targets['hazard of marriage at 22'] = (.0506949,.001038)
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
    
    targets['k then m in population at 22'] = (.0031368,.0003093)
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
    
    

    
    targets['m then k in population at 22'] = (.0058368,.0004213)
    targets['m then k in population at 23'] = (.0150030,.0005547)
    targets['m then k in population at 24'] = (.0279987,.0007066)
    targets['m then k in population at 25'] = (.0480545,.0008848)
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
    
    
    targets['k then m in sample at 22'] = (.3495573,.0264092)
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
    
    targets['hazard of new child at 22'] = (0.0136671,0.0005373)
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

    targets['ever married at 22'] = (.0648409,.0013704)
    targets['ever married at 23'] = (.1092023,.0014232)
    targets['ever married at 24'] = (.1651963,.0015907)
    targets['ever married at 25'] = (.2321176,.0017455)
    targets['ever married at 26'] = (.3186594,.0019039)
    targets['ever married at 27'] = (.3979381,.0019605)
    targets['ever married at 28'] = (.4797466,.0019814)
    targets['ever married at 29'] = (.5492900,.0019650)
    targets['ever married at 30'] = (.6053985,.0018938)
    targets['ever married at 31'] = (.6557630,.0018604)
    targets['ever married at 32'] = (.7029164,.0017864)
    targets['ever married at 33'] = (.7277533,.0017525)
    targets['ever married at 34'] = (.7562652,.0016944)
    targets['ever married at 35'] = (.7776588,.0016421)
    
    
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
    
    
    targets['mean x share'] = (0.4,0.001)
    targets['divorced if k then m and one marriage'] = (.1480121,.0024232*(1/4))
    targets['divorced if m then k and one marriage'] = (.0536316,.0004861*(1/4))
    targets['more than one mar at 40']       = (.1190139,.0012735*(1/4))
    targets['in labor force at 30 if kids'] = (.739675,.0028066*(1/4))
    
    
    
    targets['divorced and kids in population at 23'] = (.0011839, .0001569)
    targets['divorced and kids in population at 24'] = (.0010305, .0001374)
    targets['divorced and kids in population at 25'] = (.0026032, .0002107)
    targets['divorced and kids in population at 26'] = (.0043844, .00027)
    targets['divorced and kids in population at 27'] = (.0061497, .0003131)
    targets['divorced and kids in population at 28'] = (.0094241, .0003832)
    targets['divorced and kids in population at 29'] = (.0123923, .0004369)
    targets['divorced and kids in population at 30'] = (.0161872, .000489)
    targets['divorced and kids in population at 31'] = (.0208707, .0005597)
    targets['divorced and kids in population at 32'] = (.0265332, .0006283)
    targets['divorced and kids in population at 33'] = (.0304566, .0006766)
    targets['divorced and kids in population at 34'] = (.0378053, .0007527)
    targets['divorced and kids in population at 35'] = (.0428613, .0007999)
    
    targets['divorced and no kids in population at 23'] = (.0032072, .000258)
    targets['divorced and no kids in population at 24'] = (.0048785, .0002984)
    targets['divorced and no kids in population at 25'] = (.0077115, .0003617)
    targets['divorced and no kids in population at 26'] = (.0106611, .0004196)
    targets['divorced and no kids in population at 27'] = (.0133708, .00046)
    targets['divorced and no kids in population at 28'] = (.0169658, .0005122)
    targets['divorced and no kids in population at 29'] = (.0207294, .0005627)
    targets['divorced and no kids in population at 30'] = (.0255615, .0006115)
    targets['divorced and no kids in population at 31'] = (.0250619, .0006121)
    targets['divorced and no kids in population at 32'] = (.0286549, .0006522)
    targets['divorced and no kids in population at 33'] = (.0268823, .0006368)
    targets['divorced and no kids in population at 34'] = (.0316118, .0006905)
    targets['divorced and no kids in population at 35'] = (.0320216, .0006953)
    
    targets['never married and kids in population at 23'] = (.019899, .0006373)
    targets['never married and kids in population at 24'] = (.0270769, .0006952)
    targets['never married and kids in population at 25'] = (.0304564, .0007104)
    targets['never married and kids in population at 26'] = (.0349298, .0007502)
    targets['never married and kids in population at 27'] = (.0393066, .0007783)
    targets['never married and kids in population at 28'] = (.0428057, .0008028)
    targets['never married and kids in population at 29'] = (.0439285, .0008094)
    targets['never married and kids in population at 30'] = (.0463395, .0008145)
    targets['never married and kids in population at 31'] = (.0452341, .0008137)
    targets['never married and kids in population at 32'] = (.04613, .00082)
    targets['never married and kids in population at 33'] = (.0487543, .0008479)
    targets['never married and kids in population at 34'] = (.0464966, .000831)
    targets['never married and kids in population at 35'] = (.043347, .0008042)
    
    targets['never married and no kids in population at 23'] = (.8708987, .0015301)
    targets['never married and no kids in population at 24'] = (.8077268, .001688)
    targets['never married and no kids in population at 25'] = (.737426, .0018192)
    targets['never married and no kids in population at 26'] = (.6464109, .0019534)
    targets['never married and no kids in population at 27'] = (.5627553, .0019869)
    targets['never married and no kids in population at 28'] = (.4774477, .001981)
    targets['never married and no kids in population at 29'] = (.4067815, .00194)
    targets['never married and no kids in population at 30'] = (.348262, .0018459)
    targets['never married and no kids in population at 31'] = (.2990029, .0017926)
    targets['never married and no kids in population at 32'] = (.2509536, .0016949)
    targets['never married and no kids in population at 33'] = (.2234924, .0016402)
    targets['never married and no kids in population at 34'] = (.1972382, .0015704)
    targets['never married and no kids in population at 35'] = (.1789942, .0015138)
    
    
    targets['remarriage chance if no kids, 26-30'] = (.1068772,.0040101)
    targets['remarriage chance if no kids, 31-35'] = (.0971051,.0030155 )
    targets['remarriage chance if kids, 26-30']    = (.1357563,.0059429)
    targets['remarriage chance if kids, 31-35']    = (.10925,  .0030496)
    
    
    targets['divorced ratio above over below at 23'] = (.8509691, .1271747)
    targets['divorced ratio above over below at 24'] = (1.038425, .1212888)
    targets['divorced ratio above over below at 25'] = (.7918638, .0644932)
    targets['divorced ratio above over below at 26'] = (.7877995, .0524962)
    targets['divorced ratio above over below at 27'] = (.7188274, .04182)
    targets['divorced ratio above over below at 28'] = (.6949408, .0343558)
    targets['divorced ratio above over below at 29'] = (.7022896, .0315099)
    targets['divorced ratio above over below at 30'] = (.6495594, .0257173)
    targets['divorced ratio above over below at 31'] = (.7433905, .0275904)
    targets['divorced ratio above over below at 32'] = (.8449129, .0289524)
    targets['divorced ratio above over below at 33'] = (.7125299, .0247043)
    targets['divorced ratio above over below at 34'] = (.7796795, .024101)
    targets['divorced ratio above over below at 35'] = (.7753022, .0235935)
    
    
    targets['above median among divorced mothers at 23'] = (.6155104, .091935)
    targets['above median among divorced mothers at 24'] = (.5464597, .0750517)
    targets['above median among divorced mothers at 25'] = (.4653033, .0482203)
    targets['above median among divorced mothers at 26'] = (.4200451, .0363862)
    targets['above median among divorced mothers at 27'] = (.4381934, .028839)
    targets['above median among divorced mothers at 28'] = (.4157751, .0230296)
    targets['above median among divorced mothers at 29'] = (.3484932, .0193085)
    targets['above median among divorced mothers at 30'] = (.3504218, .0165505)
    targets['above median among divorced mothers at 31'] = (.4178161, .0149386)
    targets['above median among divorced mothers at 32'] = (.4334456, .0135022)
    targets['above median among divorced mothers at 33'] = (.3986076, .0125915)
    targets['above median among divorced mothers at 34'] = (.4037833, .0113524)
    targets['above median among divorced mothers at 35'] = (.4348327, .0106641)
    
    
    targets['above median among never married mothers at 23'] = (.479423, .0199036)
    targets['above median among never married mothers at 24'] = (.3803031, .0155953)
    targets['above median among never married mothers at 25'] = (.3823566, .013396)
    targets['above median among never married mothers at 26'] = (.3689366, .0124172)
    targets['above median among never married mothers at 27'] = (.3326871, .0111305)
    targets['above median among never married mothers at 28'] = (.3921138, .011011)
    targets['above median among never married mothers at 29'] = (.3314707, .0103466)
    targets['above median among never married mothers at 30'] = (.3582643, .0101812)
    targets['above median among never married mothers at 31'] = (.3751338, .0104006)
    targets['above median among never married mothers at 32'] = (.3536149, .0102537)
    targets['above median among never married mothers at 33'] = (.3700371, .0102495)
    targets['above median among never married mothers at 34'] = (.3538409, .0103363)
    targets['above median among never married mothers at 35'] = (.4144861, .0110767)
    
    targets['abortion 30s over 20s'] = (25/65,0.001)
    targets['abortion ratio'] = (0.41*0.85*186,0.2)
        
    targets['unplanned pregnancies aborted'] = (0.4,0.001)
    targets['sorting at 30'] = (0.6,0.0005)
    targets['sorting overall'] = (0.6140971,0.0014194) 
    
    all_targets['high education'] = targets.copy()
    
    
    targets = dict()
    
    
    targets['ever kids by years after marriage, 1'] = (.6371565, .0030758)
    targets['ever kids by years after marriage, 2'] = (.694856, .0028766)
    targets['ever kids by years after marriage, 3'] = (.7534413, .0027062)
    targets['ever kids by years after marriage, 4'] = (.7929944, .0026026)
    targets['ever kids by years after marriage, 5'] = (.8249145, .0024629)
    targets['ever kids by years after marriage, 6'] = (.8519912, .0023541)
    targets['ever kids by years after marriage, 7'] = (.8737004, .0022449)
    targets['ever kids by years after marriage, 8'] = (.881806, .0022094)
    targets['ever kids by years after marriage, 9'] = (.8990071, .0021148)
    targets['ever kids by years after marriage, 10'] = (.9025669, .0020562)
    
    
    targets['divorced by years after marriage, 1'] = (.0231847, .0009513)
    targets['divorced by years after marriage, 2'] = (.0601223, .0014428)
    targets['divorced by years after marriage, 3'] = (.0943086, .0017544)
    targets['divorced by years after marriage, 4'] = (.1245948, .001997)
    targets['divorced by years after marriage, 5'] = (.153307, .0021667)
    targets['divorced by years after marriage, 6'] = (.1724024, .0022962)
    targets['divorced by years after marriage, 7'] = (.1908162, .002408)
    targets['divorced by years after marriage, 8'] = (.2068432, .002495)
    targets['divorced by years after marriage, 9'] = (.216311, .0025916)
    targets['divorced by years after marriage, 10'] = (.2174107, .0025545)
    
    
    
    targets['hazard of marriage at 22'] = (.0459992, .0009848)
    targets['hazard of marriage at 23'] = (.0480085, .0010465)
    targets['hazard of marriage at 24'] = (.0507163, .0011001)
    targets['hazard of marriage at 25'] = (.0505483, .0011527)
    targets['hazard of marriage at 26'] = (.0513198, .0011989)
    targets['hazard of marriage at 27'] = (.0533364, .0012555)
    targets['hazard of marriage at 28'] = (.0534521, .0013044)
    targets['hazard of marriage at 29'] = (.0518105, .0012773)
    targets['hazard of marriage at 30'] = (.0526195, .0013629)
    targets['hazard of marriage at 31'] = (.0509349, .001358)
    targets['hazard of marriage at 32'] = (.0477095, .0013521)
    targets['hazard of marriage at 33'] = (.0461142, .0013488)
    targets['hazard of marriage at 34'] = (.0437701, .0013087)
    targets['hazard of marriage at 35'] = (.0398101, .0012805)
    
    
    
    targets['k then m in population at 22'] = (.0471337, .0008786)
    targets['k then m in population at 23'] = (.0630275, .0010425)
    targets['k then m in population at 24'] = (.0765292, .001157)
    targets['k then m in population at 25'] = (.0891691, .0012334)
    targets['k then m in population at 26'] = (.0974226, .0013089)
    targets['k then m in population at 27'] = (.1065944, .0013648)
    targets['k then m in population at 28'] = (.1122013, .0013894)
    targets['k then m in population at 29'] = (.1141492, .0014134)
    targets['k then m in population at 30'] = (.1187029, .0013998)
    targets['k then m in population at 31'] = (.1203012, .0014483)
    targets['k then m in population at 32'] = (.1200893, .0014305)
    targets['k then m in population at 33'] = (.1245074, .001461)
    targets['k then m in population at 34'] = (.1270016, .0014689)
    targets['k then m in population at 35'] = (.1215744, .0014184)
    
    

    
    targets['m then k in population at 22'] = (.0569954, .0009611)
    targets['m then k in population at 23'] = (.0745946, .0011272)
    targets['m then k in population at 24'] = (.0946999, .0012743)
    targets['m then k in population at 25'] = (.113709, .001374)
    targets['m then k in population at 26'] = (.1358313, .0015123)
    targets['m then k in population at 27'] = (.1569375, .0016087)
    targets['m then k in population at 28'] = (.1726471, .0016638)
    targets['m then k in population at 29'] = (.1872816, .001734)
    targets['m then k in population at 30'] = (.2009962, .0017343)
    targets['m then k in population at 31'] = (.2210822, .0018475)
    targets['m then k in population at 32'] = (.2294688, .0018504)
    targets['m then k in population at 33'] = (.2389265, .001887)
    targets['m then k in population at 34'] = (.2457679, .0018993)
    targets['m then k in population at 35'] = (.249713, .0018788)
    
    
    targets['k then m in sample at 22'] = (.4526467, .0062011)
    targets['k then m in sample at 23'] = (.4579749, .0054619)
    targets['k then m in sample at 24'] = (.4469404, .0049857)
    targets['k then m in sample at 25'] = (.4395206, .0045602)
    targets['k then m in sample at 26'] = (.4176675, .0043115)
    targets['k then m in sample at 27'] = (.4044838, .0040408)
    targets['k then m in sample at 28'] = (.3938983, .003864)
    targets['k then m in sample at 29'] = (.3786911, .0037764)
    targets['k then m in sample at 30'] = (.3712956, .0035625)
    targets['k then m in sample at 31'] = (.3523933, .003516)
    targets['k then m in sample at 32'] = (.3435461, .0034191)
    targets['k then m in sample at 33'] = (.3425861, .0033878)
    targets['k then m in sample at 34'] = (.3406974, .0033269)
    targets['k then m in sample at 35'] = (.3274401, .003233)
    
    
    
    targets['hazard of new child at 22'] = (.0723763, .0014004)
    targets['hazard of new child at 23'] = (.0686844, .0014701)
    targets['hazard of new child at 24'] = (.0660664, .0015141)
    targets['hazard of new child at 25'] = (.0656788, .0016069)
    targets['hazard of new child at 26'] = (.0640507, .0016735)
    targets['hazard of new child at 27'] = (.0641009, .0017645)
    targets['hazard of new child at 28'] = (.0665277, .0019017)
    targets['hazard of new child at 29'] = (.0588031, .0018288)
    targets['hazard of new child at 30'] = (.0558206, .0018416)
    targets['hazard of new child at 31'] = (.0531504, .0019242)
    targets['hazard of new child at 32'] = (.0456938, .0018372)
    targets['hazard of new child at 33'] = (.0415602, .0017985)
    targets['hazard of new child at 34'] = (.0408571, .0018053)
    targets['hazard of new child at 35'] = (.0343058, .001669)
    
    
    
    targets['ever married at 21'] = (.1198675, .0012708)
    targets['ever married at 22'] = (.1725177, .0015664)
    targets['ever married at 23'] = (.2182572, .0017721)
    targets['ever married at 24'] = (.2691566, .0019303)
    targets['ever married at 25'] = (.3205891, .0020199)
    targets['ever married at 26'] = (.369406, .0021304)
    targets['ever married at 27'] = (.4139921, .0021783)
    targets['ever married at 28'] = (.4536668, .0021917)
    targets['ever married at 29'] = (.4897026, .0022219)
    targets['ever married at 30'] = (.5310084, .0021597)
    targets['ever married at 31'] = (.5669783, .002206)
    targets['ever married at 32'] = (.5937097, .0021613)
    targets['ever married at 33'] = (.6194696, .0021484)
    targets['ever married at 34'] = (.6386038, .0021192)
    targets['ever married at 35'] = (.6639089, .0020503)
    
    
    
    targets['divorced by years after marriage if kids first, 1'] = (.0201514, .0017159)
    targets['divorced by years after marriage if kids first, 2'] = (.0471255, .0025286)
    targets['divorced by years after marriage if kids first, 3'] = (.0813251, .0032273)
    targets['divorced by years after marriage if kids first, 4'] = (.1166151, .0038275)
    targets['divorced by years after marriage if kids first, 5'] = (.1376778, .0041438)
    targets['divorced by years after marriage if kids first, 6'] = (.1580319, .0044724)
    targets['divorced by years after marriage if kids first, 7'] = (.1836251, .0047913)
    targets['divorced by years after marriage if kids first, 8'] = (.2043214, .0049916)
    targets['divorced by years after marriage if kids first, 9'] = (.2066351, .0051605)
    targets['divorced by years after marriage if kids first, 10'] = (.2115082, .0050606)
    
    
    targets['divorced by years after marriage if marriage first, 1'] = (.010328, .0022173)
    targets['divorced by years after marriage if marriage first, 2'] = (.0289688, .002371)
    targets['divorced by years after marriage if marriage first, 3'] = (.0505473, .0025391)
    targets['divorced by years after marriage if marriage first, 4'] = (.0645814, .0025787)
    targets['divorced by years after marriage if marriage first, 5'] = (.0918597, .0028409)
    targets['divorced by years after marriage if marriage first, 6'] = (.105973, .0029251)
    targets['divorced by years after marriage if marriage first, 7'] = (.113954, .0029394)
    targets['divorced by years after marriage if marriage first, 8'] = (.1333069, .0030814)
    targets['divorced by years after marriage if marriage first, 9'] = (.1407658, .0031131)
    targets['divorced by years after marriage if marriage first, 10'] = (.1443613, .0030424)
    
    
    targets['mean x share'] = (0.4,0.001)
    targets['divorced if k then m and one marriage'] = (.1727448,.0014886*(1/4))
    targets['divorced if m then k and one marriage'] = (.1394705,.000983*(1/4))
    targets['more than one mar at 40']       = (.1714891,.0015326*(1/4))
    targets['in labor force at 30 if kids'] = (.5466037,.0033879*(1/4))
    
    
    targets['divorced and kids in population at 23'] = (.0128405, .000483)
    targets['divorced and kids in population at 24'] = (.0182726, .0005829)
    targets['divorced and kids in population at 25'] = (.0256893, .0006847)
    targets['divorced and kids in population at 26'] = (.0317844, .0007743)
    targets['divorced and kids in population at 27'] = (.0382502, .0008482)
    targets['divorced and kids in population at 28'] = (.0425242, .0008883)
    targets['divorced and kids in population at 29'] = (.0536633, .0010016)
    targets['divorced and kids in population at 30'] = (.0566093, .0010001)
    targets['divorced and kids in population at 31'] = (.0690686, .0011289)
    targets['divorced and kids in population at 32'] = (.0735763, .0011489)
    targets['divorced and kids in population at 33'] = (.0803377, .0012028)
    targets['divorced and kids in population at 34'] = (.0875592, .0012469)
    targets['divorced and kids in population at 35'] = (.0941355, .0012675)
    
    
    targets['divorced and no kids in population at 23'] = (.011186, .0004512)
    targets['divorced and no kids in population at 24'] = (.0149702, .0005285)
    targets['divorced and no kids in population at 25'] = (.0189704, .0005904)
    targets['divorced and no kids in population at 26'] = (.022642, .0006566)
    targets['divorced and no kids in population at 27'] = (.0214305, .0006404)
    targets['divorced and no kids in population at 28'] = (.0240053, .0006738)
    targets['divorced and no kids in population at 29'] = (.0259961, .0007073)
    targets['divorced and no kids in population at 30'] = (.0317722, .0007591)
    targets['divorced and no kids in population at 31'] = (.0312788, .000775)
    targets['divorced and no kids in population at 32'] = (.0304979, .0007567)
    targets['divorced and no kids in population at 33'] = (.0327354, .0007874)
    targets['divorced and no kids in population at 34'] = (.0336342, .0007953)
    targets['divorced and no kids in population at 35'] = (.0399444, .00085)
    
    
    targets['never married and kids in population at 23'] = (.2555029, .0018711)
    targets['never married and kids in population at 24'] = (.2705469, .0019334)
    targets['never married and kids in population at 25'] = (.269658, .0019207)
    targets['never married and kids in population at 26'] = (.2706259, .0019611)
    targets['never married and kids in population at 27'] = (.2714562, .0019667)
    targets['never married and kids in population at 28'] = (.2683261, .0019506)
    targets['never married and kids in population at 29'] = (.2646187, .0019607)
    targets['never married and kids in population at 30'] = (.2368133, .0018398)
    targets['never married and kids in population at 31'] = (.241201, .0019046)
    targets['never married and kids in population at 32'] = (.2286718, .0018482)
    targets['never married and kids in population at 33'] = (.2173221, .001825)
    targets['never married and kids in population at 34'] = (.2097827, .0017961)
    targets['never married and kids in population at 35'] = (.1935926, .001715)
   
    
    targets['never married and no kids in population at 23'] = (.5262399, .0021421)
    targets['never married and no kids in population at 24'] = (.4602966, .0021692)
    targets['never married and no kids in population at 25'] = (.4097529, .0021285)
    targets['never married and no kids in population at 26'] = (.3599681, .0021187)
    targets['never married and no kids in population at 27'] = (.3145517, .0020535)
    targets['never married and no kids in population at 28'] = (.2780071, .0019723)
    targets['never married and no kids in population at 29'] = (.2456788, .0019134)
    targets['never married and no kids in population at 30'] = (.2321783, .0018273)
    targets['never married and no kids in population at 31'] = (.1918207, .0017529)
    targets['never married and no kids in population at 32'] = (.1776185, .0016819)
    targets['never married and no kids in population at 33'] = (.1632083, .0016353)
    targets['never married and no kids in population at 34'] = (.1516136, .0015821)
    targets['never married and no kids in population at 35'] = (.1424985, .0015173)
    
    
    targets['remarriage chance if no kids, 26-30'] = (.056675, .0027803)
    targets['remarriage chance if no kids, 31-35'] = (.0521429,.0023292)
    targets['remarriage chance if kids, 26-30']    = (.1068772,.0040101)
    targets['remarriage chance if kids, 31-35']    = (.0971051,.0030155 )
    
    
    
    targets['divorced ratio above over below at 23'] = (.801952, .0531031)
    targets['divorced ratio above over below at 24'] = (.838632, .0466088)
    targets['divorced ratio above over below at 25'] = (.9354082, .0452757)
    targets['divorced ratio above over below at 26'] = (.9460736, .0418577)
    targets['divorced ratio above over below at 27'] = (.9424763, .0396133)
    targets['divorced ratio above over below at 28'] = (.9870032, .0395609)
    targets['divorced ratio above over below at 29'] = (.94164, .0351559)
    targets['divorced ratio above over below at 30'] = (.9568118, .0331444)
    targets['divorced ratio above over below at 31'] = (1.012983, .0326977)
    targets['divorced ratio above over below at 32'] = (1.028341, .0331363)
    targets['divorced ratio above over below at 33'] = (1.033319, .0326461)
    targets['divorced ratio above over below at 34'] = (1.036322, .0306522)
    targets['divorced ratio above over below at 35'] = (1.078097, .0299674)
    
    
    targets['above median among divorced mothers at 23'] = (.5301601, .0253048)
    targets['above median among divorced mothers at 24'] = (.4818047, .0213643)
    targets['above median among divorced mothers at 25'] = (.4950256, .017925)
    targets['above median among divorced mothers at 26'] = (.4604572, .0165411)
    targets['above median among divorced mothers at 27'] = (.4860508, .0144643)
    targets['above median among divorced mothers at 28'] = (.5381645, .0136395)
    targets['above median among divorced mothers at 29'] = (.5232977, .0126617)
    targets['above median among divorced mothers at 30'] = (.5277151, .0116606)
    targets['above median among divorced mothers at 31'] = (.5202228, .0107895)
    targets['above median among divorced mothers at 32'] = (.5374461, .0103028)
    targets['above median among divorced mothers at 33'] = (.5158849, .0100985)
    targets['above median among divorced mothers at 34'] = (.508917, .0095279)
    targets['above median among divorced mothers at 35'] = (.5268524, .008963)
    
    
    targets['above median among never married mothers at 23'] = (.4899995, .006338)
    targets['above median among never married mothers at 24'] = (.4665557, .0061399)
    targets['above median among never married mothers at 25'] = (.4362643, .0060104)
    targets['above median among never married mothers at 26'] = (.4412364, .0060289)
    targets['above median among never married mothers at 27'] = (.4515305, .0061302)
    targets['above median among never married mothers at 28'] = (.4488592, .0061409)
    targets['above median among never married mothers at 29'] = (.4697128, .006266)
    targets['above median among never married mothers at 30'] = (.4841722, .0063757)
    targets['above median among never married mothers at 31'] = (.4392147, .0065432)
    targets['above median among never married mothers at 32'] = (.4400828, .0065726)
    targets['above median among never married mothers at 33'] = (.4438823, .0067927)
    targets['above median among never married mothers at 34'] = (.4294347, .0068816)
    targets['above median among never married mothers at 35'] = (.4517876, .0069976)
    
    
    targets['abortion ratio'] = (0.41*0.85*186,0.2)
    targets['abortion 30s over 20s'] = (25/65,0.001)
    targets['sorting at 30'] = (0.6,0.0005) #(0.630436,0.0022637)
    targets['sorting overall'] = (0.630436,0.0022637)
    
    targets['unplanned pregnancies aborted'] = (0.4,0.001)
    all_targets['low education'] = targets.copy()
        
    if pick is None:
        return all_targets
    else:
        return all_targets[pick]


    
