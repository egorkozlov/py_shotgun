#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:45:48 2020

@author: egorkozlov
"""


descs = ['divorced ratio above over below at ',
         'divorced ratio above over below at ']   


xs = list()

xs.append(""" |  23   .8509691   .1271747 |
              |  24   1.038425   .1212888 |
              |  25   .7918638   .0644932 |
              |---------------------------|
              |  26   .7877995   .0524962 |
              |  27   .7188274     .04182 |
              |  28   .6949408   .0343558 |
              |  29   .7022896   .0315099 |
              |  30   .6495594   .0257173 |
              |---------------------------|
              |  31   .7433905   .0275904 |
              |  32   .8449129   .0289524 |
              |  33   .7125299   .0247043 |
              |  34   .7796795    .024101 |
              |  35   .7753022   .0235935 |""")


xs.append("""     |  23    .801952   .0531031 |
                  |  24    .838632   .0466088 |
                  |  25   .9354082   .0452757 |
                  |---------------------------|
                  |  26   .9460736   .0418577 |
                  |  27   .9424763   .0396133 |
                  |  28   .9870032   .0395609 |
                  |  29     .94164   .0351559 |
                  |  30   .9568118   .0331444 |
                  |---------------------------|
                  |  31   1.012983   .0326977 |
                  |  32   1.028341   .0331363 |
                  |  33   1.033319   .0326461 |
                  |  34   1.036322   .0306522 |
                  |  35   1.078097   .0299674 |
                  |---------------------------|
                  |  36   1.150801    .031961 |
                  |  37   1.030538   .0278805 |
                  |  38   1.098337   .0284611 |
                  |  39   1.109698   .0289851 |
                  |  40   1.172143   .0278424 |""")

      
 
def parse_line(line,begin='|',end='|',delim=' '):
    i = line.index(begin) + 1
    iend = line.index(begin,i)
    assert line[i] == delim   
    parsed = list()
    while i < iend:
        if line[i] == ' ':
            i+=1
            continue
        if line[i] == end:
            break
        j = line.index(' ',i)
        parsed.append(line[i:j])
        i = j
    
    return parsed
        
        

for xin, desc in zip(xs,descs):
    for line in xin.splitlines():
        try:
            parsed = parse_line(line)
        except (ValueError, AssertionError):
            continue
        print("targets['{}{}'] = ({}, {})".format(desc,*parsed))        
        