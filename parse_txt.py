#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:45:48 2020

@author: egorkozlov
"""


descs = ['above median among divorced at ',
         'above median among never married at ']   


xs = list()

xs.append(""" |  23   .5301601   .0253048 |
              |  24   .4818047   .0213643 |
              |  25   .4950256    .017925 |
              |---------------------------|
              |  26   .4604572   .0165411 |
              |  27   .4860508   .0144643 |
              |  28   .5381645   .0136395 |
              |  29   .5232977   .0126617 |
              |  30   .5277151   .0116606 |
              |---------------------------|
              |  31   .5202228   .0107895 |
              |  32   .5374461   .0103028 |
              |  33   .5158849   .0100985 |
              |  34    .508917   .0095279 |
              |  35   .5268524    .008963 |""")


xs.append(""" |  23   .4899995    .006338 |
              |  24   .4665557   .0061399 |
              |  25   .4362643   .0060104 |
              |---------------------------|
              |  26   .4412364   .0060289 |
              |  27   .4515305   .0061302 |
              |  28   .4488592   .0061409 |
              |  29   .4697128    .006266 |
              |  30   .4841722   .0063757 |
              |---------------------------|
              |  31   .4392147   .0065432 |
              |  32   .4400828   .0065726 |
              |  33   .4438823   .0067927 |
              |  34   .4294347   .0068816 |
              |  35   .4517876   .0069976 |""")

      
 
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
        