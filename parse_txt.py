#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:45:48 2020

@author: egorkozlov
"""


desc = 'divorced by years after marriage if kids first, '     

xin =  """       |       1    .010328   .0022173 |
  |       2   .0289688    .002371 |
  |       3   .0505473   .0025391 |
  |       4   .0645814   .0025787 |
  |       5   .0918597   .0028409 |
  |-------------------------------|
  |       6    .105973   .0029251 |
  |       7    .113954   .0029394 |
  |       8   .1333069   .0030814 |
  |       9   .1407658   .0031131 |
  |      10   .1443613   .0030424 |



       """
      
 
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
        
        


for line in xin.splitlines():
    try:
        parsed = parse_line(line)
    except (ValueError, AssertionError):
        continue
    print("targets['{}{}'] = ({}, {})".format(desc,*parsed))        
    