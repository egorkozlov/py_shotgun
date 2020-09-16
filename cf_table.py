#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:52:31 2020

@author: egorkozlov
"""

# this parses moment files to table

from tiktak import filer

names = [('BL','baseline'),
         ('NS','no social stigma'),
         ('FC','costless abortion'),
         ('NR','no remar penalty'),
         #('SD','no skills depreciation'),
         ('FD','no divorce costs'),
         ('ND','infinite divorce costs'),         
         #('PG','no pay gap')
         ]

educ = 'col'

def file_name(nm): return '{} {}.pkl'.format(educ,nm)




entries = [('divorced in 10 years, kids-first','divorced by years after marriage if kids first, 10'),
           ('divorced in 10 years, marriage-first','divorced by years after marriage if marriage first, 10'),
           ('kids-first at 30','k then m in sample at 30'),
           ('unplanned pregnancies aborted','unplanned pregnancies aborted'),
           ('single mothers at 35','single mothers among mothers at 35')]

table_cols = list()
for ename, nm in names:
    fname = file_name(nm)
    mom = filer(fname,0,0)
    table_col = [r'\textbf{' + ename + r'}'] + [mom[e] for _, e in entries]
    table_cols.append(table_col)
    
    
table_left = [r'\textbf{Experiment}'] + [r'\textit{' + e + r'}' for e, _ in entries]

# build table rows

table_rows = list()
for irow in range(len(entries)+1):
    table_row = table_left[irow]
    for icol in range(len(names)):
        try:
            table_row += ' & ' + '{:02.1f}'.format(100*table_cols[icol][irow])
        except ValueError:
            table_row += ' & ' + '{}'.format(table_cols[icol][irow])
    table_row += r' \\'
    table_rows.append(table_row)
    

[print(r) for r in table_rows]