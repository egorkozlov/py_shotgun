#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:26:49 2020

@author: egorkozlov
"""

# this produces table about the child support
from tiktak import filer


names = [('No','no child support'),
         ('BL','baseline'),
         ('Full','full child support'),
         ]





entries = [(r'\% divorced in 10 years if KF','divorced by years after marriage if kids first, 10'),
           (r'\% divorced in 10 years if MF','divorced by years after marriage if marriage first, 10'),
           (r'\% KF women in population at 30','k then m in population at 30'),
           (r'\% MF women in population at 30','m then k in population at 30'),
           (r'\% single mothers among mothers at 35','single mothers among mothers at 35'),
           (r'\% unplanned pregnancies aborted','unplanned pregnancies aborted'),
           #(r'\% with kids 5 years after marriage','ever kids by years after marriage, 5')
           ]

table_cols = list()
for educ in ['col','hs']:
    
    def file_name(nm): return '{} {}.pkl'.format(educ,nm)

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
    for icol in range(2*len(names)):
        try:
            table_row += ' & ' + '{:02.1f}'.format(100*table_cols[icol][irow])
        except ValueError:
            table_row += ' & ' + '{}'.format(table_cols[icol][irow])
    table_row += r' \\'
    table_rows.append(table_row)
    

[print(r) for r in table_rows]