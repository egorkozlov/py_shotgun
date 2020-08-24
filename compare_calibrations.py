#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:06:41 2020

@author: egorkozlov
"""

from estimates import get_point

high_e = True
x, targ_mode = get_point(high_e)


from fit_plot import FitPlots


what = 'no child support'

fp = FitPlots(targ_mode=targ_mode,
               base='col {}.pkl'.format(what),
               compare='col baseline.pkl',
               base_name=what,
               compare_name='baseline',
               #graphs_title_add="Experiment: Removing Subsistence Constraint",
               moments_aux=None)#,moments_aux=moments_aux)