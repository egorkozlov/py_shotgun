#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:06:41 2020

@author: egorkozlov
"""




from fit_plot import FitPlots
fp = FitPlots(targ_mode='low education',
               compare='high school baseline.pkl',
               base='high school 2O percent child support.pkl',
               base_name='20 percent child support',
               compare_name='Baseline',
               #graphs_title_add="Experiment: Child Support for Divorced",
               moments_aux=None)#,moments_aux=moments_aux)