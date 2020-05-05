#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 12:00:54 2020

@author: egorkozlov
"""

import numpy as np

def get_age_distribution(flist=None):
    if flist is None:
        flist = [847361,3631473,5347413,5936839,6803464,6742099,6927052,
                  6926038,6845848,7243977,6927451,6881834,6715119,6617213,
                  6771984,6490318,6414025,6401949,6327301,6693319]
    age_distribution = np.array(flist)
    
    age_dist_normailzed = age_distribution/age_distribution.sum()
    return age_dist_normailzed