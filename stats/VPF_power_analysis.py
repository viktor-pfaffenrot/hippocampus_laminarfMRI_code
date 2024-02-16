#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 12:15:34 2024

@author: pfaffenrot
"""

import scipy.io
import scipy.stats
import statsmodels.stats.power as smp
import matplotlib.pyplot as plt
import numpy as np
import math

def pooled_std(s1,s2):
    return np.sqrt((s1**2 + s2**2)/2)

mat = scipy.io.loadmat(
    "/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/avg_4_power_analysis_pre_vs_post.mat"
)
Xm, Xs = (mat["Xm2"], mat["Xs2"])

subfields = ['subiculum','CA1','CA2','CA3','CA4/DG']

power_analysis = smp.TTestPower()
power = []
sample_size = []
for ii in range(0,5):
    effect_size =  (Xm[ii])/(pooled_std(Xs[ii],0))

    
    power.append(power_analysis.power(effect_size=effect_size, nobs=9, alpha=0.05))
    sample_size.append(math.ceil(power_analysis.solve_power(
    effect_size=effect_size, power=0.8, alpha=0.05
    )))
