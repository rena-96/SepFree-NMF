# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 15:16:00 2021

@author: R.Sechi
This notebook generates the syntethic dataset with added normal-distributed noise
 and whose peaks are not well-separated.
Then, the minimal memory effect is computed for the synthetic dataset
for a NMF into 5 compounds.
The generated synthetic dataset represents a time-resolved Raman spectrum.
"""
import numpy as np
import matplotlib.pyplot as plt
from artdata import artdata
from method_nmf import nmf
from memory_nmf import mme

#%% noise, interference
M, W, H, N, K, f, tspan = artdata(.7, 0.4)
#%%
#uncomment the line below if you want to obtain all the components of the 
#decomposition
#M_r, W_r, H_r, K_r, A_optimized, chi, Uitgm = nmf(M,  r = 5, params = [-0.0001,-1.,1.,0.,0.], weight=False)
#%%
#compute minimal memory effect
detS = mme(M, tspan, f, params= [-0.001,-10.,1.,0.,0.], clus_list=[5])