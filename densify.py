# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 20:33:03 2021

@author: R.Sechi
This notebook generates interference between the Lorentzians curves by 
shifting them together.
"""
import numpy as np
#
def densify(w, sep):
    """  Input: w = 1d-array, center of the Lorenztians curves
            sep = float, [0,1], if 1: well-separated, if 0: high interference
       Output: w = 1d-array, modified position of the centers of the peaks
    This function shifts the position of the peaks w. 
    """
    f = np.linspace(300, 3210, 2478)
    f_breaks = np.linspace(np.amin(f), np.amax(f), 7)
    f_mid = f_breaks[1::2]
    f_beg = f_breaks[::2]
    num_peak = w.shape[0]
    for k in range(num_peak):
   
        which = np.nonzero(f_beg > w[k, 0])
       # print(which[0][0])
        w[k, 0] = w[k, 0] + (1. - sep) * (f_mid[int(which[0][0])-1] - w[k, 0])
    return (w)    
    


