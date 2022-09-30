# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:39:14 2021

@author: R.Sechi
This notebook returns a set of artificial data as in 
doi: 10.1177/0003702816662600


"""
import numpy as np
from gendata import gendata



def artdata(sep, noise_level):
    """# % Generate 5 species example data
    'sep' value in [0,1], 0 means no separation at all, 1 reasonable
        separation
     'noise_level'  relative noise level added to data
     
     For Python version:
         Input: sep = float, [0,1]
                noise_level = float, [0,1]
        Output:
            Mtilde= matrix with time-resolved Raman specturm
            W = 2d-array, spectral amplitudes
            H = 2d-array, curves with the concentration proportions
            N = 2d-array, matrix with relative noise level, added to M
            K = 2d-array, transition probability matrix obtained from H
            f = 1d-array, frequencies
            tspan = 2d-array, the time points to each spectrum
            
         """
    H = np.loadtxt('PythonData/H.txt')

    W, H, K, f, tspan = gendata(sep, H)
    M = W.dot(H) 
    print (M.shape)
    N = np.random.random(M.shape)
    N = (noise_level / np.linalg.norm(N)) * N * np.linalg.norm(M)
    N = np.abs(N)
    Mtilde = M + N

    return(Mtilde, W, H, N, K, f, tspan)

