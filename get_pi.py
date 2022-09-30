#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:16:57 2021

@author: bzfsechi
"""

import numpy as np
from scipy.linalg import eig, pinv
from cmdtools.analysis import pcca
#%%
def get_pi(M, pi="uniform"):
    """Input:
        M: 2D array, matrix
        pi: string, uniform (uniform distribution)
            or statdistr (stationary distribution), 
            eigenvector to eigenvalue 1"""
    if pi == "uniform":
        dim = np.size(M, 1)
        pi = np.full(dim, 1./dim)
    elif pi == "statdistr":
        ews, levs = eig(M, left=True, right=False)
        inds = np.argsort(abs(ews.real))
        levs = levs[:,inds]
        pi = levs[:,-1]
 
    return pi/pi.sum()

# def rebinding(M, nclus, pi="uniform"):
#     """, paper marcus and max 12-13"""
#     dim = np.shape(M)[0]
#     pi = get_pi(M, pi=pi)
#     #print(pi)
#     chi = pcca.pcca(M, nclus)
#     num = chi.T.dot(np.diag(pi).dot(chi))
#     den = chi.T.dot(np.diag(pi).dot(np.ones((dim,1))))
#     den = den*np.eye(nclus)

#     S_c = pinv(den).dot(num)
#     return(S_c, np.linalg.det(S_c))
    
# def rebinding_nmf(H_r):#, pi="uniform"):
#     """h_r=chi, so the stiffness matrix is given by 
#     h_r.T.pi.h_r/(h_r.T.pi.e),
#     TODO: stiffness with initial distribution"""    
#     rank = H_r.shape[1] 
#     dim = H_r.shape[0]
#     pi = np.full(dim, 1./dim)   
#     print(rank, dim)
#     num = H_r.T.dot(np.diag(pi).dot(H_r))
#     den = H_r.T.dot(np.diag(pi).dot(np.ones((dim,1))))
#     den = den*np.eye(rank)
#     print(num.shape, den.shape)
#     S_c = pinv(den).dot(num)
#     return(S_c, np.linalg.det(S_c), den, num)
        
    