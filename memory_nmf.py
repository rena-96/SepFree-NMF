#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 15:44:21 2021

@author: R.Sechi
The difference between this function and memory_nmf_new is that the final plot
 in mme here does not contain the sum of the values of H for each timestep. 
 mme= Minimal Memory Effect
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import string
import seaborn as sns
from scipy.linalg import pinv
from method_nmf import nmf




def memory_nmf(H_r):#, pi="uniform"):
    """Compute the minimal memory effect as determinant of the 
    overlap matrix S. S is obtained from H. Uniform distribution is used.
    Input: 
        H = 2d-array
    Output:
        S= overlap matrix ,2d-array
        det(S) = minimal memory effect, scalar
   """ 
    rank = H_r.shape[1] 
    dim = H_r.shape[0]
    pi = np.full(dim, 1./dim)   
#  for options, you can also compute the transition matrix K and then the 
# stationary distribution
    num = H_r.T.dot(np.diag(pi).dot(H_r))
    den = H_r.T.dot(np.diag(pi).dot(np.ones((dim,1))))
    den = den*np.eye(rank)  
    S_c = pinv(den).dot(num)
    return(S_c, np.linalg.det(S_c))

def mme(data, timer, wn, params=[1.,1.,1.,1.,1.], clus_list=[3]):
    """Compute the minimal memory effect and plot the result of the NMF
    without separability assumption as 2 plots, one for H, and one for W.

    One can try different number of clusters for the rank of the nmf.
    Input:
        data= 2d-array, dataset to decompose
        timer= 1d-array, time values for the data
        wn = 1d-array wavenumbers
        params = list, parameters for the objective function
        clus_list = list, list fof ranks for the NMF
    Output:
        plots= figures
        list_memory= 2d-array, number of clusters and corresponding
        minimal memory effect , det(S)."""
    list_memory = []
    maxclus = max(clus_list)
    for nclus in clus_list:
        
        
        M_a, W_a, H_a, P_a, _,Chi_a, _ = nmf(data, r=nclus, params=params )
         #memory 
        Sa, detSa = memory_nmf(H_a.T)
        list_memory.append((nclus, detSa))
        
        #%% plot
        plt.figure(figsize=(15,4))
        plt.suptitle('SepFree NMF', fontsize=16)
        labels= list(string.ascii_uppercase)

        cmap = ListedColormap(sns.color_palette("hls", int(maxclus+1)).as_hex())#plt.get_cmap("tab20")
        plt.subplot(1, 2, 1)
        plt.title("Plot $H$")
        for i in range(len(H_a)):
            plt.plot(timer,H_a[i], '--',linewidth=2, color=cmap((i)),label=labels[i])
            #plt.xticks(np.arange(len(wl), step=120))#,labels=(np.round(wl[/1000)))
        plt.grid()
        plt.legend()

        plt.xlabel("t[s]")#"$\lambda$/nm")
        plt.ylabel("concentration proportion")
        plt.subplot(1, 2, 2)
        plt.title("Plot $W$")
        for i in range(nclus):

            plt.plot(wn,W_a.T[i], "-",linewidth=2, color=cmap(i), label=labels[i])
    #    plt.xticks(np.arange(len(wn), step=300),labels=(np.round(wn[::300])))

        plt.legend()
        plt.grid()

        plt.xlabel(r"$\nu$[cm-1]")
        plt.ylabel("intensity [a.u.]")
       # plt.gca().invert_xaxis() 
      
     #   plt.savefig("figure.pdf")
        # if savefigs==True:
        current_values = plt.gca().get_yticks()
        plt.gca().set_yticklabels(['{:,.0f}'.format(x) for x in current_values])  
        plt.savefig("figure.pdf")
        plt.show()
        #%% sentence with the minimal memory effect
        print("The minimal memory effect for %d components is %f."% (nclus, detSa))
    list_memory = np.array(list_memory, dtype=[("compounds", int), ("detS", float)])
    if len(clus_list)>1:
        print("The best result is ", np.sort(list_memory, order='detS')[:-1]  )
    return(list_memory)