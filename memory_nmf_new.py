#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 15:44:21 2021

@author: R. Sechi
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

def mme_new(data, timer, wn, params=[1.,1.,1.,1.,1.], clus_list=[3]):
    """Compute the minimal memory effect and plot the result of the NMF
    without separability assumption as 3 plots, one for H, one for the
    sum of the H values for each timestep, and one for W.
    We want the sum of H to be as close to one as possible.
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
        
        
        M_a, W_a, H_a, K_a, _,Chi_a, _ = nmf(data, r=nclus, params=params )
        
         #memory 
        Sa, detSa = memory_nmf(H_a.T)
        list_memory.append((nclus, detSa))
        
        #%% plot
        fig = plt.figure(figsize=(15,4))
        plt.suptitle('SepFree NMF', fontsize=16)
        labels= list(string.ascii_uppercase)
        cmap = ListedColormap(sns.color_palette("hls", int(maxclus+1)).as_hex())#plt.get_cmap("tab20")
        
        ax1 = fig.add_subplot(122)
        ax1.set_title("Plot $W$")
        for i in range(nclus):
            ax1.plot(wn,W_a.T[i], "-", linewidth=2,color=cmap(i), label=labels[i])
    #    plt.xticks(np.arange(len(wn), step=300),labels=(np.round(wn[::300])))

        plt.legend()
        plt.grid()

        ax1.set_xlabel(r"$\nu$[cm-1]")
        ax1.set_ylabel("intensity [a.u.]")
        # other column
        gs = fig.add_gridspec(3, 2, hspace=0)
        ax2 = fig.add_subplot(gs[:2, 0])
        ax2.set_title("Plot $H$")
        for i in range(len(H_a)):
           # ax2.plot(timer,H_a[i]/abs(H_a).sum(axis=0), '--', color=cmap((i)),label=labels[i])
            ax2.plot(timer,H_a[i], '--', linewidth=2,color=cmap((i)),label=labels[i])
        
            #plt.xticks(np.arange(len(wl), step=120))#,labels=(np.round(wl[/1000)))
        plt.grid()
        plt.legend(loc=1)
       
        ax2.set_ylabel("concentration proportion")
        ax3 = fig.add_subplot(gs[2, 0])
        ax3.plot(timer,H_a.sum(axis=0),linewidth=2, label= "$\Sigma$ $H_t$")
        ax3.set_ylim(0.8, 1.2)
        plt.grid()
        plt.legend()
        ax3.set_xlabel("t[s]")#"$\lambda$/nm")
        ax3.label_outer()
        ax2.label_outer()
        
       # plt.gca().invert_xaxis() 
      
     #   plt.savefig("methanol_3_sec_%i.pdf"%nclus)
        # if savefigs==True:
            
        plt.show()
        #%% sentence with the minimal memory effect
        print("The minimal memory effect for %d components is %f."% (nclus, detSa))
    list_memory = np.array(list_memory, dtype=[("compounds", int), ("detS", float)])
    if len(clus_list)>1:
        print("The best result is ", np.sort(list_memory, order='detS')[-1]  )
    return(list_memory, W_a, H_a )