# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 15:16:00 2021

@author: R.Sechi
This notebook generates the syntethic dataset with no noise nor interference.
Then, the minimal memory effect is computed for the synthetic dataset
for a NMF 5 compounds.
The generated synthetic dataset represents a time-resolved Raman spectrum.
"""
import numpy as np
import matplotlib.pyplot as plt
from artdata import artdata
from method_nmf import nmf
#from memory_nmf_new import mme_new
from memory_nmf import mme 
from matplotlib.colors import ListedColormap
import string
import seaborn as sns

#%%1 no noise, no interference
#generates dataset
M, W, H, N, K, f, tspan = artdata(1., 0.)

#plot of W and H of the synthetic dataset
plt.figure(figsize=(15,4))
plt.suptitle('Synthetic dataset', fontsize=15)
labels= list(string.ascii_uppercase)

cmap = ListedColormap(sns.color_palette("hls", int(5+1)).as_hex())#plt.get_cmap("tab20")
plt.subplot(1, 2, 1)
plt.title("Plot $H$")
for i in range(len(H)):
    plt.plot(tspan,H[i], '--',linewidth=2, color=cmap((i)),label=labels[i])
    #plt.xticks(np.arange(len(wl), step=120))#,labels=(np.round(wl[/1000)))
plt.grid()
plt.legend()

plt.xlabel("t[s]")#"$\lambda$/nm")
plt.ylabel("concentration proportion")
plt.subplot(1, 2, 2)
plt.title("Plot $W$")
for i in range(5):

    plt.plot(f,W.T[i], "-",linewidth=2, color=cmap(i), label=labels[i])
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
plt.savefig("trueHW.pdf")
plt.show()






#uncomment the line below if you want to obtain all the components of the 
#decomposition
#M_r, W_r, H_r, K_r, A_optimized, chi, Uitgm = nmf(M,  r = 5, params = [-0.0001,-1.,1.,0.,0.], weight=False)
#%%
#compute minimal memory effect
# save figure under the name "figure.pdf"
detS = mme(M, tspan, f, params= [-0.001,-10.,1.,0.,0.], clus_list=[5])
#%%
#plot of the generated dataset as bidimensional heatmap
plt.figure(figsize=(10,10))
plt.imshow(M.T, aspect="auto")
plt.yticks(np.arange(0,tspan.shape[0], step=20), np.around(tspan[::20],1))
plt.xticks(np.arange(0,f.shape[0], step=300), labels= np.around(f[::300]))
plt.colorbar()
plt.xlabel(r"$\nu$ [cm-1]")
plt.ylabel("time delay [s]")
plt.title("Synthetic Raman spectrum")
#plt.savefig("synthetic_data.pdf")
plt.show()
#%%

#uncommenting the following lines, you will see the comparison between the D 
#species of the artificial generated dataset and D and E obtained
# with the algorithm.
# plt.plot(f,W[:,3],label = "D artificial \n data" )
# plt.plot(f, (detS[1])[:,3]+1000,label = "D reconstructed" )
# plt.plot(f, ( detS[1])[:,4]+2000, label = "E reconstructed ")
# plt.xlabel(r"$\nu$ [cm-1]")
# plt.ylabel("Raman intensity")
# plt.legend()
# #%%
# plt.plot(tspan, H[3,:],label = "D artificial \n data" )
# plt.plot(tspan, (detS[2])[3,:],label = "D reconstructed" )
# plt.plot(tspan, ( detS[2])[4,:], label = "E reconstructed ")
# plt.xlabel(r"time delay ")
# plt.ylabel("state proportion ")
# plt.legend()