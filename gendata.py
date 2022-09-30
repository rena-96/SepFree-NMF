# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 20:49:47 2021

@author: R.Sechi
"""
import numpy as np
import matplotlib.pyplot as plt
from densify import densify
from lorentz_parafun import lorentz_parfun
#%%

def gendata(sep, H):
    """ Generates synthetic dataset as in doi: 10.1177/0003702816662600, 
    without adding relative noise.
    Input: 
        sep= float [1,0], separation grade of the Raman peaks, 
        1-> well-separated peaks
        0 -> high interference
        H= 2d-array, matrix with the kinetic curves
    Output:  
        W= 2d-array, matrix with the peaks amplitudes
        H= 2d-array, same as input
        K= 2d-array, Koopman transition probability matrix for timestep 
        between the measurement of the spectra
        f= 1d-array, frequencies
        tspan= 1d-array, time values of the spectra"""
    f = np.linspace(300, 3210, 2478)
    
#function [W, H, K, f, tspan] = gendata(sep)
#f = linspace(300, 3210, 2478);

    w1_params = 1e4 *np.array( [    0.307422784540460  , 0.000392574246463,   1.455928866193969,
   0.156096159915473  , 0.000303129027706,   0.746998370806083,
   0.127220337949040  , 0.000286920793291,   1.135618326874843,
   0.084557996408843  , 0.000401671696353,   0.301080121906525,
   0.315090291702677  , 0.000650154475367,   0.185087973214349,
   0.169869658279571  , 0.000813013193150,   0.103636040755365,
   0.076091619902954  , 0.000974211385697,   0.023916315960571 ])
    w1_params = np.reshape(w1_params, (7, 3))
    w1_params [:, 2] = w1_params[:,2]*0.5
    w2_params = 1e6 *np.array( [
   0.001001997197771 ,  0.000001583113834  , 0.032730098776435,
   0.003067719328406 ,  0.000003512331970  , 0.007758016743085,
   0.002930879425168 ,  0.000002513680550  , 0.012180441795462,
   0.000783939422285 ,  0.000002116995928  , 0.010083679048010,
   0.001030186834284 ,  0.000002054031885  , 0.009137234154511,
   0.001210925390462 ,  0.000002715914269  , 0.008070027466481,
   0.003037634593434 ,  0.000003769524125 ,  0.004005125640032,
   0.000518469745310 ,  0.000002344625153 ,  0.003030720838585,
   0.001607087033096 ,  0.000002281414662 ,  0.002282018828316,
   0.003056806406564 ,  0.000002694841690 ,  0.018484666236737,
   0.002983268076895 ,  0.000002817633371 , 0.002897358120990,
   0.001003844748654 ,  0.000000067631098 ,  0.057694208562437,
   0.001380834573228 ,  0.000005338844681 ,  0.001801580464834,
   0.003004121661074 ,  0.000005898130686 ,  0.001703208251031,
   0.002737568107921 ,  0.000005778499949 ,  0.001057380893449,
   0.001005928860975 ,  0.000000807063190 ,  0.007113431873642,
   0.000786000499083 ,  0.000001309151247 ,  0.006769313017552,
   0.003170980382966 ,  0.000006323096927 ,  0.000785221929041,
   0.001164191963878 ,  0.000026408785803 ,  0.000477796702030,
   0.002949538851919 ,  0.000003663267553 ,  0.001056383849359,
   0.000620012356222 ,  0.000004060900255 ,  0.000694034684541,
   0.000990670161182 ,  0.000001697274533  , 0.000630434517332,
   0.003208397199561 ,  0.000001000014058  , 0.000522190463200])
    w2_params = np.reshape(w2_params, (23,3))

#% This scaling makes the first species "small", and hence presumably more
#% difficult to detect.
    #w1_params(:,3) = w1_params(:,3) * 0.5;

    w3_params = 1e4 * np.array([
   0.293945369626981 ,  0.000531123203595 ,  0.521364334996143,
   0.287198137622590 ,  0.001033758426059 ,  1.790983102603328,
   0.298034442860098 ,  0.001224548674735 ,  0.470767015452515,
   0.280767999478192 ,  0.000893403586983 ,  0.267857293104895,
   0.290209937186533 ,  0.001080953013198 ,  0.303430225300156,
   0.043669478664494 ,  0.000453747630641 ,  0.120499210353156,
   0.285635317497733 ,  0.000807314848781 ,  0.356568115385855,
   0.295409971766293 ,  0.000675458576517 ,  0.209740173500541,
   0.274362001446827 ,  0.005148656399637 ,  0.071922657361501,
   0.145753835330227 ,  0.001746520520333 ,  0.082008369871866,
   0.084025878472849 ,  0.000937600786885 ,  0.070440008434370,
   0.293180001037005 ,  0.000845126156928 ,  1.309204886313470  ])
    w3_params = np.reshape(w3_params, (12,3) )


    w4_params = 1e4 * np.array([
   0.308164955649732,   0.000554930617972 ,  1.116752983140727,
   0.158739318597909,   0.000387956740385  , 0.586217679753232,
   0.070769588971476,   0.000427541652898  , 0.494252445204085,
   0.118275013268578,   0.000541676202714  , 0.313050740196735,
   0.316831647441829,   0.000806453748109 ,  0.230150657611183,
   0.158981550906583,   0.000170494408898 ,  0.866012663163046,
   0.040285979894797,   0.000790086394983 ,  0.112112578929518,
   0.169864773760108,   0.000872055378499 ,  0.104588075197453,
   0.056400352796060,   0.000897629570574 ,  0.065709315810195,
   0.071171110476781,   0.000271058181285 ,  0.340032030763652,
   0.159188586492578,   0.000115734767329 ,  0.529726929535237])
    w4_params = np.reshape(w4_params, (11,3) )

    w5_params = 1e4 * np.array( [
   0.091888102684603 ,  0.000314344091520 ,  0.301945700885340,
   0.137695726793342 ,  0.000866225203442 ,  0.084252090010115,
   0.225563055838154 ,  0.000388589961307 ,  2.560132887473927,
   0.229587731006686 ,  0.000577642944510 ,  0.141080809848103,
   0.273521443100234 ,  0.000459848668289 ,  0.059880517019103,
   0.294330053137586 ,  0.000359290008367 ,  1.179144636951206,
   0.300561698572789 ,  0.000445580666864 ,  0.172202166910520])
    w5_params = np.reshape(w5_params, (7,3) )


    W = np.zeros((len(f), 5))
 #   print(densify(w1_params, sep))
    W[:,0] = lorentz_parfun(densify(w1_params, sep), f)
    W[:,1] = lorentz_parfun(densify(w2_params, sep), f)
    W[:,2] = lorentz_parfun(densify(w3_params, sep), f)
    W[:,3] = lorentz_parfun(densify(w4_params, sep), f)
    W[:,4] = lorentz_parfun(densify(w5_params, sep), f)

    K = np.array([[-0.53,  0.02,  0.0,   0.,  0.0], 
     [ 0.53, -0.66,  0.25,  0,  0.0],
     [ 0.0,   0.43, -0.36,  0,  0.1],
     [ 0.0,   0.21,  0.0,   0,  0.0],
     [ 0.0,   0.0,   0.11,  0, -0.1]])
    tspan = np.linspace(0,50,151)
    
   # H = compute_exact_kinetic(tspan, K);
   
    # tmp = W
    # W = H.T
    # H = tmp.T
    return( W, H,  K, f,  tspan)
