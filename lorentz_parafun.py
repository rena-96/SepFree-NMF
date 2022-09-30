# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 18:38:54 2021

@author: R.Sechi
"""
import numpy as np
from lorentz import lorentz 

def lorentz_parfun(params, x):
    """% y = lorentz_parfun(params, x)
%
% Evaluate a sum of Lorentzians defined through the params array at x.
%
% Each row of 'params' has length 3 with the following meaning:
%
%   params(:,1)  --  x0, the base point
%   params(:,2)  --  gamma, the half-width at half height
%   params(:,3)  --  I, the intensity
%
% See also: lorentz"""
    num_1 , n = np.shape(params)
    assert(n==3)
    y = np.zeros(np.shape(x))
    for k in range(num_1):

        x0 = params[k, 0]
        gamma = params[k, 1]
        I = params[k, 2]
        f = lorentz(x, x0, gamma,I)
        y = y + f
    return y
