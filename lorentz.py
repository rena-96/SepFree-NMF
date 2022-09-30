# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 14:57:26 2021

@author: R.Sechi
"""

import numpy as np
def lorentz(x, x0, gamma, I):
    """% Return a function handle to a Lorentzian function with given parameters
    %Input:
      x0    The center of the peak
      gamma HWHM, the half width at half maximum
      I     The intensity"""
    
    return I*gamma**2/( (x-x0)**2 + gamma**2 )
    
