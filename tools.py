#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:32:33 2020

@author: bzfsechi
"""

import numpy as np


def pi_pcca(lambdas):
    """Compute the weights of the time-resolved spectrum for a specific 
    value of lambda to consider in the PCCA+ algorithm. The smallest 
    wavelength difference is taken as unity and the differences between 
    walengths are scaled to that smallest value. Return the density pi. """
    diff = abs(lambdas[:-1]-lambdas[1:])
    min_lambda = np.amin(diff)
    pi = np.zeros((len(diff)+1))
    pi[0] = diff[0]/min_lambda
    pi[-1] = diff[-1]/min_lambda
    for i in range(1, len(pi)-1):
        pi[i] = (diff[i-1]+diff[i])/(2*min_lambda)
    return pi/np.sum(pi)
