# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 12:37:05 2017

@author: Brittany
"""

import numpy as np
import scipy as sp
from scipy import linalg, special
import math
import matplotlib.pyplot as plt
import warnings
import produceW

class MyException(Exception):
    pass

def create_P(N, L_left, L_right, p_AVG):
    print("creating P")
    
    threshold = 1e100
    
    p = np.zeros(N)
    # p= first row of the matrix P, when i=1
    # P(i, j)= probability of connection from j to i.  
    #  This depends on the value of i-j mod N, 
    #  since we are using a ring model to define the p's
    # Then P is circulant, so we can construct P from its first row.
    # We'll set the max possible value of any P(i, j) to be p_max (see below)
    
    for i in range(1,N):
        if L_right == 0:
            v = (N-i)/L_left
        elif L_left == 0:
            v = i/L_right
        else:
            v = min((N-i)/L_left, i/L_right)
            
        if v <= threshold:
            p[i] = math.exp(-v)
    
    p_avg_current = (1/(N-1))*np.sum(p) #current average value for p
    
    #scale values to make average p value be p_AVG
    for i in range(N):
        p[i] = (p_AVG*p[i])/p_avg_current
    
    p_avg_check = (1/(N-1))*np.sum(p) #check that new average p value is about p_AVG
    
    if p_avg_check>1:
        raise ValueError("Average p value is greater than 1")
    
    p_max = np.ndarray.max(p)
    p_min = np.ndarray.min(p)
    print("p_max = {0}".format(p_max))
    
    p_valid = True #make sure our entries in P are between 0 and 1
    if (p_max>1) or (p_min<0):
        valid_p = False
        raise ValueError("Some probability ouside of [0,1]")
    
    # Make the P matrix from the p vector
    p_short = p[1::]
    p_new = np.concatenate(([p[0]], p_short[::-1]), 0)
    P = linalg.toeplitz(p_new, p)
    
    
    ## Generate the sonet matrix with spatial structure from P
    if p_valid:
        return(P)
    else:
        raise MyException("P was not created.")