# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:32:09 2017

@author: Brittany
"""
import numpy as np

def addMat(double[:,:] X, double[:,:] Y):
    cdef int n = X.shape[0]
    cdef int m = X.shape[1]
    cdef double[:,:] result = np.zeros((n,m))
    cdef int i, j
    
    for i in range(n):
        for j in range(m):
            result[i,j] = X[i,j] + X[i,j]
    return result