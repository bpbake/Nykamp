# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:32:09 2017

@author: Brittany
"""
import numpy as np

def produceW(double[:,:] A, double[:,:] B, double[:,:] M_tilde, double[:,:] M_theta, double c, double d, double e, int N):
    cdef int i, j, n
    #
    cdef double[:,:] X = np.random.standard_normal((N,N)) # random matrix of iid standard normal entries
    # cdef double[:,:] Z = np.zeros((N,N)) # Z = SX, a random vector with entries iid standard normal
    # Z is Gaussian with covariance matrix approximately Sigma (ones)
    cdef double[:,:] W = np.zeros((N,N)) # dichomatized Z

    cdef double[:,:] MX = np.multiply(M_tilde, X)
    cdef double cMij, dMij, eMij
    cdef double result
    
    for i in range(N):
        for j in range(N):
            cMij = c*M_tilde[i,j]
            dMij = d*M_tilde[i,j]
            eMij = e*M_tilde[i,j]

            result = A[i,j]*X[i,j] + B[i,j]*X[j,i]

            for n in range(N):
                if (n != i) and (n != j):
                    result += cMij*MX[i,n] + dMij*MX[n,j] + eMij*(MX[j,n] + MX[n,i])
            if result >= M_theta[i,j]:
                W[i,j]=1
            else:
                W[i,j]=0
        if (i%100 == 0):
            print("row={0} out of N={1}".format(i,N))
    print("\n")
    return W