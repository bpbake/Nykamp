# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:32:09 2017

@author: Brittany
"""
import numpy as np
import sys

def produceW(double[:,:] A, double[:,:] B, double[:,:] M_tilde, double[:,:] M_theta, double c, double d, double e, int N):
    cdef int i, j, n
    #
    cdef double[:,:] X = np.random.standard_normal((N,N)) # random matrix of iid standard normal entries
    # cdef double[:,:] Z = np.zeros((N,N)) # Z = SX, a random vector with entries iid standard normal
    # Z is Gaussian with covariance matrix approximately Sigma (ones)
    cdef double[:,:] W = np.zeros((N,N)) # dichomatized Z

    cdef double[:,:] MX = np.multiply(M_tilde, X) # entry ij is sigma_tilde_ij*X_ij
    # cdef double[:,:] cM = np.multiply(c,M_tilde)
    # cdef double[:,:] dM = np.multiply(d,M_tilde)
    # cdef double[:,:] eM = np.multiply(c,M_tilde)
    cdef double[:] MXrowsum = np.sum(MX, axis=1) # vector of sums of rows
    cdef double[:] MXcolsum = np.sum(MX, axis=0) # vector of sums of columns
    # cdef double[:,:] sumCDE = np.zeros((N,N))

    cdef double cMij, dMij, eMij
    cdef double sumCDEij
    cdef double result
    
    for i in range(N):
        for j in range(N):
            cMij = c*M_tilde[i,j]
            dMij = d*M_tilde[i,j]
            eMij = e*M_tilde[i,j]

            sumCDEij = cMij*MXrowsum[i] + dMij*MXcolsum[j] + eMij*(MXrowsum[j] + MXcolsum[i])

            result = A[i,j]*X[i,j] + B[i,j]*X[j,i] - (cMij*MX[i,j] + dMij*MX[i,j] + 2*eMij*MX[j,i]) + sumCDEij

            # for n in range(N):
            #     if (n != i) and (n != j):
            #         result += cMij*MX[i,n] + dMij*MX[n,j] + eMij*(MX[j,n] + MX[n,i])

            if result >= M_theta[i,j]:
                W[i,j]=1
            else:
                W[i,j]=0
                
        if (i%100 == 0):
            print("row={0} out of N={1}".format(i,N))
            sys.stdout.flush()
    print("\n")
    return W