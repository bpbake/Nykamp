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

class MyException(Exception):
    pass

def createW(N, L, p_AVG, alpha_recip, alpha_conv, alpha_div, alpha_chain):
    threshold = 1000
    
    p = np.zeros(N)
    # p= first row of the matrix P, when i=1
    # P(i, j)= probability of connection from j to i.  
    #  This depends on the value of i-j mod N, 
    #  since we are using a ring model to define the p's
    # Then P is circulant, so we can construct P from its first row.
    # We'll set the max possible value of any P(i, j) to be p_max (see below)
    
    for i in range(N):
        M = (-i)%N
        if 1 <= M and M < (N/2): # going left/counter-clockwise (i.e., when i = 1, j=501-100)
            if abs(M/L) > threshold:
                p[i] = 0
            else:
                p[i] = math.exp(-M/L)
        elif N/2 <= M and M < N: # going right/clockwise (i.e., when i = 1, j=2-500)
            if abs((M-N)/L) > threshold:
                p[i] = 0
            else:
                p[i] = math.exp(-(N-M)/L)
    
    p_avg_current = (1/(N-1))*np.sum(p) #current average value for p
    
    #scale values to make average p value be p_AVG
    for i in range(N):
        p[i] = (p_AVG*p[i])/p_avg_current
    
    p_avg_check = (1/(N-1))*np.sum(p) #check that new average p value is about p_AVG
    
    if p_avg_check>1:
        raise ValueError("Average p value is greater than 1")
    
    p_max = np.ndarray.max(p)
    p_min = np.ndarray.min(p)
    
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
        c1=0.2658  #coefficients for best fit equation of rho as a function of p1, p2, alpha
        c2=0.0470  #rho(p1,p2,alpha)~c3*(p1+c1)*(p2+c1)*(alpha+c2)
        c3=1.9531
        
        M_theta = np.zeros((N,N)) # M_theta(i,j) = theta_ij 
        # matrix of thresholds.  If Z_ij > theta_ij then W_ij = 1
        for i in range(N):
            for j in range(N):
                M_theta[i,j] = math.sqrt(2)*special.erfinv(1 - (2*P[i,j]))

        Msigma = np.ones((N,N)) - np.identity(N)# Matrix of sigmas.

        M_tilde = np.multiply(Msigma, np.add(P,c1))# M_tilde(i,j) = sigma_ij^tilde
        # we multiply sigma_ij by (p_ij + c1) to get sigma_ij^tilde

        m = np.mean(np.sum(np.square(M_tilde), axis=0)) # average row/column sum of M_tilde (the sigma tilde squared)
        # Each row/column sum should be the same, but we're averaging just because
        
        # now, we will solve for c, d, and e using (AA)*(AA)^T = E
        E = np.array([[(c3*(alpha_conv +c2)), (c3*(alpha_chain +c2))], [(c3*(alpha_chain +c2)), (c3*(alpha_div +c2))]])
    
        D, V = np.linalg.eig(E)
        # D is a vector of the eigenvalues, each repeated according to its multiplicity. The eigenvalues are not necessarily ordered. 
        #   The resulting array will be of complex type, unless the imaginary part is zero in which case it will be cast to a real type. 
        #   When a is real the resulting eigenvalues will be real (0 imaginary part) or occur in conjugate pairs
        # V is an array of the corresponding normalized (unit “length”) eigenvectors, 
        #   such that the column V[:,i] is the eigenvector corresponding to the eigenvalue D[i]
        
        BB = np.zeros((2, 2))
        
        for i in range(2):
            if (D[i] < 0) and (D[i] > -1e-12):
                D[i] = 0 # If eigen value is really small and negative, we can assume there's a roundoff error and it should be zero
        
        if (D[0] >= 0) and (D[1] >= 0):
            BB[0,0] = math.sqrt(D[0])
            BB[1,1] = math.sqrt(D[1])
        else:
            raise MyException("Some eigenvalues in D are negative.") #If eigen value is smaller (more negative) than -1e-12, then we claim to have a problem
        
        AA = np.matmul(np.matmul(V,BB), np.linalg.inv(V)) 
        c = AA[0,0]/np.sqrt(m)
        d = AA[1,1]/np.sqrt(m)
        e = AA[0,1]/np.sqrt(m)
        
        # Now we'll define F, G, and H
        # Lets's make things easier by defining these:
        M_square = np.square(M_tilde) # M_square[i,j] = sigma_ij^tilde^2
        m_sq = np.subtract(m,M_square) # m_sq[i,j] = m - sigma_ij^tilde^2
        m_sqT = np.subtract(m,np.transpose(M_square)) # m_sqT[i,j] = m - sigma_ji^tilde^2
        
        F = np.subtract(np.square(Msigma), np.multiply(M_square,np.add(np.multiply(math.pow(c,2), m_sq),np.add(np.multiply(math.pow(d,2), m_sq), np.multiply(2*math.pow(e,2),m_sqT)))))   
        #F[i,j] = A[i,j]^2+B[i,j]^2 = Msigma[i,j]^2 - M[i,j]^2*((c^2 + d^2)*(m-M[i,j]^2) + 2*e^2*(m-M[j,i]^2));
        # note this requires F[i,j] >=0
    
        G = np.transpose(F)#G[i,j] = A[j,i]^2 + B[i,j]^2 = Msigma[j,i]^2 - M[j,i]^2*((c^2 + d^2)*(m-M[j,i]^2) + 2*e^2*(m - M[i,j]^2));
        # note this requires G[i,j]>=0
        
        H = np.multiply(np.multiply(M_square, np.transpose(M_square)),np.subtract(c3*(alpha_recip +c2), np.multiply(e*(c+d),np.add(m_sq, m_sqT))))
        # H[i,j] = 2*A[i,j]*B[i,j] = M[i,j]*M[j,i]*rho_recip - M[i,j]*M[j,i]*(c*e*(m-M[i,j])^2) + c*e*(m-M[j,i]^2) + d*e*(m-M[i,j]^2) + d*e*(m-M[j,i]^2);
        # note sign(H[i,j]) = sign(B[i,j]) since A[i,j]>=0
        
        # note that F[i,j], G[i,j], and H[i,j] do not depend on m because of how c, d, and e depend on m
        # note that F[i,j]=G[j,i] and H[i,j]=H[j,i]
        
        
        A = np.zeros((N, N)) # A[i,j] = a_ij
        B = np.zeros((N, N)) # B[i,j] = b_ij
        
        quit_fg0 = 0 #check: if F[i,j]=0 and G[i,j]=0, then we should have H[i,j]=0
        quit_f0 = 0 #check: if F[i,j]=0 then we should have H[i,j]=0
        quit_g0 = 0 #check: if G[i,j]=0 then we should have H[i,j]=0
        quit_a2f = 0 #check: if A[i,j]=sqrt(F[i,j]) then we should have H[i,j]=0
        quit_neg = 0 #check: discriminant should be non-negative
        # note that above, we already ruled out the cases where F[i,j]<0 or G[i,j]<0
    
        if (F < np.zeros(N)).any():
            raise MyException("F had some negative entries.  A was not created.")
            
        coeff_a = np.square(np.subtract(F,G)) + 4*np.square(H)
        coeff_b = -2*(np.multiply(np.square(H),np.subtract(F,G))) - 2*(np.multiply(F,np.square(np.subtract(F,G)))) - 4*(np.multiply(np.square(H),F))
        coeff_c = np.square(np.square(H)) + 2*(np.multiply(np.multiply(np.square(H),F),np.subtract(F,G))) + np.multiply(np.square(np.subtract(F,G)),np.square(F))
        Disc = np.square(coeff_b) - 4*np.multiply(coeff_a, coeff_c)
    
    
        Asq = np.divide(np.add(-1*coeff_b, np.sqrt(Disc)),2*coeff_a)
        A = np.sqrt(Asq)
        B = np.divide(H, np.add(A,np.transpose(A)))
            
        np.fill_diagonal(A,0)
        np.fill_diagonal(B,0)
            
        if  (Disc < np.zeros(N)).any():
            warnings.warn("discriminant<=0 for some i,j.  alpha_recip not accounted for.")
            for i in range(N):
                for j in range(i):
                    if Disc[i,j] < 0:
                        B[i,j] = 0
                        A[i,j] = math.sqrt(F[i,j])
                    if Disc[j,i] < 0:
                        B[j,i] = 0
                        A[j,i] = math.sqrt(F[j,i])
                        
        checkFmat = np.amax(np.absolute(np.subtract(np.add(np.square(A),np.square(B)),F)))
        print("|A^2 + B^2 - F| leq {0}".format(checkFmat)) 
    
        # Now we get to the fun part where we calculate W
        X = np.random.standard_normal((N,N)) # random matrix of iid standard normal entries
        Z = np.zeros((N,N)) # Z = SX, a random vector with entries iid standard normal
        # Z is Gaussian with covariance matrix approximately Sigma (ones)
        W = np.zeros((N,N)) # dichomatized Z
        
        for i in range(N):
            for j in range(N):
                Z[i,j] = A[i,j]*X[i,j] + B[i,j]*X[j,i]
                for n in range(N):
                    if (n != i) and (n != j):
                        Z[i,j] = Z[i,j] + c*M_tilde[i,j]*M_tilde[i,n]*X[i,n] + d*M_tilde[i,j]*M_tilde[n,j]*X[n,j] + e*M_tilde[i,j]*(M_tilde[j,n]*X[j,n] + M_tilde[n,i]*X[n,i])
                if Z[i,j] >= M_theta[i,j]:
                    W[i,j]=1
                else:
                    W[i,j]=0
            if (i%50 == 0):
                print("row={0} out of N={1}".format(i,N))
        print("\n")
        
    return(W)