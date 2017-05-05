# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:41:12 2017

@author: Brittany
"""

import numpy as np
from useP_to_create_W_withC import *
from create_P import *
try:
   import cPickle as pickle
except:
   import pickle

import sys

N = 1000 #minimum of 1000
p_AVG = 0.04

if len(sys.argv) >= 3:
   start_index = int(sys.argv[1])
   num_matrices = int(sys.argv[2])
else:
   start_index = int(input("enter a starting index: "))
   num_matrices = int(input("how many matrices to generate? "))

for i in range(start_index, start_index+num_matrices): #so i=start_index, start_index+1,start_index+2,...,start_index+num_matrices-1
    np.random.seed(i)
    
    trying = True
    while trying:  
        try:
            print("\nmaking matrix {0}".format(i))
            
            #generate Ls, alphas
            L_left = math.exp(np.random.uniform(1.2, 5))*(N/100)
            L_right = math.exp(np.random.uniform(1.2, 5))*(N/100)
            alpha_recip = np.random.uniform(0, 0.3)
            alpha_conv = np.random.uniform(0, 0.3)
            alpha_div = np.random.uniform(0, 0.3)
            alpha_chain = np.random.uniform(-0.4, 0.3)

            P = create_P(N, L_left, L_right, p_AVG)
            
            #call other program to create W (and return W)
            W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
            
            #save the W
            W_filename = "matrices/W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,i)
            with open(W_filename, 'wb') as fp:
                pickle.dump(W, fp)
                
            # generate statistics of W
            p_hat = np.sum(W)/(N*(N-1))
            alpha_recip_hat = (np.trace(np.matmul(W,W))/(N*(N-1)*math.pow(p_hat,2)))-1
            alpha_conv_hat = ((np.sum(np.matmul(np.transpose(W),W)) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_div_hat = ((np.sum(np.matmul(W,np.transpose(W))) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
            alpha_chain_hat = ((np.sum(np.matmul(W,W)) - np.trace(np.matmul(W,W))) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
                
            #create dictionary of stats and save
            stats = dict([('L_left', L_left), ('L_right', L_right), ('p_AVG', p_AVG), ('alpha_recip', alpha_recip), ('alpha_conv', alpha_conv), ('alpha_div', alpha_div), ('alpha_chain', alpha_chain), ('p_hat', p_hat), ('alpha_recip_hat', alpha_recip_hat), ('alpha_conv_hat', alpha_conv_hat), ('alpha_div_hat', alpha_div_hat), ('alpha_chain_hat', alpha_chain_hat)])
            
            stat_filename = "matrices/Stats_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,i) #pickle the dictionary of stats for each W
            with open(stat_filename, "wb") as f:
                pickle.dump(stats, f)
                
            trying = False
                
        except Exception as error:
            print(error)
            continue
        
        
