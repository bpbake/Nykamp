# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:24:58 2017

@author: Brittany
"""
try:
   import cPickle as pickle
except:
   import pickle
   
import matplotlib.pyplot as plt
import numpy as np
import math

N = 1000
print("N={0}".format(N))

W_filename = "W_N{0}_{1}.pickle".format(N,1)
print("filename: {0} \n".format(filename))

with open(W_filename, 'rb') as fp:
    try:
        W = pickle.load(fp)
    except (EOFError):
        print("unpickling error")
        
stat_filename = "Stats_W_N{0}_{1}.pickle".format(N,1)

with open(stat_filename, 'rb') as fs:
    try:
        stats = pickle.load(fs)
    except (EOFError):
        print("unpickling error")

print(stats)

plt.matshow(W)
    
# generate statistics of matrix
#p_hat = np.sum(W)/(N*(N-1))
#alpha_recip_hat = (np.trace(np.matmul(W,W))/(N*(N-1)*math.pow(p_hat,2)))-1
#alpha_conv_hat = ((np.sum(np.matmul(np.transpose(W),W)) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_div_hat = ((np.sum(np.matmul(W,np.transpose(W))) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_chain_hat = ((np.sum(np.matmul(W,W)) - np.trace(np.matmul(W,W))) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1