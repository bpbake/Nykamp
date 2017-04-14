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
p = 0.04
i = 1
print("N={0}".format(N))

#read in pickled W matrix
W_filename = "matrices\W_N{0}_p{1}_{2}.pickle".format(N,p,i)

with open(W_filename, 'rb') as wf:
    try:
        W = pickle.load(wf)
    except (EOFError):
        print("unpickling error")

# reading in pickled stats file:        
stat_filename = "matrices\Stats_W_N{0}_p{1}_{2}.pickle".format(N,p,i)
with open(stat_filename, 'rb') as sf:
    try:
        stats = pickle.load(sf)
    except (EOFError):
        print("unpickling error")

print(stats)

# add a new entry to the stats dictionary
stats["hi"] = 5 #key= "hi", value=5

print(stats)

# save updated stats dictionary as a picle file
stats_with_results_filename = "matrices\StatsWithSync_W_N{0}_p{1}_{2}.pickle".format(N,p,i)
with open(stats_with_results_filename, 'wb') as f:
    pickle.dump(stats, f)

# plot the W matrix
plt.matshow(W)
    
# generate statistics of matrix
#p_hat = np.sum(W)/(N*(N-1))
#alpha_recip_hat = (np.trace(np.matmul(W,W))/(N*(N-1)*math.pow(p_hat,2)))-1
#alpha_conv_hat = ((np.sum(np.matmul(np.transpose(W),W)) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_div_hat = ((np.sum(np.matmul(W,np.transpose(W))) - np.sum(W)) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1
#alpha_chain_hat = ((np.sum(np.matmul(W,W)) - np.trace(np.matmul(W,W))) / (N*(N-1)*(N-2)*math.pow(p_hat,2))) -1