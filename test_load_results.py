# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 00:02:45 2017

@author: rhino
"""
N = 1000
p_AVG = .04
w_index = 1


results_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)
with open(results_filename, 'rb') as rf:
    try:
        results = pickle.load(rf) 
    except (EOFError):
        print("unpickling error")
            
print("\n{0}".format(results))