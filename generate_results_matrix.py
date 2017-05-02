# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 00:02:45 2017

@author: rhino
"""
N = 1000
p_AVG = .04

start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

result_matrix = [[0 for col in range(13)] for row in range(400)]
                                                  
for w_index in range(start_index, end_index+1):
   
    results_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)
    with open(results_filename, 'rb') as rf:
        try:
            results = pickle.load(rf) 
        except (EOFError):
            print("unpickling error")
                
   # print("\n{0}".format(results))
    result_matrix[w_index][0] = w_index
    result_matrix[w_index][1] = results['synchrony']
    result_matrix[w_index][2] = results['L']
    result_matrix[w_index][3] = results['alpha_chain_hat']
    result_matrix[w_index][4] = results['alpha_conv_hat']
    result_matrix[w_index][5] = results['alpha_div_hat']
    result_matrix[w_index][6]= results['alpha_recip_hat']
    result_matrix[w_index][7] = results['alpha_chain']
    result_matrix[w_index][8] = results['alpha_conv']
    result_matrix[w_index][9] = results['alpha_div']
    result_matrix[w_index][10] = results['alpha_recip']
    result_matrix[w_index][11] = results['p_hat']
    result_matrix[w_index][12] = results['p_AVG']


np.savetxt('result_matrix.txt', result_matrix, delimiter=',') 
    
    


