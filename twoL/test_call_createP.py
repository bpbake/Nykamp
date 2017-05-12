# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:44:10 2017

@author: Brittany
"""

import numpy as np
from useP_to_create_W_withC import *
from create_P import *   
import matplotlib.pyplot as plt

N = 3000 #minimum of 1000
p_AVG = 50/N

#L_left = math.exp(np.random.uniform(1.2, 5))*(N/100)
#L_right = math.exp(np.random.uniform(1.2, 5))*(N/100)
L_left = 20000#*(N/100)
L_right = 0#1e-100

#alpha_recip = np.random.uniform(0, 0.3)
#alpha_conv = np.random.uniform(0, 0.3)
#alpha_div = np.random.uniform(0, 0.3)
#alpha_chain = np.random.uniform(-0.4, 0.3)
alpha_recip = .2
alpha_conv = .2
alpha_div = .2
alpha_chain = .2

P = create_P(N, L_left, L_right, p_AVG)
plt.matshow(P)
plt.show()


#W = create_W(N, P, alpha_recip, alpha_conv, alpha_div, alpha_chain)
#plt.matshow(W)