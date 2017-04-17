# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

try:
   import cPickle as pickle
except:
   import pickle


from analyze import analyze_autocor

from brian2 import *

import numpy
numpy.set_printoptions(threshold=numpy.nan)

import matplotlib.pyplot as plt
import math

start_scope()

N = 1000 #Number of excitatory neurons

tau = 10*ms
Er = -60*mV
tauS = 5*ms
Ee = 0*mV
#eqs= '''
#dv/dt = (-v+Er)/tauI : volt
#'''

eqs= '''
dv/dt = (-v+Er)/tau: volt
'''

vthreshold = -55*mV
vreset = -65*mV
refract = 10*ms
transienttime = 100*ms
simulationtime = 400*ms

#Neuron Groups

G = NeuronGroup(N, eqs, threshold='v>-55*mV', reset='v=-65*mV', refractory='refract', method='euler')

G.v='vreset+(vthreshold-vreset)*rand()'


ext_rate=300*Hz
ext_mag=1*mV

P = PoissonGroup(N, ext_rate)
S = Synapses(P,G, on_pre="v+=ext_mag")
S.connect(j='i')

#When source neuron fires a spike the target neuron will jump below value

j = 0.1*mV

#Weight of neuron connection

p_avg =0.04

delta=2*ms

S = Synapses(G, G,"w:volt",on_pre='v_post +=w') #Synapse from inhibitory neuron to inhibitory neuron
S.connect()
#SII.connect(condition='i!=j', p=pii)

######Load in matrix
#L = math.inf
#alpha_conv=0.3
#alpha_div=0.3
#alpha_chain = 0.3
start_index = int(input("enter a starting index: "))
end_index = int(input("enter end index"))

for i in range(start_index, end_index):
    
    filename = "W_N{0}_p{1}_{2}.pickle".format(N,p_avg,i)
    #"ws/W_N{0}_L{1}_c{2}_d{3}_ch{4}.pickle".format(NI,L,alpha_conv,alpha_div,alpha_chain)
    #filename = "ws/W_N{0}_L{1}_c{2}_d{3}_ch{4}.txt".format(NI,L,alpha_conv,alpha_div,alpha_chain)
    
    with open(filename, 'rb') as fp:
        try:
            W = pickle.load(fp)
            #W=numpy.loadtxt(filename)
        except (EOFError):
            print("unpickling error")
    
    stat_filename = "Stats_W_N{0}_p{1}_{2}.pickle".format(N,p_avg,i)
    
    with open(stats_filename, 'rb') as fp:
        try:
            stats = pickle.load(fp)
            #W=numpy.loadtxt(filename)
        except (EOFError):
            print("unpickling error")
    
            
    S.w=W.transpose().flatten()*j
    

    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    PRMi = PopulationRateMonitor(G)
    
    # run(transienttime)
    
    # if spikemon.num_spikes > (transienttime*N/refract*0.5):
    #     print("\nnetwork saturated, skipping matrix\n")
    
    # if spikemon.num_spikes < (2*N):
    #     print("\nnetwork not spiking, skipping matrix\n")
        
    run(simulationtime)
    
    synchrony = analyze_autocor(PRMi.rate)
    
    
    figure(figsize=(8,5))
    #subplot(122)
    #plot(statemon.t/ms, statemon.v[0])
    #xlabel('Time (ms)')
    #ylabel('v')
    
    #subplot(121)
    plot(spikemon.t/ms,spikemon.i, '.k')
    xlabel('Time (ms)')
    ylabel('Neuron index')
    plt.tight_layout()
    
    print("\nExcitatory Synchrony = {}".format(synchrony))
   
    results = dict([('stats',stats),('synchrony', synchrony),('Population Rate Monitor',PRMi)])   
    result_filename = "Results_W_N{0}_p{1}_{2}.pickle".format(N,p,i) 
    with open(result_filename, "wb") as f:
        pickle.dump(results, f)    
   