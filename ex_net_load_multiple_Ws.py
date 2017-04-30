# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:06:45 2017

@author: rhino
"""

try:
   import cPickle as pickle
except:
   import pickle

#import dill #pickle works fine


from analyze import analyze_autocor

try:
    del input
except:
    pass

input_orig = input

from brian2 import *

import numpy
numpy.set_printoptions(threshold=numpy.nan)

import matplotlib.pyplot as plt
import math

start_scope()

N = 1000 # Number of excitatory neurons
p_AVG =0.04

tau = 10*ms
Er = -60*mV
tauS = 5*ms
Ee = 0*mV

eqs= '''
dv/dt = (-v+Er)/tau: volt
'''

vthreshold = -55*mV
vreset = -65*mV
refract = 1*ms

transienttime = 500*ms
simulationtime = 1000*ms

#Neuron Groups
G = NeuronGroup(N, eqs, threshold='v>-55*mV', reset='v=-65*mV', refractory='refract', method='euler')
G.v='vreset+(vthreshold-vreset)*rand()'


ext_rate=300*Hz
ext_mag=1*mV

P = PoissonGroup(N, ext_rate)
Sp = Synapses(P,G, on_pre="v+=ext_mag")
Sp.connect(j='i')
#When source neuron fires a spike the target neuron will jump below value

j = 0.2*mV #Weight of neuron connection

S = Synapses(G, G,"w:volt",on_pre='v_post +=w') #Synapse from inhibitory neuron to inhibitory neuron
S.connect()

statemon = StateMonitor(G, 'v', record=0)
spikemon = SpikeMonitor(G)

PRM = PopulationRateMonitor(G)

store() # record state of simulation for future reference

#SII.connect(condition='i!=j', p=pii)

######Load in matrices one at a time
start_index = int(input_orig("enter a starting index: "))
end_index = int(input_orig("enter end index: "))

for w_index in range(start_index, end_index+1):
    
    restore() # set the state back to what it was when the store() command was called
    
    W_filename = "matrices\W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)
    #W_filename = "ws/W_N{0}_L{1}_c{2}_d{3}_ch{4}.txt".format(N,L,alpha_conv,alpha_div,alpha_chain)
    with open(W_filename, 'rb') as wf:
        try:
            W = pickle.load(wf) # load in W matrix
        except (EOFError):
            print("unpickling error")
    
    stats_filename = "matrices\Stats_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)
    with open(stats_filename, 'rb') as sf:
        try:
            stats = pickle.load(sf) # load in the stats for the W matrix (L, p_hat, alpha values, alpha_hat values)
        except (EOFError):
            print("unpickling error")
    
    S.w=W.transpose().flatten()*j
    
    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    PRM = PopulationRateMonitor(G)
    run(transienttime)
    
    
    if spikemon.num_spikes > (transienttime*N/refract*0.5): # if the number of spikes it too large, assume it's saturated
        print("\nnetwork saturated, skipping matrix {0}\n".format(w_index))
        stats['saturated'] = True # add to the stats dict
        result_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index) 
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) #pickle the new stats dict 
        continue # go to next matrix
    if spikemon.num_spikes < (2*N): # if the number of spikes is too small, we assume it's not spiking
        print("\nnetwork not spiking, skipping matrix {0}\n".format(w_index))
        stats['not spiking'] = True #add to the stats dict
        result_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index)  
        with open(result_filename, "wb") as rf:
            pickle.dump(stats, rf) # pickle the new stats file
        continue # go to next matrix
    print("\nnumber of spikes in transient: {0}\n".format(spikemon.num_spikes))
    
    
    for k,v in sorted(stats.items()):
        print(k+":{0}".format(v)) # print the statistics for W
    print("\n")
        
    
    # reset before run(simulationtime)
    statemon = StateMonitor(G, 'v', record=0)
    spikemon = SpikeMonitor(G)
    
    PRM = PopulationRateMonitor(G)    
    run(simulationtime)
    
    synchrony = analyze_autocor(PRM.rate)
    
    #plot the results of the simulation
    figure(figsize=(20,10))
    #subplot(122)
    #plot(statemon.t/ms, statemon.v[0])
    #xlabel('Time (ms)')
    #ylabel('v')
    
    subplot(211)
    plot(spikemon.t/ms,spikemon.i, '.k')
    xlabel('Time (ms)')
    ylabel('Neuron index')
    plt.tight_layout()
    
    subplot(212)
    plot(PRM.t/ms,PRM.smooth_rate(window='flat', width=0.5*ms)/Hz)
    
    print("\nExcitatory Synchrony = {0}\n".format(synchrony))
    
    # add to the stats dict
    stats['synchrony'] = synchrony
    stats['PRM rate'] = PRM.rate/hertz
    stats['PRM time'] = PRM.t/ms
    stats['spikemon times'] = spikemon.t/ms
    stats['spikemon index'] = spikemon.i/1
    
    result_filename = "matrices\Results_W_N{0}_p{1}_{2}.pickle".format(N,p_AVG,w_index) 
    with open(result_filename, "wb") as rf:
       pickle.dump(stats, rf) # pickle the new stats dict