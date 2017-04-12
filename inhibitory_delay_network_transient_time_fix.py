# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 15:51:29 2017

@author: nykamp
"""

from brian2 import *
import numpy
try:
   import cPickle as pickle
except:
   import pickle
   
numpy.set_printoptions(threshold=numpy.nan)

start_scope()

NI = 1000 #Number of inhibitory neurons #should match W matrix value

tauI = 10*ms
Er = -60*mV

muext = 25*mV
sigmaext = 1*mV
tau = 20*ms

eqs = """
dv/dt = (-v+muext + sigmaext * sqrt(tau) * xi)/tau : volt
"""
#old equation
#eqs= '''
#dv/dt = (-v+Er)/tauI : volt
#'''

vthreshold = 20*mV
vreset = 10*mV
refract = 2*ms
transienttime = 100*ms
simulationtime = 100*ms

#Neuron Groups

GI = NeuronGroup(NI, eqs, threshold='v>vthreshold',
                    reset='v=vreset', refractory=refract, method='euler')

#GI.v='vreset+(vthreshold-vreset)*rand()'
GI.v=vreset

ext_rate=300*Hz
ext_mag=1*mV
#
#P = PoissonGroup(NI, ext_rate)
#S = Synapses(P,GI, on_pre="v+=ext_mag")
#S.connect(j='i')

#When source neuron fires a spike the target neuron will jump below value

jii = .7*mV
print("jii={0}".format(jii))

#Weight of neuron connection

pii =0.1

delta=5*ms

SII = Synapses(GI, GI,"w:volt",on_pre='v +=-w',delay=delta) #Synapse from inhibitory neuron to inhibitory neuron
SII.connect()
#SII.connect(condition='i!=j', p=pii)

######Load in matrix
N = 1000
L = 100
alpha_chain =0.39

filename = "W_N{0}_L{1}_Alpha{2}.pickle".format(N,L,alpha_chain)

with open(filename, 'rb') as fp:
    try:
        W = pickle.load(fp)
    except (EOFError):
        print("unpickling error")
        
SII.w=W.transpose().flatten()*jii

statemon = StateMonitor(GI, 'v', record=0)
spikemon = SpikeMonitor(GI)

PRMi = PopulationRateMonitor(GI)

run(transienttime)
print("\nnumber of spikes: %s\n" % spikemon.num_spikes)
#print("\npopulation rate: %s\n" %PRMi.rate)

if spikemon.num_spikes > (transienttime*NI/refract*0.5):
    print("\nnetwork saturated, skipping matrix\n")
    NewPRM = PopulationRateMonitor(GI)
    
if spikemon.num_spikes < (1*NI):
    print("\nnetwork not spiking, skipping matrix\n")
    NewPRM = PopulationRateMonitor(GI)

run(simulationtime)

if spikemon.num_spikes > (transienttime*NI/refract*0.5) or spikemon.num_spikes < (1*NI):
    i_synchrony = analyze_autocor(NewPRM.rate)
else:
    i_synchrony = analyze_autocor(PRMi.rate)


#figure(figsize=(12,4))
#subplot(122)
#plot(statemon.t/ms, statemon.v[0])
#xlabel('Time (ms)')
#ylabel('v')

subplot(211)
plot(spikemon.t/ms,spikemon.i, '.k')
xlabel('Time (ms)')
ylabel('Neuron index')

subplot(212)
plot(PRMi.t/ms,PRMi.smooth_rate(window='flat', width=0.5*ms)/Hz)

#subplot(213)
#plot(NewPRM.t/ms,NewPRM.smooth_rate(window='flat', width=0.5*ms)/Hz)

print("\nInhibitory Synchrony = {}".format(i_synchrony))

