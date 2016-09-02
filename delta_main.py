# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 20:28:52 2016

@author: Nicholas
"""
import delta_functions
import numpy as np
import matplotlib.pyplot as plt


# model parameters for the synthetic data set
N_STATES = 3 # number of states
MEMORY = 5 # memory
TIME = 50 # length of the time sequence
TIME_SERIES = np.arange(0, TIME + 1, 1) #Time vector
deltaT = 40.8 # time interval between consecutive transitions in [sec]
TRANSITION_RATES = np.array([[-0.01, 0.008, 0.001], [0.0085, -0.015, 0.015], [0.0015, 0.007, -0.016]])
TRANSITION_PROBS = delta_functions.rate_to_prob(TRANSITION_RATES, deltaT)
    
#TRANSITION_PROBS = np.array([[0.7, 0.3, 0.2], [0.2, 0.5, 0.4], [0.1, 0.2, 0.4]]) # transition matrix
EMISSIONS = np.array([0, 4, 8], dtype = float) # emission values
NOISE = 0.3 # gaussian noise
PI0 = np.array([0, 1, 0], dtype = float)

##Call function to generate synthetic data
fluo = np.append(np.array(0),delta_functions.synthetic(TIME, N_STATES, MEMORY, TRANSITION_PROBS, EMISSIONS, NOISE)) 
#print fluo        
##Call function to draw sample points from data        
SAMPLE = delta_functions.sample (fluo, TIME + 1, 0)
#print SAMPLE
plt.plot(SAMPLE[0,:], SAMPLE[1,:], TIME_SERIES, fluo, linewidth=2.0) 

##Call function to estimate change rates
deltas = delta_functions.delta(SAMPLE, MEMORY)

#print deltas[1,:]
plt.plot(deltas[0,:], deltas[1,:], linewidth=2.0) 

#plt.hist(deltas[1,:])