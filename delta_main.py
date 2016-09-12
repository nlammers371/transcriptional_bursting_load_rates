# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 20:28:52 2016

@author: Nicholas
"""
import delta_synth_functions as d_synth
import numpy as np
import matplotlib.pyplot as plt
#import scipy.interpolate as interp
#import scipy.signal as sig

# model parameters for the synthetic data set
N_STATES = 3 # number of states
MEMORY = 5 # memory (in terms of # states)
V_SCALE = 1 #number of time steps per state step (assume state dwell ~1min)
TIME = 600 # time steps in sequence (assume time res ~10s and set series length ~20min)
TIME_SERIES = np.arange(0, TIME+V_SCALE, 1) #Time vector
TRANSITION_PROBS = np.array([[.3, .3, .5], \
                             [.2, .3, .1], \
                             [.5, .4, .4]])
#TRANSITION_PROBS = np.array([[0.7, 0.3, 0.2], [0.2, 0.5, 0.4], [0.1, 0.2, 0.4]]) # transition matrix
EMISSIONS = np.array([0, 4, 8], dtype = float) # emission values
V_NOISE = 0.4 # gaussian noise added to underlying states
#noise signal vector. R1 gives signal period in time steps, R2 gives magnitude relative to max(c)
EXP_NOISE = np.array([[2.], [.001]], dtype=float)  
PI0 = np.array([0, 1, 0], dtype = float)

##Call function to generate synthetic data
#synthetic (t_steps, v_scale, K, w, A, v, v_noise, exp_noise)
fluo =d_synth.synthetic(TIME, V_SCALE, N_STATES, MEMORY, TRANSITION_PROBS, EMISSIONS, V_NOISE, EXP_NOISE) 
fluo_zero = np.append(np.zeros(V_SCALE),fluo)
'''
    Try simplest approach first: take differences directly.
    Assess results using line plots + histogram
    
    Note: unassisted differentiation method begins to fail as noise approaches
          1/10 the separation of states.
'''
##take differences
diff = np.append(np.diff(fluo_zero),0)

#use knowledge of memory to disentangle loading rate from unloading rate
mem = MEMORY*V_SCALE

l_rates = np.zeros(len(diff), dtype = float)
l_rates[0:mem] = diff[0:mem]
#print l_rates
#iteratively calcualte loading rates
for i in range(mem,len(l_rates)):
    u = l_rates[i-mem]
    l_rates[i] = max(u,0) + diff[i] #assume l_rates < 0 are erroneous

#print fit
#plt.plot(TIME_SERIES, fluo_zero, TIME_SERIES, l_rates, TIME_SERIES, diff, linewidth=2.0) 
#histogram
plt.hist(l_rates)