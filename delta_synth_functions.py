# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 20:06:03 2016

@author: Nicholas
"""
#import matplotlib.pyplot as plt
import numpy as np
#import itertools
#import timeit


'''
    Code to generate synthetic data for validation 
    of FL trace differentiation model
'''

def synthetic (t_steps, v_scale, K, w, A, v, v_noise, exp_noise):
    """ Generates a fluorescence sequence using the given model parameters.

    Args:
        t_steps (int): total time steps in sequence 
        v_scale (int): number time steps making up each promoter state in seq
        K (int): number of naive states
        w (int): memory
        A (float matrix): transition matrix
        v (float array): emission values
        v_noise (float): STD of Gaussian noise added to "v" states
        exp_noise (float matrix): specifies scale and frequencies of additional
                                  noise signals to add

    Returns:
        Fluorescence sequence t_samp sampling steps, generated randomly using
        the model parameters. Gaussian noise is added to all loading rates, vs.
    """
    #determine number of promoter states based on t_steps and v_scale
    v_steps = np.divide(t_steps,v_scale)
    
    naive_states = np.zeros(v_steps, dtype=int)
    # assigns naive state 1 as the first state in the sequence
    naive_states[0] = 1

    # generates a sequence of naive states using the transition probabilities
    for t in xrange(1,v_steps):
        naive_states[t] = np.random.choice(K, 1, p=A[:,naive_states[t-1]])

    fluo_values = np.zeros(v_steps)

    # calculates the compound flurescence using the naive states
    for t in xrange(w):
        state_ls = naive_states[xrange(t+1)]
        fluo_values[t] = sum([v[i] for i in state_ls])

    for t in xrange(w, v_steps):
        state_ls = naive_states[xrange(t-w+1, t+1)]
        fluo_values[t] = sum([v[i] for i in state_ls])

    # adds gaussian noise
    fluo_values += np.random.normal(loc=0.0, scale=v_noise, size=v_steps)

    #generate array with length equal to number of desired time steps
    fluo_samp = np.zeros(t_steps, dtype=float)

    #fill with fluo_values
    for s in range(0,v_steps):
        fluo_samp[v_scale*s:v_scale*(s+1)] = fluo_values[s]
    
    #now add noise signals specified by exp_noise input 
    #restricted to Gaussian noise for now
    n_signals = len(exp_noise[1,:])   
    
    e_noise = np.zeros((t_steps), dtype=float)
    
    for i in range(0,n_signals):
        #determine # noise signals using scale paramterters stored in
        #first row of exp_nois array
        noise_scale = int(exp_noise[0,i])
        n_steps = t_steps//noise_scale
        factor = exp_noise[1,i]
        #generate noise vector
        rand = np.random.normal(loc=0.0, scale=factor, size=n_steps)
        
        #now add each signal to e_noise vector
        for k in range(0,n_steps):
            
            e_noise[k*noise_scale:(k+1)*(noise_scale)] += rand[k]
    
    #add noise to fluo_samp
    fluo_noise = fluo_samp + e_noise
            
    return fluo_noise