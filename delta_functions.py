# -*- coding: utf-8 -*-
"""
Created on Thu Sep 01 20:06:03 2016

@author: Nicholas
"""
#import matplotlib.pyplot as plt
import numpy as np
#import itertools
#import timeit


"Synthetic Data Generator Adapted from Vahe's Code"

def synthetic (seq_length, K, w, A, v, noise):
    """ Generates a fluorescence sequence using the given model parameters.

    Args:
        seq_length: length of the sequence
        K: number of naive states
        w: memory
        A: transition matrix
        v: emission values
        noise: Gaussian noise

    Returns:
        Fluorescence sequence of length seq_length, generated randomly using
        the model parameters. Gaussian noise is added at all time points.
    """
    naive_states = np.zeros(seq_length, dtype=int)
    # assigns naive state 1 as the first state in the sequence
    naive_states[0] = 1

    # generates a sequence of naive states using the transition probabilities
    for t in xrange(1,seq_length):
        naive_states[t] = np.random.choice(K, 1, p=A[:,naive_states[t-1]])

    fluo_values = np.zeros(seq_length)

    # calculates the compound flurescence using the naive states
    for t in xrange(w):
        state_ls = naive_states[xrange(t+1)]
        fluo_values[t] = sum([v[i] for i in state_ls])

    for t in xrange(w, seq_length):
        state_ls = naive_states[xrange(t-w+1, t+1)]
        fluo_values[t] = sum([v[i] for i in state_ls])

    # adds a gaussian noise
    fluo_values += np.random.normal(loc=0.0, scale=noise, size=seq_length)

    return fluo_values

"function borrowed from Vahe to convert transition rates to transition probs"

def rate_to_prob (R, deltaT):
    """ Calculates the transition probability matrix using the transition rates

    Args:
        R: transition rate matrix
        deltaT: time interval between consecutive transitions

    Returns:
        Transition probability matrix.
    """
    RT = R*deltaT
    RT_eigen, RT_right_vec = np.linalg.eig(RT)
    RT_left_vec = np.linalg.inv(RT_right_vec)
    RT_eigen_exp_diag = np.diag(np.exp(RT_eigen))
    A = reduce(np.dot, [RT_right_vec, RT_eigen_exp_diag, RT_left_vec])
    return A

"""
    function to sample points from data
    
    trace (float array): Fluorescence Data 
    points (int):  # sample points to draw
    rand (binary):  If 1 randomly draw sample points without replacement
                   Else, points uniformly selected
"""

def  sample (trace, points, rand):
    l = len(trace)
    #generate array of indices
    #if random option off, use regular intervals
    if rand == 0:
        indices = np.arange(0, l, l//points, dtype=int)
        times = np.arange(0, l, l//points, dtype=float)
    #if rand option on, randomly select points without replacement
    if rand == 1:
    #select random indices without replacement
        choices = np.arange(0, l, 1)
        indices = np.sort(np.random.choice(choices, points, replace=False))
        times = np.array(np.sort(np.random.choice(choices, points, replace=False)), dtype=float)
    #select points
    samp = trace[indices]
    #define sample matrix to contain time points and sample points
    sample = np.zeros([2,points])
    sample[0, :] = times
    sample[1, :] = samp
    return sample
    
"""
    function to estimate change rates from sample points
    
    samp (float array): array of sample points and times 
    length (int): desired span for change estimate (defaualted to 2)
    window (int): dwell time of system
"""
def delta (samp, window, length=2):
    #subtract nth array element from (n+2)th element
    diffs = np.diff(samp)
    #divide fl change by delta t
    diffs = np.append(np.array(0),np.divide(diffs[1,:], diffs[0,:]))
    #define time array.
    t_vec = np.append(np.array(samp[0,0]-samp[0,1])/2,np.add(samp[0,1:],samp[0,0:(-1)])/2)
    #assume regular sampling interval to calc #steps per window
    w = int(window/(t_vec[1]-t_vec[0]))    
    #for deltas ranging from length of window to end of sequence, add 
    #rate of change at t = T-w if greater than 0. This accounts for fact that
    #these polymerases are now terminating transcriptions    
    fix_deltas = diffs
    for i in range(w,len(t_vec)):
        fix_deltas[i] = diffs[i] + max(diffs[i-w],0)
    #define matrix to hold outputs
    deltas = np.zeros([2, len(fix_deltas)])
    deltas[0,:] = t_vec
    deltas[1,:] = fix_deltas 
    return deltas
