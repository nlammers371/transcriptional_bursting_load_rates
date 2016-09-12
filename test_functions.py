# -*- coding: utf-8 -*-
"""
Created on Sun Sep 04 10:01:31 2016

@author: Nicholas
"""

import numpy as np
import matplotlib.pyplot as plt
import math 

M = 1000.
N = np.arange(1,int(M)+1, dtype=float)
P = np.zeros(int(M), dtype=float)
F = np.zeros(int(M), dtype=float)

for i in range(0,int(M)):
    f = N[i]/M
    F[i] = f    
    logP = -M*math.log(2) - M*f*math.log(f*((1-f)**(1/f-1)))
    P[i] = math.exp(logP)

plt.plot(F,P)
    