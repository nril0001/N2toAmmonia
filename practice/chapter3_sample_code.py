# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:36:10 2023

@author: Nathan
This is the code for a single electron transfer system
A + e- <--> B
"""

import math
import numpy as np

# Specify simulation parameters
theta_i = 20.0
theta_v = -20.0
sigma = 100.0
deltaX = 2e-4
deltaTheta = 0.02

# Calculate other parameters
deltaT = deltaTheta / sigma
maxT = 2 * abs(theta_v - theta_i) / sigma
maxX = 6 * math.sqrt(maxT)
n = int(maxX / deltaX)  # number of spacesteps
m = int(maxT / deltaT)  # number of timesteps

# Calculate Thomas coefficients
lambda_ = deltaT / (deltaX**2)
alpha = -lambda_
beta = 2.0*lambda_ + 1.0
gamma = -lambda_

# Create containers
g_mod = np.zeros(n)
C = np.ones(n)  # concentration profile

# Modify gamma coefficients
g_mod[0] = 0  # boundary condition
for i in range(1, n-1):
    g_mod[i] = gamma / (beta - g_mod[i-1] * alpha)

# Open file to output CV
CV = open("CV_Output.txt", "w")

# BEGIN SIMULATION
Theta = theta_i
for k in range(m):
    if k < m / 2:
        Theta -= deltaTheta
    else:
        Theta += deltaTheta
    
    # Forward sweep - create modified deltas
    C[0] = 1.0 / (1.0 + math.exp(-Theta))
    for i in range(1, n-1):
        C[i] = (C[i] - C[i-1]*alpha) / (beta - g_mod[i-1] * alpha)

    # Back Substitution
    C[n-1] = 1.0
    for i in range(n-2, -1, -1):
        C[i] = C[i] - g_mod[i] * C[i+1]

    # Output current
    flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*deltaX)
    CV.write(str(Theta) + "\t" + str(flux) + "\n")

# END SIMULATION
CV.close()