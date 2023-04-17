# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 07:11:30 2023

@author: Nathan

Chapter 4 code snippets for dynamic distance ranging
"""

# Ajust grid
import math
import numpy as np

omega = 1.1
h = 1e-4
X = [0.0]

maxX = 6 * math.sqrt(maxT)

while X[-1] < maxX:
    X.append(X[-1] + h)
    h *= omega

n = len(X)

# Calculate Thomas coefficients
alpha = [0] * n
beta = [0] * n
gamma = [0] * n

for i in range(1, n-1):
    delX_m = X[i] - X[i-1]
    delX_p = X[i+1] - X[i]
    alpha[i] = -(2 * deltaT) / (delX_m * (delX_m + delX_p))
    gamma[i] = -(2 * deltaT) / (delX_p * (delX_m + delX_p))
    beta[i] = 1 - alpha[i] - gamma[i]
