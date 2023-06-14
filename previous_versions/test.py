# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 15:26:19 2023

@author: Nathan
"""

import pybamm
import numpy as np
import matplotlib.pyplot as plt

# Non-dimensionalise parameters
E0 = 0

E_start = 0.5
E_reverse = -0.5
t_reverse = E_start - E_reverse
v = 0.05

#creating time scale and non-dimensionalizing
Tmax_nd = ((abs(E_start - E_reverse) / v * 2))

#number of time steps
m = 20000

#length of time step, nondimensional
deltaT_nd = Tmax_nd / m
#kinetic constants
k0 = 10 #no units

Cdl = 10 #no units
Ru = 10 #no units

alpha = 0.5

# Input voltage protocol

b1 = (pybamm.t <= 10000)
b2 = (pybamm.t > 10000)
Edc_forward = -pybamm.t
Edc_backwards = pybamm.t - 2 * t_reverse
Eapp = E_start + \
    (pybamm.t <= t_reverse) * Edc_forward + \
    (pybamm.t > t_reverse) * Edc_backwards 

Eeff = Eapp                          


# Create model
model = pybamm.BaseModel()
O = pybamm.Variable("O")
R = pybamm.Variable("R")
Eeff = pybamm.Variable("E")
i = pybamm.Variable("i")

i_f = (k0 * (R) * np.e**((1-alpha) * (Eeff))) - (k0 * (O) * np.e**(-alpha * (Eeff)))

model.rhs = {O: i_f, 
             R: -i_f, 
             Eeff: (pybamm.t <= t_reverse) * Edc_forward + (pybamm.t > t_reverse) * Edc_backwards } # du/dt = -v

model.algebraic = {i: O - i} # i = i_f

model.initial_conditions = {O: pybamm.Scalar(1), 
                            R: pybamm.Scalar(0), 
                            i: pybamm.Scalar(0), 
                            Eeff: E_start}

model.variables = {"O": O, 
                   "R": R, 
                   "i": i, 
                   "E": Eeff}

# Discretise using default discretisation
disc = pybamm.Discretisation()
disc.process_model(model);

# Solve #################################
t_eval = np.linspace(0, Tmax_nd, m)
dae_solver = pybamm.CasadiSolver(mode="fast", atol=1e-8, rtol=1e-8)
solution = dae_solver.solve(model, t_eval)
#########################################

# Extract u and v
t_sol = solution.t
u = solution["O"](t_eval)
i = solution["i"](t_eval)
e = solution["E"](t_eval)

#QUICK PLOTS, IMPROVE# 
plt.cla()
plt.plot(e, i)
plt.xlabel("Eapp [V]")
plt.ylabel("current [A]")

plt.tight_layout()
plt.show()