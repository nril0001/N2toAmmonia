# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 11:14:40 2023

@author: Nathan
"""
## Main function
import catalyticmodel06 as cm
import time
import matplotlib.pylab as plt
import os
from datetime import date
import numpy as np
from scipy.optimize import curve_fit


def main(data, x):
    area = 1
    radius = np.sqrt(area/np.pi)
    k0 = x
    Ru = 0
    Cdl = 0
    Gamma = 1
    atol = 1e-7
    rtol = 1e-7
    t_steps = 2**(10)
    x_steps = 750
    solver = "Casadi"
    # solver = "Scikits"
    
    const_parameters = {
        "Faraday Constant [C mol-1]": 96485.3328959,
        "Gas constant [J K-1 mol-1]": 8.314459848,
        "Temperature [K]": 298.2,
        "Far-field concentration of S(soln) [mol cm-3]": 1,
        "Far-field concentration of P(soln) [mol cm-3]": 0,
        "Surface coverage of P [mol cm-2]": 0,
        "Electrode Coverage [mol cm-2]": Gamma,
        "Diffusion Coefficient of S [cm2 s-1]": 1,
        "Diffusion Coefficient of P [cm2 s-1]": 1,
        "Diffusion Layer Thickness [cm]": 1,
        "Electrode Area [cm2]": 1,
        "Electrode Radius [cm]": radius,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": 0.0,
        "Scan Rate [V s-1]": 2.569689e-04,
        "Uncompensated Resistance [Ohm]": Ru,
        "Capacitance [F]": Cdl, #1e-8,
    }
        
    #conditions that will change often over the course of testing
    input_parameters = {
        "G": 0,
        "G'": 0,
        "Reversible Potential 1 [V]": 0,
        "Electrosorption Rate [mol-1 cm3 s-1]": k0,
        "Symmetry factor [non-dim]": 0.5,
        
    }    
        
    #setting main model to reference CatalyticModel class
    cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps, x_steps, solver)
    current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
    xss = x_steps
    tss = t_steps
    while len(E_nd) == tss/2:
        xss = xss + 2
        cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, tss, xss, solver)
        current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        
    ##redimensionalizing here for now. Messy to do in main, move later
    I_d = current 
    E_d = E_nd / cmodel._E_0
    # print("complete in time: " + str((time.time()-ti)/60) + " minutes") 
    return I_d[:1024]

if __name__ =='__main__':
    #list of options to pass into the model
    seioptions = ()
    
    #retrieve date
    today = date.today()
    folder = today.strftime("%Y%m%d")
    
    #folder pathway - needs to be set to read all txt files that want analysed
    #use same naming conventions for text files as files in GammaTemp Variation
    output = "output/"+folder+"/"
    
    if os.path.isdir(output) == 0:
        os.mkdir(output, 0o666)
        os.mkdir(output+"SC", 0o666)
        os.mkdir(output+"CV", 0o666)
    
    files = ['output/20230928/test.txt']
    
    for i in files:
        # Open the input file and read lines
        with open(i, 'r') as file:
            lines = file.readlines()
        
        print(i)
        
        volt = []
        curr = []
        row = []
        for row in lines:
            row = row.split(" ")
            # print(row[1])
            # print(float(row[1]))
            volt.append(float(row[0]))
            curr.append(float(row[1]))
    
    #Timer to measure performance
    ti = time.time()
    # e, i = main(1, 100000000)
    
    # plt.cla()
    plt.plot(volt, curr)  
    # plt.xlabel("Eapp [V]")
    # plt.ylabel("Current [A]")
    # # plt.legend(loc='upper left')
    plt.grid()
    # plt.savefig(output+"fitting.png", dpi=1000)
    
    params = curve_fit(main, volt, curr, p0=[0.0000001])
    print(params)
    
    
    
    
