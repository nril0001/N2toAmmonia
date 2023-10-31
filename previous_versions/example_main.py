## Main function
import model07 as cm
import time
import re
import matplotlib.pylab as plt
import os
from datetime import date
import numpy as np


def main():
    #list of options to pass into the model
    seioptions = ()
    
    k0 = 100
    Ru = 800000
    Cdl = 100
    atol = 1e-5
    rtol = 1e-5
    t_steps = 2**12
    x_steps = 100
    
    #constants that can vary, but generally won't change expt to expt
    const_parameters = {
        "Faraday Constant [C mol-1]": 96485.3328959,
        "Gas constant [J K-1 mol-1]": 8.314459848,
        "Far-field concentration of S(soln) [mol cm-3]": 1e-6,
        "Far-field concentration of P(soln) [mol cm-3]": 0,
        "Diffusion Coefficient of S [cm2 s-1]": 2e-5,
        "Diffusion Coefficient of P [cm2 s-1]": 1.5e-5,
        "Electrode Area [cm2]": 1,
        "Temperature [K]": 298.2,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": .0,
        "Scan Rate [V s-1]": 0.05,
        "Electrode Coverage [mol cm-2]": 1e-12,
    }
    
    #conditions that will change often over the course of testing
    input_parameters = {
        "Reversible Potential [V]": 0,
        "Redox Rate [s-1]": k0,
        "Catalytic Rate For [cm2 mol-l s-1]": 1000*1000,
        "Catalytic Rate Back [cm2 mol-l s-1]": 1000*1000,
        "Symmetry factor [non-dim]": 0.5,
        "Uncompensated Resistance [Ohm]": Ru,
        "Capacitance [F]": Cdl, #1e-8,
    }
    
    #setting main model to reference CatalyticModel class
    cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps, x_steps)
    #e, c = aa.AnalyticalModel(100)
    current, E_nd, T_nd = cmodel.simulate(input_parameters)
    # simulating analytical solution
    #I_ana_nd = amodel.simulate(E_nd)
    ##redimensionalizing here for now. Messy to do in main, move later
    I_d = current * cmodel._I_0
    E_d = E_nd * cmodel._E_0

    
    #QUICK PLOTS, IMPROVE# 
    
    plt.cla()
    plt.plot(T_nd, E_d, color = 'Red', label = 'Pybamm', linestyle='dashed')
    plt.xlabel("Eapp [V]")
    plt.ylabel("Current [A]")
    plt.legend()
    plt.grid()
    
    
    return

if __name__ =='__main__':
    #Timer to measure performance
    ti = time.time()

    main()
    

    #literally just to test that that main is working properly (delete later)
    print("complete in time: " + str((time.time()-ti)/60) + " minutes")

    #TODO: add multiple fxns for running models at the same time
    #TODO: improve plotting for pybamm based models: investigate using pybamm's plotting features??
    #TODO: LOW PRIORITY make compatible with interactive plottersz
