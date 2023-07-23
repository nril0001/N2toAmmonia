## Main function
import catalyticmodel05_3 as cm
import time
import re
import matplotlib.pylab as plt
import os
from datetime import date
import numpy as np


def main():
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
    
    files = [['DigiElech/2023-06-26 solution/CV/CV_k0_10000_R_800_Cd_1e-9_Ds_2e-5_Dp_1.5e-5.txt', 1e-3, 1e-5]]
    
    # retrieve variables from file pathway
    for t in files:
        variables = re.findall(
            "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", t[0])
        for l in range(len(variables)):
            variables[l] = float(variables[l])
            t.append(variables[l])
    
    for i in files:
        print(i)
        
        volt = []
        curr = []
        row = []
        count = 0 
        f = open(i[0],'r')
        for row in f:
            count += 1
            row = row.split("\t")
            if count < 3:
                continue
            else:
                volt.append(float(row[0]))
                curr.append(float(row[1]))
    
    curr = np.array(curr)
    
    k0 = 100
    Ru = 800
    Cdl = 1e-9
    atol = 1e-12
    rtol = 1e-8
    t_steps = 2**14
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
    #setting solved answers to ones usable here
    # current, E_nd, O_nd, R_nd, S_nd, P_nd, T_nd = cmodel.simulate(input_parameters)
    current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
    # current, E_nd, T_nd = cmodel.simulate(input_parameters)
    # simulating analytical solution
    #I_ana_nd = amodel.simulate(E_nd)
    ##redimensionalizing here for now. Messy to do in main, move later
    I_d = current * cmodel._I_0
    E_d = E_nd * cmodel._E_0

    #I_ana_d = I_ana_nd *cmodel._I_0
    #print(k0[0])
    # print(len(curr))
    
    #QUICK PLOTS, IMPROVE# 
    
    # Plot current
    plt.cla()
    plt.plot(volt, np.array(-curr), color = 'orange', linestyle = 'dashdot', label = 'Digielch')
    plt.plot(E_d[1:], I_d[1:], color = 'Red', label = 'Pybamm', linestyle='dashed')
    plt.xlabel("Eapp [V]")
    plt.ylabel("Current [A]")
    plt.legend()
    plt.grid()
    plt.savefig(output+"CV_cat04_dim.png", dpi=600)
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, O_nd, R_nd))))
    
    # plt.cla()
    # plt.plot(T_nd, capa, color = 'Red', label = 'Pybamm', linestyle='dashed')
    # plt.xlabel("Time [s]")
    # plt.ylabel("Eapp [V]")
    # plt.legend()
    # plt.grid()
    
    #Plot concentration profiles
    # plt.cla()
    # plt.plot(volt[:398], np.array(curr[:398])/curr[0], color = 'orange', linestyle = 'dashdot', label = 'Digielch - S')
    # plt.plot(volt[399:], np.array(curr[399:])/curr[0], color = 'green', linestyle = 'dashdot', label = 'Digielch - P')
    # plt.plot(E_d[1:], O_nd[1:], color = 'Red', label = 'S', linestyle='dashed')
    # plt.plot(E_d[1:], R_nd[1:], color = 'Blue', label = 'P', linestyle='dashed')
    # plt.xlabel("Eapp [V]")
    # plt.ylabel("Surface Coverage [non-dim]")
    # plt.legend()
    # plt.grid()
    #plt.savefig(output+"CV_cat04_dim.png")
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, O_nd, R_nd))))
    
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
