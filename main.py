## Main function
import catalyticmodel06 as cm
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
    
    files = [['DigiElech/2023-08-23 diffusion only comparison/CV/0 Ru/CV_k0_1e3_R_0_Cd_0_Ds_1e-5_v_4e-1.txt', 1e-3, 1e-5]]
    
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
    area = 1
    radius = np.sqrt(area/np.pi)
    k0 = [1]
    kf = 1e-10
    kb = 1e-10
    Ru = 1
    Cdl = 0
    atol = 1e-7
    rtol = 1e-7
    t_steps = [2**(10)]
    x_steps = [750]
    # solver = "Casadi"
    solver = "Scikits"
    
    DS_d = 1
    DP_d = 1
    DY_d = 1
    CS_d = 1
    Gamma = 1
    F = 96485.3328959
    R = 8.314459848
    T = 298.2
    # r = np.sqrt(area/np.pi)  
    r = 100
    G = 0
    G_ = 0
    
    # Create scaling factors for non-dimensionalisation
    T_0 = DS_d/r**2 # time, units in s-1
    X_0 = 1/r # distance, units in cm-1
    C_0 = 1/CS_d # concentration, units in mol cm-3 
    G_0 = 1/Gamma # surface concentration units in mol-1 cm2
    K_0 = (Gamma*r)/DS_d # adsorption rate constant, units in mol-3 cm s
    D_0 = 1/DS_d # diffusion, units in s cm-2
    V_0 = (r**2/DS_d)*(F / (R * T)) # scan rate, units in s V-1
    B_0 = (r*CS_d)/Gamma #Saturation coefficient
    E_0 = F / (R * T) # potential, units in V-1
    I_0 = 1/(np.pi*r*F*DS_d*CS_d) # current, units in A-1
    
    E_ds = []
    I_ds = []
    Z_ds = []
    
    for o in k0:   
        print(o*K_0)
        # print(B_0)
        # print(2.569689992599293e-04*V_0)
        #constants that can vary, but generally won't change expt to expt
        const_parameters = {
            "Faraday Constant [C mol-1]": F,
            "Gas constant [J K-1 mol-1]": R,
            "Temperature [K]": T,
            "Far-field concentration of S(soln) [mol cm-3]": CS_d,
            "Far-field concentration of P(soln) [mol cm-3]": 0,
            "Far-field concentration of Y(soln) [mol cm-3]": CS_d,
            "Surface coverage of P [mol cm-2]": 0,
            "Surface coverage of X [mol cm-2]": 0,
            "Surface coverage of Z [mol cm-2]": 0,
            "Electrode Coverage [mol cm-2]": Gamma,
            "Diffusion Coefficient of S [cm2 s-1]": DS_d,
            "Diffusion Coefficient of P [cm2 s-1]": DP_d,
            "Diffusion Coefficient of Y [cm2 s-1]": DY_d,
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
            "G": G,
            "G'": G_,
            "Reversible Potential 1 [V]": 0,
            "Electrosorption Rate [mol-1 cm3 s-1]": o,
            "Desorption Rate [s-1]": 0,
            "Catalytic Rate For [mol-1 cm3 s-1]": kf*1000,
            "Catalytic Rate Back [s-1]": kb*1000,
            "Symmetry factor [non-dim]": 0.5,
            
        }
        
        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps[0], x_steps[0], solver)
        #e, c = aa.AnalyticalModel(100)
        #setting solved answers to ones usable here
        current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        xss = x_steps[0]
        tss = t_steps[0]
        # while len(E_nd) == tss/2:
        #     xss = xss + 2
        #     tss = tss + 2
        #     cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, tss, xss, solver)
        #     current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current 
        E_d = E_nd / cmodel._E_0
        E_ds.append(E_d)
        I_ds.append(I_d)
        # Z_ds.append(Z_nd)
        np.savetxt(output+"test.txt", np.column_stack((E_d, I_d)))
        
        print("complete in time: " + str((time.time()-ti)/60) + " minutes") 
    
             

    #I_ana_d = I_ana_nd *cmodel._I_0
    #print(k0[0])
    # print(len(curr))
    
    #QUICK PLOTS, IMPROVE# 
    
    # Plot current
    plt.cla()
    for u in range(len(E_ds)):    
        plt.plot(E_ds[u], I_ds[u], label=str(k0[u]*K_0) + "- S - Red")  
        # plt.plot(E_ds[u], Z_ds[u], label=str(k0[u]*K_0) + "- S - Z")   
    # plt.plot(volt, np.array(-curr), color = 'orange', linestyle = 'dashdot', label = 'Digielch')
    # plt.plot(E_d, O_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, R_nd, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, current, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(T_nd, R_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
        np.savetxt(output+"test.txt", np.column_stack((E_ds[u], I_ds[u])))
    plt.xlabel("Eapp [V]")
    plt.ylabel("Current [A]")
    plt.legend(loc='upper left')
    plt.grid()
    plt.savefig(output+"CV_cat04_dim.png", dpi=1000)
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, O_nd, R_nd))))
    
    
    # #Plot concentration profiles
    # plt.cla()
    # # # plt.plot(volt[:398], np.array(curr[:398])/curr[0], color = 'orange', linestyle = 'dashdot', label = 'Digielch - S')
    # # # plt.plot(volt[399:], np.array(curr[399:])/curr[0], color = 'green', linestyle = 'dashdot', label = 'Digielch - P')
    # plt.plot(E_d[1:], O_nd[1:], color = 'Red', label = 'S', linestyle='dashed')
    # plt.plot(E_d[1:], R_nd[1:], color = 'Blue', label = 'P', linestyle='dashed')
    # plt.xlabel("Eapp [V]")
    # plt.ylabel("Surface Coverage [non-dim]")
    # plt.legend()
    # plt.grid()
    # plt.savefig(output+"CV_cat04_dim.png")
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, O_nd, R_nd))))
    
    return

if __name__ =='__main__':
    #Timer to measure performance
    ti = time.time()

    main()
    
