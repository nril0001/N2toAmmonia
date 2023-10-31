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
    
    files = [['DigiElech/2023-10-02 adsorption/CV_R_130_v_2e-2_k0_1e-3_G_0.txt']]
    
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
    #edit to control parameters
    area = 0.05
    radius = np.sqrt(area/np.pi)
    srate = [0.02, 0.1, 0.2, 0.4, 1]
    k0 = 1e-6
    Ru = 130
    Cdl = 0
    atol = 1e-7
    rtol = 1e-7
    t_steps = [2**(12)]
    x_steps = [750]
    solver = "Casadi"
    # solver = "Scikits"
    
    D_thickness = 3.882e-3
    DS_d = 2.27e-6
    DP_d = 2.27e-6
    CS_d = 0.002
    Gamma = 2.2e-9
    F = 96485.3328959
    R = 8.314459848
    T = 298.2
    G = 0
    G_ = 0
    
    E_ds = []
    I_ds = []
    Z_ds = []
    
    for o in srate:   
        print("k0 = " + str(k0) + ", " + str(o*1000) + " mV/s")
        #constants that can vary, but generally won't change expt to expt
        const_parameters = {
            "Faraday Constant [C mol-1]": F,
            "Gas constant [J K-1 mol-1]": R,
            "Temperature [K]": T,
            "Standard Unity Concentration [mol cm-3]": 0.001,
            "Far-field concentration of S(soln) [mol cm-3]": CS_d,
            "Far-field concentration of P(soln) [mol cm-3]": 0,
            "Surface coverage of P [mol cm-2]": 0,
            "Electrode Coverage [mol cm-2]": Gamma,
            "Diffusion Coefficient of S [cm2 s-1]": DS_d,
            "Diffusion Coefficient of P [cm2 s-1]": DP_d,
            "Diffusion Layer Thickness [cm]": D_thickness,
            "Electrode Area [cm2]": 0.05,
            "Electrode Radius [cm]": radius,
            "Voltage start [V]": 2.0,
            "Voltage reverse [V]": -0.55,
            "Voltage amplitude [V]": 0.0,
            "Scan Rate [V s-1]": o,
            "Uncompensated Resistance [Ohm]": Ru,
            "Capacitance [F]": Cdl, #1e-8,
        }
        
        #conditions that will change often over the course of testing
        input_parameters = {
            "G": G,
            "G'": G_,
            "Reversible Potential 1 [V]": 0,
            "Electrosorption Rate [mol-1 cm3 s-1]": k0,
            "Symmetry factor [non-dim]": 0.63,
            
        }
        
        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps[0], x_steps[0], solver)
        #setting solved answers to ones usable here
        current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        xss = x_steps[0]
        tss = t_steps[0]
        while len(E_nd) == tss/2:
            xss = xss + 2
            cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, tss, xss, solver)
            current, E_nd, O_nd, R_nd, T_nd = cmodel.simulate(input_parameters)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current / cmodel._I_0
        E_d = E_nd / cmodel._E_0
        E_ds.append(E_d)
        I_ds.append(I_d)
        # Z_ds.append(Z_nd)
        # np.savetxt(output+"test.txt", np.column_stack((E_d, I_d)))
        
        print("complete in time: " + str((time.time()-ti)/60) + " minutes") 
    
             

    #I_ana_d = I_ana_nd *cmodel._I_0
    #print(k0[0])
    # print(len(curr))
    
    #QUICK PLOTS, IMPROVE# 
    
    # Plot current
    plt.cla()
    for u in range(len(E_ds)):    
        plt.plot(E_ds[u], (I_ds[u] / area * 1000), label="v = " + str(srate[u]*1000) + " mV/s")
    # plt.plot(volt, np.array(-curr), linestyle = 'dashdot', label = 'Digielch - ks = 1e-3, Gmax = 1e-9')
    # plt.plot(E_d, O_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, R_nd, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, current, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(T_nd, R_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
        # np.savetxt(output+"test.txt", np.column_stack((E_ds[u], I_ds[u])))
    plt.xlabel("Eapp (V)")
    # plt.ylabel("Current (A)")
    plt.ylabel("j / mA cm-2")
    plt.title("k0 = " + str(k0) + " cm3 mol-1 s-1")
    plt.legend(loc='upper right')
    plt.grid()
    path = "CV_cat06_casadi1_k0_"
    ext = ".png"
    ext1 = ".dat"
    count = 1
    while os.path.exists(output+path+str(k0)+"_Ru_"+str(Ru)+ext):
        path = path+"_"+str(count)
        count += 1
    plt.savefig(output+path+str(k0)+"_Ru_"+str(Ru)+ext, dpi=600)
    np.savetxt(output+path+str(k0)+"_Ru_"+str(Ru)+ext1, np.transpose(np.vstack((E_ds[-1], I_ds))))
    
    
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
    
