## Main function - confirms that model 7 acts in line with paper
import catalyticmodel07_b as cm
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
    # files = [['DigiElech/2023-10-02 adsorption/CV_R_130_v_2e-2_k0_1e-3_G_0_Gamma_1.txt']]

    # retrieve variables from file pathway
    # for t in files:
    #     variables = re.findall(
    #         "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", t[0])
    #     for l in range(len(variables)):
    #         variables[l] = float(variables[l])
    #         t.append(variables[l])
    
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
    area = 1
    radius = np.sqrt(area/np.pi)
    srate = [1]
    k0 = 1e4
    kads = [2e-4, 8e-4, 3.3e-3, 1.3e-2, 5.2e-2, 2.1e-1, 6.7]
    # kads = [1]
    kdes = 1
    Ru = 0
    Cdl = 0
    atol = 1e-14
    rtol = 1e-14
    t_steps = [2**(13)]
    x_steps = [750]
    solver = "Casadi"
    # solver = "Scikits"
    
    D_thickness = 1
    DS_d = 1
    CS_d = 1
    Gamma = 1
    F = 96485.3328959
    R = 8.314459848
    T = 298.2
    
    E_ds = []
    I_ds = []
    Z_ds = []
    
    for o in kads:   
        print("k0 = " + str(k0) + ", " + str(o*1000) + " mV/s")
        #constants that can vary, but generally won't change expt to expt
        const_parameters = {
            "Faraday Constant [C mol-1]": F,
            "Gas constant [J K-1 mol-1]": R,
            "Temperature [K]": T,
            "Far-field concentration of S(soln) [mol cm-3]": CS_d,
            "Far-field concentration of P(soln) [mol cm-3]": 0,
            "Surface coverage of S [mol cm-2]": 0,
            "Surface coverage of P [mol cm-2]": 0,
            "Electrode Coverage [mol cm-2]": Gamma,
            "Diffusion Coefficient of S [cm2 s-1]": DS_d,
            "Diffusion Layer Thickness [cm]": D_thickness,
            "Electrode Area [cm2]": 1,
            "Electrode Radius [cm]": radius,
            "Voltage start [V]": 30,
            "Voltage reverse [V]": -30,
            "Voltage amplitude [V]": 0.0,
            "Scan Rate [V s-1]": srate[0],
            "Uncompensated Resistance [Ohm]": Ru,
            "Capacitance [F]": Cdl, #1e-8,
        }
        
        #conditions that will change often over the course of testing
        input_parameters = {
            "Reversible Potential 1 [V]": 0,
            "Redox Rate (ads) [s-1]": k0,
            "Adsorption Rate [mol-1 cm3 s-1]": o,
            "Desorption Rate [s-1]": kdes,
            "Symmetry factor [non-dim]": 0.5,
            
        }
        
        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, t_steps[0], x_steps[0], solver)
        #setting solved answers to ones usable here
        current, E_nd, O_nd, R_nd, T_nd, X = cmodel.simulate(input_parameters)
        xss = x_steps[0]
        tss = t_steps[0]
        while len(E_nd) == tss/2:
            xss = xss + 2
            tss = tss + 2
            cmodel = cm.CatalyticModel(const_parameters,seioptions, atol, rtol, tss, xss, solver)
            current, E_nd, O_nd, R_nd, T_nd, X = cmodel.simulate(input_parameters)
        ##redimensionalizing here for now. Messy to do in main, move later
        I_d = current 
        E_d = E_nd 
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
        plt.plot(E_ds[u], (I_ds[u]), label="v = " + str(kads[u]*1000) + " mV/s")
    # plt.plot(volt, np.array(-curr), linestyle = 'dashdot', label = 'Digielch - ks = 1e-3, Gmax = 1e-9')
    # plt.plot(E_d, O_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, R_nd, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(E_d, current, color = 'Blue', label = 'Pybamm', linestyle='dashed')
    # plt.plot(T_nd, R_nd, color = 'Red', label = 'Pybamm', linestyle='dashed')
        # np.savetxt(output+"test.txt", np.column_stack((E_ds[u], I_ds[u])))
    plt.xlabel("Eapp (V)")
    # plt.ylabel("Current (A)")
    plt.ylabel("j / mA cm-2")
    plt.title("k0 = " + str(k0) + " s-1")
    plt.legend(loc='upper left')
    plt.grid()
    plt.savefig(output+"CV_cat07_multi_srate_casadi1.png", dpi=600)
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
    
