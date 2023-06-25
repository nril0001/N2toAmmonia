## Main function
import surfacemodel01 as cm1
import surfacemodel02 as cm2
import re
import time
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
        
    
    files = [['DigiElech/2023-06-13 Surface only/RC ONLY/CV_k0_1e6_kf_0_kb_0_R_3e5_Cd_5.4e-7_Ds_0_Dp_0.txt']]
    
    #retrieve variables from file pathway
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
    
    curr = np.array(curr) * (-1)
        
    #constants that can vary, but generally won't change expt to expt
    const_parameters = {
        "Faraday Constant [C mol-1]": 96485.3328959,
        "Gas constant [J K-1 mol-1]": 8.314459848,
        "Far-field concentration of S(soln) [mol cm-3]": 1e-6,
        "Far-field concentration of P(soln) [mol cm-3]": 1e-6,
        "Diffusion Coefficient of S [cm2 s-1]": 1e-6,
        "Diffusion Coefficient of P [cm2 s-1]": 0,
        "Electrode Area [cm2]": 1,
        "Temperature [K]": 298.2,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": 0.0,
        "Scan Rate [V s-1]": 0.05,
        "Electrode Coverage [mol cm-2]": 1e-12,
    }
    
    k0 = files[0][5]
    Ru = 0
    Cdl = 0
    atol = 1e-10
    rtol = 1e-13
    steps = 80000

    #conditions that will change often over the course of testing
    input_parameters = {
        "Reversible Potential [V]": 0,
        "Redox Rate [s-1]": k0,
        "Catalytic Rate For [cm2 mol-l s-1]": 0,
        "Catalytic Rate Back [cm2 mol-l s-1]": 0,
        "Symmetry factor [non-dim]": 0.5,
        #28 Mar 2023: not fully implemented
        "Capacitance [F]": Cdl, #1e-8,
        "Uncompensated Resistance [Ohm]": Ru
    }
    
    
    #setting main model to reference CatalyticModel class
    # cmodel1 = cm1.CatalyticModel(const_parameters,seioptions)
    cmodel2 = cm2.CatalyticModel(const_parameters,seioptions, atol, rtol, steps)
    #e, c = aa.AnalyticalModel(100)
    #setting solved answers to ones usable here
    #current, E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
    # c1, E_nd1, O_nd1, R_nd1, times_nd1, Eeff1, E1, solution = cmodel1.simulate(input_parameters)
    c2, E_nd2, O_nd2, R_nd2, times_nd2, Eeff2 = cmodel2.simulate(input_parameters)

    # simulating analytical solution
    #I_ana_nd = amodel.simulate(E_nd)
    ##redimensionalizing here for now. Messy to do in main, move later
    # I_d1 = c1 * cmodel1._I_0
    # E_d1 = E_nd1 * cmodel1._E_0
    I_d2 = c2 * cmodel2._I_0
    E_d2 = E_nd2 * cmodel2._E_0

    #I_ana_d = I_ana_nd *cmodel._I_0
    #print(k0[0])
    
    #QUICK PLOTS, IMPROVE# 
    plt.cla()
    # plt.plot(E_d1[1:], I_d1[1:], color = "red", label = "Surface Model 1")
    plt.plot(E_d2[1:], I_d2[1:], color = "blue", label = "Surface Model 2")
    plt.plot(volt, curr, color = 'purple', label = 'Digielch')
    plt.title("K0: "+str(k0)+", Ru: "+str(Ru)+", Cdl: "+str(Cdl)+", atol: "+str(atol)+", rtol: "+str(rtol))
    plt.xlabel("Eapp [V]")
    plt.ylabel("current [A]")
    plt.grid()
    plt.legend()
    # plt.savefig(output+"CV_surf02_dim_ko_"+str(k0)+"_Ru_"+str(Ru)+"_Cdl_"+str(Cdl)+"_atol_"+str(atol)+"_rtol_"+str(rtol)+".png", dpi=600)
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d1, I_d1))))

    # plt.cla()
    # #plt.plot(E_d, I_d, color = "red", label = "current after model")
    # plt.plot(E_d[:10001], O_nd[:10001], color = "red", label = "O Maxmium peak (Oxidation)")
    # plt.plot(E_d[10001:], O_nd[10001:], color = "orange", label = "O Minimum peak (Reduction)")
    # plt.plot(E_d[:10001], R_nd[:10001], color = "purple", label = " R Maxmium peak (Oxidation)")
    # plt.plot(E_d[10001:], R_nd[10001:], color = "pink", label = "R Minimum peak (Reduction)")
    # plt.plot(volt, np.array(curr)/curr[0], color = 'g', label = 'Digielch')
    # #plt.plot(E_d, i, color = "blue", label = "current in model")
    # plt.xlabel("Eapp [non-dim]")
    # plt.ylabel("current [non-dim]")
    # plt.legend()
    #plt.savefig(output+"CV_cat04_nondim.png")
    #np.savetxt(output+"CV_cat04_nondim.dat", np.transpose(np.vstack((E_nd, current))))
    
    # plt.cla()
    # plt.plot(T_nd, E_nd)
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("Eapp [non-dim]")
    # plt.savefig(output+"Eappvstime_cat04_nondim.png")
    
    # plt.cla()
    # plt.plot(T_nd, cat_conc, label="Cat_conc")
    # plt.plot(T_nd, i_f, label="i_f")
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("rates [non-dim]")
    # plt.legend()
    # plt.savefig(output+"ratesvstime_cat04.png")
    # np.savetxt(output+"ratesvstime_cat04.dat", np.transpose(np.vstack((T_nd, O_nd, R_nd))))
    
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