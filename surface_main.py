## Main function
import surfacemodel01 as cm
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
    
    k0 = 267

    #conditions that will change often over the course of testing
    input_parameters = {
        "Reversible Potential [V]": 0,
        "Redox Rate [s-1]": k0,
        "Catalytic Rate For [cm2 mol-l s-1]": 0,
        "Catalytic Rate Back [cm2 mol-l s-1]": 0,
        "Symmetry factor [non-dim]": 0.5,
        #28 Mar 2023: not fully implemented
        "Capacitance [F]": 0, #1e-8,
        "Uncompensated Resistance [Ohm]": 0
    }
    
    files = [['C:/Users/natha/Desktop/Code/DigiElech/2023-06-06 Solution only/SC/SC_k0_1e-3_Ds_1e-5_Dp_1e-5.txt', 1e-3, 1e-5]]
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
            if row[0] == '':
                continue
            else:
                volt.append(float(row[0]))
                curr.append(float(row[1]))
    
    curr = np.array(curr)
    
    #setting main model to reference CatalyticModel class
    cmodel = cm.CatalyticModel(const_parameters,seioptions)
    #e, c = aa.AnalyticalModel(100)
    #setting solved answers to ones usable here
    #current, E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
    c, E_nd, O_nd, R_nd, times_nd, Eeff, E = cmodel.simulate(input_parameters)
    # simulating analytical solution
    #I_ana_nd = amodel.simulate(E_nd)
    ##redimensionalizing here for now. Messy to do in main, move later
    I_d = c * cmodel._I_0
    E_d = E_nd * cmodel._E_0

    #I_ana_d = I_ana_nd *cmodel._I_0
    #print(k0[0])
    
    #QUICK PLOTS, IMPROVE# 
    plt.cla()
    plt.plot(E_d, I_d)
    plt.title("K0: "+str(k0))
    plt.xlabel("Eapp [V]")
    plt.ylabel("current [A]")
    plt.grid()
    #plt.savefig(output+"CV_cat04_dim_ko_"+str(k0)+".png")
    # np.savetxt(output+"current_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, I_d))))

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
    # plt.plot(T_nd, current)
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("current [non-dim]")
    # plt.savefig(output+"Curr_vs_Time_cat04_nondim.png")
    
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