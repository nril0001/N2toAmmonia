## Main function
import catalyticmodel04 as cm
import numpy as np
import time
import matplotlib.pylab as plt


def main(const_parameters,input_parameters):
    #list of options to pass into the model
    seioptions = ()
    #setting main model to reference CatalyticModel class
    cmodel = cm.CatalyticModel(const_parameters,seioptions)
    #setting solved answers to ones usable here
    current, E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
    # simulating analytical solution
    #I_ana_nd = aa.simulate(E_nd)
    ##redimensionalizing here for now. Messy to do in main, move later
    I_d = current * cmodel._I_0
    E_d = E_nd * cmodel._E_0

    #I_ana_d = I_ana_nd *cmodel._I_0
    print(k0[0])
    
    #QUICK PLOTS, IMPROVE# 
    plt.cla()
    plt.plot(E_d, I_d)
    plt.xlabel("Eapp [V]")
    plt.ylabel("current [A]")
    plt.savefig("output/CurrentvsEappdim_cat04_vs_temp.png")
    # np.savetxt("output/cu rent_dim_pybamm_kf_1.dat", np.transpose(np.vstack((E_d, I_d))))

    plt.cla()
    plt.plot(E_nd, current)
    plt.xlabel("Eapp [non-dim]")
    plt.ylabel("current [non-dim]")
    plt.savefig("output/currentvsEapp_cat04_vs_temp.png")
    # np.savetxt("output/current_nondim_pybamm_kf_1.dat", np.transpose(np.vstack((E_nd, current))))

    # plt.cla()
    # plt.plot(T_nd, current)
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("current [non-dim]")
    # plt.savefig("output/currentvstime_cat03.png")
    
    # plt.cla()
    # plt.plot(T_nd, E_nd)
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("Eapp [non-dim]")
    # plt.savefig("output/Eappvstime_cat03.png")
    
    # plt.cla()
    # plt.plot(T_nd, cat_conc, label="Cat_conc")
    # plt.plot(T_nd, i_f, label="i_f")
    # plt.xlabel("time [non-dim]")
    # plt.ylabel("rates [non-dim]")
    # plt.legend()
    # plt.savefig("output/ratesvstime_cat03.png")
    # np.savetxt("output/ratesvstime_cat01.dat", np.transpose(np.vstack((T_nd, O_nd, R_nd))))
    
    return

if __name__ =='__main__':
    #Timer to measure performance
    ti = time.time()

    #constants that can vary, but generally won't change expt to expt
    const_parameters = {
        "Faraday Constant [C mol-1]": 96485.3328959,
        "Gas constant [J K-1 mol-1]": 8.314459848,
        "Far-field concentration of S(soln) [mol cm-3]": 1e-2,
        "Far-field concentration of P(soln) [mol cm-3]": 0e-6,
        "Diffusion Coefficient of S [cm2 s-1]": 1e-5,
        "Diffusion Coefficient of P [cm2 s-1]": 1e-5,
        "Electrode Area [cm2]": 1,
        "Temperature [K]": 298.2,
        "Voltage start [V]": 0.5,
        "Voltage reverse [V]": -0.5,
        "Voltage amplitude [V]": 0.0,
        "Scan Rate [V s-1]": 0.05,
        "Electrode Coverage [mol cm-2]": 1e-12,
    }

    #conditions that will change often over the course of testing
    input_parameters = {
        "Reversible Potential [V]": 0.0,
        "Redox Rate [s-1]": 1e5,
        "Catalytic Rate For [cm2 mol-l s-1]": 1e4,
        "Catalytic Rate Back [cm2 mol-l s-1]": 1e-3,
        "Symmetry factor [non-dim]": 0.5,
        #28 Mar 2023: not fully implemented
        "Capacitance [F]": 0, #1e-8,
        "Uncompensated Resistance [Ohm]": 0.0
    }

    main(const_parameters,input_parameters)
    

    #literally just to test that that main is working properly (delete later)
    print("complete in time: " + str((time.time()-ti)/60) + " minutes")

    #TODO: add multiple fxns for running models at the same time
    #TODO: improve plotting for pybamm based models: investigate using pybamm's plotting features??
    #TODO: LOW PRIORITY make compatible with interactive plottersz