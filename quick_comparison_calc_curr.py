## Main function
import catalyticmodel04 as cm
import numpy as np
import time
import matplotlib.pylab as plt


def main(const_parameters,input_parameters):
    #list of options to pass into the model
    seioptions = ()
    files = [#['digielchcomp/surfaceonly/CV_k0_1.5.txt',2, 'digielchcomp/surfaceonly/SC_k0_1.5.txt'], 
             ['digielchcomp/surfaceonly/CV_k0_10.txt', 10, 'digielchcomp/surfaceonly/SC_k0_10.txt'],
             ['digielchcomp/surfaceonly/CV_k0_500.txt', 500, 'digielchcomp/surfaceonly/SC_k0_500.txt'],] 
             #['digielchcomp/surfaceonly/CV_k0_1000.txt', 1e5, 'digielchcomp/surfaceonly/SC_k0_1000.txt'],
             #['digielchcomp/surfaceonly/CV_k0_1e7.txt', 1e7, 'digielchcomp/surfaceonly/SC_k0_1e7.txt'], 
             #['digielchcomp/surfaceonly/CV_k0_1e9.txt', 1e9, 'digielchcomp/surfaceonly/SC_k0_1e9.txt'],]
    
    
    for i in files:
        input_parameters = {
            "Reversible Potential [V]": 0.0,
            "Redox Rate [s-1]": i[1],
            "Catalytic Rate For [cm2 mol-l s-1]": 0,
            "Catalytic Rate Back [cm2 mol-l s-1]": 0,
            "Symmetry factor [non-dim]": 0.5,
            "Capacitance [F]": 0, #1e-8,
            "Uncompensated Resistance [Ohm]": 0.0
        }
        
        voltage = []
        curr = []
        row = []
         
        f = open(i[0],'r')
        for row in f:
            row = row.split("\t")
            voltage.append(float(row[0]))
            curr.append(float(row[1]))

        #flipping current values from DigiElech
        curr = np.array(curr) * (-1)
        
        potential = []
        surf_cov = []
        f = open(i[2],'r')
        for row in f:
            row = row.split("\t")
            potential.append(float(row[0]))
            surf_cov.append(float(row[1])/100)
        
        
        #setting main model to reference CatalyticModel class
        cmodel = cm.CatalyticModel(const_parameters,seioptions)
        #setting solved answers to ones usable here
        E_nd, O_nd, R_nd, S_nd, P_nd, cat_conc, i_f, k0, T_nd = cmodel.simulate(input_parameters)
        ##redimensionalizing here for now. Messy to do in main, move later
        E_d = E_nd * cmodel._E_0
        
        current = []
        current.append(0)
        for v in range(1, len(O_nd)):
            current.append(O_nd[v]-O_nd[v-1])
            
        current = np.array(current)
            
        I_d = current * cmodel._I_0 * 2.52e2

        offset_E = []
        maxx = E_d[np.where(I_d==max(I_d))[0][0]]
        minn = E_d[np.where(I_d==min(I_d))[0][0]]
        for l in range(0, 10001):
            offset_E.append(E_d[l] + maxx)
        for l in range(10001, 20000):
            offset_E.append(E_d[l] + minn)
            
        k00 = str(i[1])
    
        #QUICK PLOTS# 
        plt.cla()
        plt.plot(E_d, I_d, color = "b", label = "PyBamm")
        plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        plt.title("Voltammogram k0: " + k00)
        plt.xlabel("Eapp [V]")
        plt.ylabel("current [A]")
        plt.grid()
        plt.legend()
        plt.savefig("output/CurrentvsEappdim_k0_"+k00+".png", dpi=600)
        np.savetxt("output/CurrentvsEappdim_k0_"+k00+".dat", np.transpose(np.vstack((E_d, I_d))))
        
        # plt.cla()
        # plt.plot(offset_E[:10001], (I_d[:10001]), color = "red", label = "Pybamm Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], (I_d[10001:]), color = "orange", label = "Pybamm Minimum peak (Reduction)")
        # plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.title("Offset Voltammogram k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig("output/CurrentvsOffsetEappdim_k0_"+k00+".png", dpi=600)
        
        # plt.cla()
        # plt.plot(offset_E[:10001], abs(I_d[:10001]), color = "red", label = "Maxmium peak (Oxidation)")
        # plt.plot(offset_E[10001:], abs(I_d[10001:]), color = "orange", label = "Minimum peak (Reduction)")
        # plt.title("k0: " + k00)
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.grid()
        # plt.legend()
        # plt.savefig("output/offset_Evscurrent_k0_"+k00+".png", dpi=600)
        # np.savetxt("output/offset_Evscurrent_k0_"+k00+".dat", np.transpose(np.vstack((E_d, I_d))))
        
        # plt.cla()
        # plt.plot(E_d, O_nd, color = "red", label="PyBamm - Ox")
        # plt.plot(E_d, R_nd, color = "orange", label="PyBamm - Red")
        # plt.plot(potential[:20000], surf_cov[:20000], color = "blue", label="DigiElch - Ox")
        # plt.plot(potential[20000:], surf_cov[20000:], color = "green", label="DigiElch - Red")
        # plt.title("Concentration Profile: " + k00)
        # plt.xlabel("Potential [V]")
        # plt.ylabel("Surface Coverage [non-dim]")
        # plt.legend()
        # plt.grid()
        # plt.savefig("output/SurfCoveragevsE_k0_"+k00+".png", dpi=600)
    
        # plt.cla()
        # plt.plot(E_nd, current)
        # plt.xlabel("Eapp [non-dim]")
        # plt.ylabel("current [non-dim]")
        # plt.grid()
        # plt.savefig("output/currentvsEapp_cat01.png", dpi=600)
        # np.savetxt("output/current_nondim_pybamm_kf_1.dat", np.transpose(np.vstack((E_nd, current))))
    
        # plt.cla()
        # plt.plot(T_nd, current)
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("current [non-dim]")
        # plt.grid()
        # plt.savefig("output/currentvstime_cat03.png", dpi=600)
        
        # plt.cla()
        # plt.plot(T_nd, E_nd)
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("Eapp [non-dim]")
        # plt.grid()
        # plt.savefig("output/Eappvstime_cat03.png", dpi=600)
        
        # plt.cla()
        # plt.plot(T_nd, cat_conc, label="Cat_conc")
        # plt.plot(T_nd, i_f, label="i_f")
        # plt.xlabel("time [non-dim]")
        # plt.ylabel("rates [non-dim]")
        # plt.legend()
        # plt.grid()
        # plt.savefig("output/ratesvstime_cat03.png", dpi=600)
        # np.savetxt("output/ratesvstime_cat01.dat", np.transpose(np.vstack((T_nd, O_nd, R_nd))))
        
        
        # plt.cla()
        # #plt.plot(voltage, curr, color = 'g', label = 'Digielch')
        # plt.plot(E_d, I_d, color = 'b', label = 'Pybamm')
        # #plt.title(title)
        # plt.legend()
        # plt.xlabel("Eapp [V]")
        # plt.ylabel("current [A]")
        # plt.savefig(i[0]+".png")
        
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
        "Redox Rate [s-1]": 1.5,
        "Catalytic Rate For [cm2 mol-l s-1]": 0,
        "Catalytic Rate Back [cm2 mol-l s-1]": 0,
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